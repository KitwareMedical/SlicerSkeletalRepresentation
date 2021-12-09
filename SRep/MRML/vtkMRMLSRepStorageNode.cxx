#include "vtkMRMLSRepStorageNode.h"
#include "vtkMRMLEllipticalSRepNode.h"
#include "vtkMRMLRectangularGridSRepNode.h"
#include "vtkMRMLMessageCollection.h"
#include "vtkMRMLScene.h"
#include "vtkStringArray.h"

#include "srep/Util.h"

// Relax JSON standard and allow reading/writing of nan and inf
// values. Such values should not normally occur, but if they do then
// it is easier to troubleshoot problems if numerical values are the
// same in memory and files.
// kWriteNanAndInfFlag = 2,        //!< Allow writing of Infinity, -Infinity and NaN.
#define RAPIDJSON_WRITE_DEFAULT_FLAGS 2
// kParseNanAndInfFlag = 256,      //!< Allow parsing NaN, Inf, Infinity, -Inf and -Infinity as doubles.
#define RAPIDJSON_PARSE_DEFAULT_FLAGS 256

#include "rapidjson/document.h"     // rapidjson's DOM-style API
#include "rapidjson/prettywriter.h" // for stringify JSON
#include "rapidjson/filereadstream.h"
#include "rapidjson/filewritestream.h"

using srep::util::finally;

namespace {

namespace keys {
  const char * const RectangularGridSRep = "RectangularGridSRep";
  const char * const EllipticalSRep = "EllipticalSRep";
  const char * const Rows = "Rows";
  const char * const Cols = "Cols";
  const char * const CrestPoints = "CrestPoints";
  const char * const Steps = "Steps";
  const char * const Skeleton = "Skeleton";
  const char * const UpSpoke = "UpSpoke";
  const char * const DownSpoke = "DownSpoke";
  const char * const CrestSpoke = "CrestSpoke";
  const char * const Direction = "Direction";
  const char * const SkeletalPoint = "SkeletalPoint";
  const char * const Value = "Value";

  const char * const Display = "Display";
  const char * const Visibility = "Visibility";
  const char * const PiecewiseVisibility = "PiecewiseVisibility";
  const char * const Opacity = "Opacity";
  const char * const Colors = "Colors";
  const char * const SkeletalSheet = "SkeletalSheet";
  const char * const CrestCurve = "CrestCurve";
  const char * const SkeletonToCrestConnection = "SkeletonToCrestConnection";
  const char * const RelativeThickness = "RelativeThickness";
  const char * const AbsoluteThickness = "AbsoluteThickness";
  const char * const UseAbsoluteThickness = "UseAbsoluteThickness";

  const char * const CoordinateSystem = "CoordinateSystem";
}

void writeCoordinateSystem(rapidjson::PrettyWriter<rapidjson::FileWriteStream>& writer, int storageCoord) {
  if (storageCoord == vtkMRMLStorageNode::CoordinateSystemLPS) {
    writer.String("LPS");
  } else if (storageCoord == vtkMRMLStorageNode::CoordinateSystemRAS) {
    writer.String("RAS");
  } else {
    throw std::invalid_argument("Unknown storage node coordinate system type: " + std::to_string(storageCoord));
  }
}

int readCoordinateSystem(rapidjson::Value& json) {
  if (!json.IsString()) {
    throw std::invalid_argument("Expect string for coordinate system");
  }

  const std::string value = json.GetString();
  if (value == "LPS") {
    return vtkMRMLStorageNode::CoordinateSystemLPS;
  } else if (value == "RAS") {
    return vtkMRMLStorageNode::CoordinateSystemRAS;
  } else {
    throw std::invalid_argument("Unknown srep coordinate system type: " + value);
  }
}

std::array<double, 3> FromRASToCoord(const std::array<double, 3>& arr, int storageCoord) {
  if (storageCoord == vtkMRMLStorageNode::CoordinateSystemLPS) {
    return std::array<double, 3>{-arr[0], -arr[1], arr[2]};
  } else if (storageCoord == vtkMRMLStorageNode::CoordinateSystemRAS) {
    return arr;
  } else {
    throw std::invalid_argument("Unknown coordinate system type: " + std::to_string(storageCoord));
  }
}

std::array<double, 3> FromCoordToRAS(const std::array<double, 3>& arr, int storageCoord) {
  // same transformation both ways, but two functions helps keep things straight.
  return FromRASToCoord(arr, storageCoord);
}

constexpr size_t BufferSize = 65535;

rapidjson::Value::MemberIterator SafeFindMember(rapidjson::Value& json, const char* name) {
  auto iter = json.FindMember(name);
  if (iter == json.MemberEnd()) {
    throw std::invalid_argument(std::string("Error finding json member '") + name + "'");
  }
  return iter;
}

double readDouble(rapidjson::Value& json) {
  if (!json.IsDouble()) {
    throw std::invalid_argument("Expected a JSON double.");
  }
  return json.GetDouble();
}

bool readBool(rapidjson::Value& json) {
  if (!json.IsBool()) {
    throw std::invalid_argument("Expected a JSON bool.");
  }
  return json.GetBool();
}

int readInt(rapidjson::Value& json) {
  if (!json.IsInt()) {
    throw std::invalid_argument("Expected a JSON int.");
  }
  return json.GetInt();
}

template<size_t N>
void writeSingleLineArray(rapidjson::PrettyWriter<rapidjson::FileWriteStream>& writer, const std::array<double, N>& arr) {
  writer.SetFormatOptions(rapidjson::kFormatSingleLineArray);
  writer.StartArray();
  for (const auto d : arr) {
    writer.Double(d);
  }
  writer.EndArray();
  writer.SetFormatOptions(rapidjson::kFormatDefault);
}

void write(rapidjson::PrettyWriter<rapidjson::FileWriteStream>& writer, const srep::Point3d& point, int coordinateSystem) {
  writer.StartObject();
  writer.Key(keys::CoordinateSystem); writeCoordinateSystem(writer, coordinateSystem);
  writer.Key(keys::Value);
  writeSingleLineArray(writer, FromRASToCoord(point.AsArray(), coordinateSystem));
  writer.EndObject();
}

std::array<double, 3> read3dArray(rapidjson::Value& json) {
  if (!json.IsArray()) {
    throw std::invalid_argument("Attempting to read an array that is not a json array");
  }
  auto jsonArray = json.GetArray();

  if (jsonArray.Size() != 3) {
    throw std::invalid_argument("Attempting to read a 3D array that doesn't have 3 dimensions");
  }

  std::array<double, 3> array;
  std::transform(jsonArray.Begin(), jsonArray.End(), array.begin(), readDouble);
  return array;
}

srep::Point3d readPoint3d(rapidjson::Value& json) {
  if (!json.IsObject()) {
    throw std::invalid_argument("Expected a json object for srep::Point3d");
  }
  return srep::Point3d(FromCoordToRAS(
    read3dArray(SafeFindMember(json, keys::Value)->value),
    readCoordinateSystem(SafeFindMember(json, keys::CoordinateSystem)->value)));
}

void write(rapidjson::PrettyWriter<rapidjson::FileWriteStream>& writer, const srep::Vector3d& vector, int coordinateSystem) {
  writer.StartObject();
  writer.Key(keys::CoordinateSystem); writeCoordinateSystem(writer, coordinateSystem);
  writer.Key(keys::Value);
  writeSingleLineArray(writer, FromRASToCoord(vector.AsArray(), coordinateSystem));
  writer.EndObject();
}

srep::Vector3d readVector3d(rapidjson::Value& json) {
  if (!json.IsObject()) {
    throw std::invalid_argument("Expected a json object for srep::Point3d");
  }
  return srep::Vector3d(FromCoordToRAS(
    read3dArray(SafeFindMember(json, keys::Value)->value),
    readCoordinateSystem(SafeFindMember(json, keys::CoordinateSystem)->value)));
}

void write(rapidjson::PrettyWriter<rapidjson::FileWriteStream>& writer, const srep::Spoke& spoke, int coordinateSystem) {
  writer.StartObject();
  writer.Key(keys::SkeletalPoint);
  write(writer, spoke.GetSkeletalPoint(), coordinateSystem);
  writer.Key(keys::Direction);
  write(writer, spoke.GetDirection(), coordinateSystem);
  writer.EndObject();
}

srep::Spoke readSpoke(rapidjson::Value& json) {
  auto skeletalIter = SafeFindMember(json, keys::SkeletalPoint);
  auto directionIter = SafeFindMember(json, keys::Direction);
  return srep::Spoke(readPoint3d(skeletalIter->value), readVector3d(directionIter->value));
}

// raw write, no concept of what the rows and cols mean
void write(rapidjson::PrettyWriter<rapidjson::FileWriteStream>& writer, const std::vector<std::vector<srep::SkeletalPoint>>& grid, int coordinateSystem) {
  writer.StartArray();
  for (size_t i  = 0; i < grid.size(); ++i) {
    writer.StartArray();
    for (size_t j = 0; j < grid[i].size(); ++j) {
      const auto& skeletalPoint = grid[i][j];

      writer.StartObject();
      writer.Key(keys::UpSpoke);
      write(writer, skeletalPoint.GetUpSpoke(), coordinateSystem);
      writer.Key(keys::DownSpoke);
      write(writer, skeletalPoint.GetDownSpoke(), coordinateSystem);
      if (skeletalPoint.IsCrest()) {
        writer.Key(keys::CrestSpoke);
        write(writer, skeletalPoint.GetCrestSpoke(), coordinateSystem);
      }
      writer.EndObject();
    }
    writer.EndArray();
  }
  writer.EndArray();
}

// raw read, no concept of what the rows and cols mean
std::vector<std::vector<srep::SkeletalPoint>> read2DSkeletalPointVec(rapidjson::Value& json) {
  if (!json.IsArray()) {
    throw std::invalid_argument("Expected a JSON array.");
  }

  std::vector<std::vector<srep::SkeletalPoint>> grid;
  for (auto& row : json.GetArray()) {
    grid.push_back(std::vector<srep::SkeletalPoint>{});
    if (!row.IsArray()) {
      throw std::runtime_error("Error parsing in read2DSkeletalPointVec. Row is not array.");
    }
    for (auto& object : row.GetArray()) {
      auto upIter = SafeFindMember(object, keys::UpSpoke);
      auto downIter = SafeFindMember(object, keys::DownSpoke);

      auto upSpoke = readSpoke(upIter->value);
      auto downSpoke = readSpoke(downIter->value);

      auto crestIter = object.FindMember(keys::CrestSpoke);
      if (crestIter != object.MemberEnd()) {
        auto crestSpoke = readSpoke(crestIter->value);
        grid.back().push_back(srep::SkeletalPoint(upSpoke, downSpoke, crestSpoke));
      } else {
        grid.back().push_back(srep::SkeletalPoint(upSpoke, downSpoke));
      }
    }
  }

  return grid;
}

void write(rapidjson::PrettyWriter<rapidjson::FileWriteStream>& writer, vtkMRMLRectangularGridSRepNode& mrmlSRep, int coordinateSystem)
{
  const srep::RectangularGridSRep* srep = mrmlSRep.GetRectangularGridSRep();
  //writing out rows and cols more for people who look at the JSON
  const auto& rows = srep ? srep->GetNumRows() : 0;
  const auto& cols = srep ? srep->GetNumCols() : 0;
  const auto& grid = srep ? srep->GetSkeletalPoints() : srep::RectangularGridSRep::SkeletalGrid{};

  writer.Key(keys::RectangularGridSRep);
  writer.StartObject();
  {
    writer.Key(keys::Rows); writer.Uint(rows);
    writer.Key(keys::Cols); writer.Uint(cols);

    writer.Key(keys::Skeleton);
    write(writer, grid, coordinateSystem);
  }
  writer.EndObject();
}

void read(rapidjson::Value& json, vtkMRMLRectangularGridSRepNode* rectangularGridSRep) {
  if (!rectangularGridSRep) {
    throw std::invalid_argument("Node is not a vtkMRMLRectangularGridSRepNode");
  }

  auto skeletonIter = SafeFindMember(json, keys::Skeleton);
  auto grid = read2DSkeletalPointVec(skeletonIter->value);

  // RectangularGridSRep constructor will throw if there are bad things like crest spokes not on the crest
  rectangularGridSRep->SetRectangularGridSRep(std::unique_ptr<srep::RectangularGridSRep>(new srep::RectangularGridSRep(std::move(grid))));
}

void write(rapidjson::PrettyWriter<rapidjson::FileWriteStream>& writer, vtkMRMLEllipticalSRepNode& mrmlSRep, int coordinateSystem) {
  const auto srep = mrmlSRep.GetEllipticalSRep();
  //writing out numFoldPoints and numSteps more for people who look at the JSON
  const auto numFoldPoints = srep ? srep->GetSkeleton().size() : 0;
  const auto numSteps = srep ? (srep->GetSkeleton().empty() ? 0 : srep->GetSkeleton()[0].size() - 1) : 0; // -1 because "step 0" is the spine and doesn't count
  const auto grid = srep ? srep->GetSkeleton() : srep::EllipticalSRep::UnrolledEllipticalGrid{};

  writer.Key(keys::EllipticalSRep);
  writer.StartObject();
  {
    writer.Key(keys::CrestPoints); writer.Uint(numFoldPoints);
    writer.Key(keys::Steps); writer.Uint(numSteps);
    writer.Key(keys::Skeleton);
    write(writer, grid, coordinateSystem);
  }
  writer.EndObject();
}

void read(rapidjson::Value& json, vtkMRMLEllipticalSRepNode* ellipticalSRep) {
  if (!ellipticalSRep) {
    throw std::invalid_argument("Node is not a vtkMRMLEllipticalSRepNode");
  }

  auto skeletonIter = SafeFindMember(json, keys::Skeleton);
  auto grid = read2DSkeletalPointVec(skeletonIter->value);

  // RectangularGridSRep constructor will throw if there are bad things like crest spokes not on the crest
  ellipticalSRep->SetEllipticalSRep(std::unique_ptr<srep::EllipticalSRep>(new srep::EllipticalSRep(grid)));
}

void write(rapidjson::PrettyWriter<rapidjson::FileWriteStream>& writer, const vtkColor3ub& color) {
  writer.SetFormatOptions(rapidjson::kFormatSingleLineArray);
  writer.StartArray();
  writer.Int(color[0]);
  writer.Int(color[1]);
  writer.Int(color[2]);
  writer.EndArray();
  writer.SetFormatOptions(rapidjson::kFormatDefault);
}

vtkColor3ub readVTKColor3ub(rapidjson::Value& json) {
  if (!json.IsArray()) {
    throw std::invalid_argument("Attempting to read an array that is not a json array");
  }
  auto jsonArray = json.GetArray();

  if (jsonArray.Size() != 3) {
    throw std::invalid_argument("Attempting to read a 3D array that doesn't have 3 dimensions");
  }

  std::array<int, 3> array;
  std::transform(jsonArray.Begin(), jsonArray.End(), array.begin(), readInt);
  return vtkColor3ub(array[0], array[1], array[2]);
}

void writeDisplayNodeColors(rapidjson::PrettyWriter<rapidjson::FileWriteStream>& writer, vtkMRMLSRepDisplayNode& displayNode) {
  writer.StartObject();
  writer.Key(keys::UpSpoke);
  write(writer, displayNode.GetUpSpokeColor());
  writer.Key(keys::DownSpoke);
  write(writer, displayNode.GetDownSpokeColor());
  writer.Key(keys::CrestSpoke);
  write(writer, displayNode.GetCrestSpokeColor());
  writer.Key(keys::SkeletalSheet);
  write(writer, displayNode.GetSkeletalSheetColor());
  writer.Key(keys::CrestCurve);
  write(writer, displayNode.GetCrestCurveColor());
  writer.Key(keys::SkeletonToCrestConnection);
  write(writer, displayNode.GetSkeletonToCrestConnectionColor());
  writer.EndObject();
}

void readDisplayNodeColors(rapidjson::Value& json, vtkMRMLSRepDisplayNode& displayNode) {
  if (!json.IsObject()) {
    throw std::invalid_argument("Attempting to read vtkMRMLSRepDisplayNode colors but json is not an object");
  }
  using SetColorFunc = void (vtkMRMLSRepDisplayNode::*)(const vtkColor3ub&);
  const auto readAndSetIfExists = [&displayNode, &json](const char* key, SetColorFunc setColor) {
    auto iter = json.FindMember(key);
    if (iter != json.MemberEnd()) {
      (displayNode.*setColor)(readVTKColor3ub(iter->value));
    }
  };

  readAndSetIfExists(keys::UpSpoke, &vtkMRMLSRepDisplayNode::SetUpSpokeColor);
  readAndSetIfExists(keys::DownSpoke, &vtkMRMLSRepDisplayNode::SetDownSpokeColor);
  readAndSetIfExists(keys::CrestSpoke, &vtkMRMLSRepDisplayNode::SetCrestSpokeColor);
  readAndSetIfExists(keys::CrestCurve, &vtkMRMLSRepDisplayNode::SetCrestCurveColor);
  readAndSetIfExists(keys::SkeletalSheet, &vtkMRMLSRepDisplayNode::SetSkeletalSheetColor);
  readAndSetIfExists(keys::SkeletonToCrestConnection, &vtkMRMLSRepDisplayNode::SetSkeletonToCrestConnectionColor);
}

void writeDisplayNodePiecewiseVisibilities(rapidjson::PrettyWriter<rapidjson::FileWriteStream>& writer, vtkMRMLSRepDisplayNode& displayNode) {
  writer.StartObject();
  writer.Key(keys::UpSpoke);
  writer.Bool(displayNode.GetUpSpokeVisibility());
  writer.Key(keys::DownSpoke);
  writer.Bool(displayNode.GetDownSpokeVisibility());
  writer.Key(keys::CrestSpoke);
  writer.Bool(displayNode.GetCrestSpokeVisibility());
  writer.Key(keys::SkeletalSheet);
  writer.Bool(displayNode.GetSkeletalSheetVisibility());
  writer.Key(keys::CrestCurve);
  writer.Bool(displayNode.GetCrestCurveVisibility());
  writer.Key(keys::SkeletonToCrestConnection);
  writer.Bool(displayNode.GetSkeletonToCrestConnectionVisibility());
  writer.EndObject();
}

void readDisplayNodePiecewiseVisibilities(rapidjson::Value& json, vtkMRMLSRepDisplayNode& displayNode) {
  if (!json.IsObject()) {
    throw std::invalid_argument("Attempting to read vtkMRMLSRepDisplayNode colors but json is not an object");
  }
  using SetVisibilityFunc = void (vtkMRMLSRepDisplayNode::*)(bool);
  const auto readAndSetIfExists = [&displayNode, &json](const char* key, SetVisibilityFunc setVisibility) {
    auto iter = json.FindMember(key);
    if (iter != json.MemberEnd()) {
      (displayNode.*setVisibility)(readBool(iter->value));
    }
  };

  readAndSetIfExists(keys::UpSpoke, &vtkMRMLSRepDisplayNode::SetUpSpokeVisibility);
  readAndSetIfExists(keys::DownSpoke, &vtkMRMLSRepDisplayNode::SetDownSpokeVisibility);
  readAndSetIfExists(keys::CrestSpoke, &vtkMRMLSRepDisplayNode::SetCrestSpokeVisibility);
  readAndSetIfExists(keys::CrestCurve, &vtkMRMLSRepDisplayNode::SetCrestCurveVisibility);
  readAndSetIfExists(keys::SkeletalSheet, &vtkMRMLSRepDisplayNode::SetSkeletalSheetVisibility);
  readAndSetIfExists(keys::SkeletonToCrestConnection, &vtkMRMLSRepDisplayNode::SetSkeletonToCrestConnectionVisibility);
}

void write(rapidjson::PrettyWriter<rapidjson::FileWriteStream>& writer, vtkMRMLSRepDisplayNode& displayNode) {
  writer.Key(keys::Display);
  writer.StartObject();
  writer.Key(keys::Visibility);
  writer.Bool(static_cast<bool>(displayNode.GetVisibility()));
  writer.Key(keys::Opacity);
  writer.Double(displayNode.GetOpacity());
  writer.Key(keys::RelativeThickness);
  writer.Double(displayNode.GetRelativeThickness());
  writer.Key(keys::AbsoluteThickness);
  writer.Double(displayNode.GetAbsoluteThickness());
  writer.Key(keys::UseAbsoluteThickness);
  writer.Bool(displayNode.GetUseAbsoluteThickness());
  writer.Key(keys::Colors);
  writeDisplayNodeColors(writer, displayNode);
  writer.Key(keys::PiecewiseVisibility);
  writeDisplayNodePiecewiseVisibilities(writer, displayNode);
  writer.EndObject();
}

void read(rapidjson::Value& json, vtkMRMLSRepDisplayNode& displayNode) {
  if (!json.IsObject()) {
    throw std::invalid_argument("Attempting to read vtkMRMLSRepDisplayNode but json is not an object");
  }

  // we will assume if an item is not there, we just don't set it.
  auto visibilityIter = json.FindMember(keys::Visibility);
  if (visibilityIter != json.MemberEnd()) {
    displayNode.SetVisibility(readBool(visibilityIter->value));
  }
  auto opacityIter = json.FindMember(keys::Opacity);
  if (opacityIter != json.MemberEnd()) {
    displayNode.SetOpacity(readDouble(opacityIter->value));
  }
  auto colorIter = json.FindMember(keys::Colors);
  if (colorIter != json.MemberEnd()) {
    readDisplayNodeColors(colorIter->value, displayNode);
  }
  auto piecewiseVisibilitiesIter = json.FindMember(keys::PiecewiseVisibility);
  if (piecewiseVisibilitiesIter != json.MemberEnd()) {
    readDisplayNodePiecewiseVisibilities(piecewiseVisibilitiesIter->value, displayNode);
  }
  auto relativeThicknessIter = json.FindMember(keys::RelativeThickness);
  if (relativeThicknessIter != json.MemberEnd()) {
    displayNode.SetRelativeThickness(readDouble(relativeThicknessIter->value));
  }
  auto absoluteThicknessIter = json.FindMember(keys::AbsoluteThickness);
  if (absoluteThicknessIter != json.MemberEnd()) {
    displayNode.SetAbsoluteThickness(readDouble(absoluteThicknessIter->value));
  }
  auto useAbsoluteThicknessIter = json.FindMember(keys::UseAbsoluteThickness);
  if (useAbsoluteThicknessIter != json.MemberEnd()) {
    displayNode.SetUseAbsoluteThickness(readBool(useAbsoluteThicknessIter->value));
  }
}

// must not return nullptr. Throws on error.
std::unique_ptr<rapidjson::Document> CreateJsonDocumentFromFile(const char* filePath) {
  FILE* fp = fopen(filePath, "r");
  if (!fp) {
    throw std::runtime_error("Error opening file");
  }
  const auto closeFp = finally([fp](){
    fclose(fp);
  });

  std::array<char, BufferSize> buffer;
  rapidjson::FileReadStream fs(fp, buffer.data(), buffer.size());
  std::unique_ptr<rapidjson::Document> jsonRoot(new rapidjson::Document);

  if (jsonRoot->ParseStream(fs).HasParseError()) {
    throw std::runtime_error(std::string("Error parsing file: ") + filePath);
  }

  //TODO: verify schema

  return jsonRoot;
}

}

//------------------------------------------------------------------------------
vtkMRMLNodeNewMacro(vtkMRMLSRepStorageNode);

//----------------------------------------------------------------------------
vtkMRMLSRepStorageNode::vtkMRMLSRepStorageNode()
  : vtkMRMLStorageNode()
  , CoordinateSystemWrite(vtkMRMLStorageNode::CoordinateSystemLPS)
{
  this->DefaultWriteFileExtension = "srep.json";
}

//----------------------------------------------------------------------------
vtkMRMLSRepStorageNode::~vtkMRMLSRepStorageNode() = default;

//----------------------------------------------------------------------------
void vtkMRMLSRepStorageNode::CoordinateSystemWriteRASOn() {
  this->SetCoordinateSystemWrite(vtkMRMLStorageNode::CoordinateSystemRAS);
}

//----------------------------------------------------------------------------
void vtkMRMLSRepStorageNode::CoordinateSystemWriteLPSOn() {
  this->SetCoordinateSystemWrite(vtkMRMLStorageNode::CoordinateSystemLPS);
}

//----------------------------------------------------------------------------
void vtkMRMLSRepStorageNode::SetCoordinateSystemWrite(vtkMRMLSRepStorageNode::SRepCoordinateSystemType system) {
  this->CoordinateSystemWrite = system;
  this->Modified();
}

//----------------------------------------------------------------------------
vtkMRMLSRepStorageNode::SRepCoordinateSystemType vtkMRMLSRepStorageNode::GetCoordinateSystemWrite() const {
  return this->CoordinateSystemWrite;
}

//----------------------------------------------------------------------------
bool vtkMRMLSRepStorageNode::CanReadInReferenceNode(vtkMRMLNode *refNode)
{
  return refNode->IsA("vtkMRMLSRepNode");
}

//----------------------------------------------------------------------------
std::string vtkMRMLSRepStorageNode::GetSRepType() {
  const char* filePath = this->GetFileName();
  if (!filePath)
    {
    vtkErrorMacro("vtkMRMLSRepStorageNode::GetSRepType failed: invalid filename");
    return "";
    }

  try {
    auto jsonRoot = CreateJsonDocumentFromFile(filePath);

    if (jsonRoot->HasMember(keys::RectangularGridSRep)) {
      return "vtkMRMLRectangularGridSRepNode";
    } else if (jsonRoot->HasMember(keys::EllipticalSRep)) {
      return "vtkMRMLEllipticalSRepNode";
    } else {
      vtkErrorMacro("vtkMRMLSRepStorageNode::GetSRepType failed: unable to find valid srep type");
      return "";
    }
  } catch (const std::exception& e) {
    vtkErrorMacro("vtkMRMLSRepStorageNode::GetSRepType failed: " << e.what());
    return "";
  }
}

//----------------------------------------------------------------------------
int vtkMRMLSRepStorageNode::ReadDataInternal(vtkMRMLNode * refNode)
{
  constexpr int success = 1;
  constexpr int failure = 0;

  if (!refNode)
    {
    vtkErrorMacro("vtkMRMLSRepStorageNode:ReadDataInternal: null reference node!");
    return failure;
    }

  auto srepNode = vtkMRMLSRepNode::SafeDownCast(refNode);
  if (!srepNode)
    {
    vtkErrorMacro("vtkMRMLSRepStorageNode:ReadDataInternal: refNode not a vtkMRMLSRepNode");
    return failure;
    }

  const char* filePath = this->GetFileName();
  if (!filePath)
    {
    vtkErrorMacro("vtkMRMLSRepStorageNode::ReadDataInternal failed: invalid filename");
    return failure;
    }

  try {
    auto jsonRootPtr = CreateJsonDocumentFromFile(filePath);
    auto& jsonRoot = *jsonRootPtr;

    if (jsonRoot.HasMember(keys::RectangularGridSRep)) {
      read(jsonRoot[keys::RectangularGridSRep], vtkMRMLRectangularGridSRepNode::SafeDownCast(srepNode));
    } else if (jsonRoot.HasMember(keys::EllipticalSRep)) {
      read(jsonRoot[keys::EllipticalSRep], vtkMRMLEllipticalSRepNode::SafeDownCast(srepNode));
    } else {
      vtkErrorMacro("vtkMRMLSRepStorageNode::ReadDataInternal failed because no known srep found");
      return failure;
    }

    auto displayIter = jsonRoot.FindMember(keys::Display);
    if (displayIter != jsonRoot.MemberEnd()) {
      if (!srepNode->GetDisplayNode()) {
        srepNode->CreateDefaultDisplayNodes();
      }
      read(displayIter->value, *srepNode->GetSRepDisplayNode());
    }
    return success;
  } catch (const std::exception& e) {
    vtkErrorMacro("vtkMRMLSRepStorageNode::ReadDataInternal failed: " << e.what());
    return failure;
  }
}

//----------------------------------------------------------------------------
int vtkMRMLSRepStorageNode::WriteDataInternal(vtkMRMLNode *refNode)
{
  constexpr int success = 1;
  constexpr int failure = 0;

  const std::string fullName = this->GetFullNameFromFileName();
  if (fullName.empty())
    {
    vtkErrorMacro("vtkMRMLSRepJsonStorageNode::WriteDataInternal: Writing srep node file failed: file name not specified.");
    return failure;
    }
  vtkDebugMacro("WriteDataInternal: have file name " << fullName.c_str());

  auto* srepNode = vtkMRMLSRepNode::SafeDownCast(refNode);
  if (!srepNode) {
    vtkErrorMacro("vtkMRMLSRepJsonStorageNode::WriteDataInternal: Writing srep node file failed: unable to cast input node "
      << refNode->GetID() << " to a known srep node.");
    return failure;
  }

  FILE* fp = fopen(fullName.c_str(), "wb");
  const auto closeFp = finally([fp](){
    fclose(fp);
  });

  // Prepare JSON writer and output stream.
  std::array<char, BufferSize> writeBuffer;
  rapidjson::FileWriteStream os(fp, writeBuffer.data(), writeBuffer.size());
  rapidjson::PrettyWriter<rapidjson::FileWriteStream> writer(os);

  //TODO:
  // writer.StartObject();
  // writer.Key("@schema"); writer.String(SREP_SCHEMA.c_str());

  writer.StartObject();

  // cast the input node
  if (auto rectangularGridNode = vtkMRMLRectangularGridSRepNode::SafeDownCast(refNode)) {
    write(writer, *rectangularGridNode, this->CoordinateSystemWrite);
  } else if (auto ellipticalNode = vtkMRMLEllipticalSRepNode::SafeDownCast(refNode)) {
    write(writer, *ellipticalNode, this->CoordinateSystemWrite);
  } else {
    vtkErrorMacro("vtkMRMLSRepJsonStorageNode::WriteDataInternal: Writing srep node file failed: unable to cast input node "
      << refNode->GetID() << " to a known srep node.");
    return failure;
  }

  auto* displayNode = srepNode->GetSRepDisplayNode();
  if (displayNode) {
    write(writer, *displayNode);
  }

  writer.EndObject();
  return success;
}

//----------------------------------------------------------------------------
void vtkMRMLSRepStorageNode::InitializeSupportedReadFileTypes()
{
  this->SupportedReadFileTypes->InsertNextValue("SRep JSON (.srep.json)");
}

//----------------------------------------------------------------------------
void vtkMRMLSRepStorageNode::InitializeSupportedWriteFileTypes()
{
  this->SupportedWriteFileTypes->InsertNextValue("SRep JSON (.srep.json)");
}

//----------------------------------------------------------------------------
vtkMRMLSRepNode* vtkMRMLSRepStorageNode::CreateSRepNode(const char* nodeName) {
  vtkMRMLScene* scene = this->GetScene();
  if (!scene) {
    vtkErrorMacro("vtkMRMLMarkupsJsonStorageNode::CreateSRepNode failed: invalid scene");
    return nullptr;
  }

  const auto srepType = this->GetSRepType();
  if (srepType.empty()) {
    return nullptr;
  }

  std::string newNodeName;
  if (nodeName && strlen(nodeName) > 0) {
    newNodeName = nodeName;
  } else {
    newNodeName = scene->GetUniqueNameByString(this->GetFileNameWithoutExtension(this->GetFileName()).c_str());
  }

  auto* srepNode = vtkMRMLSRepNode::SafeDownCast(scene->AddNewNodeByClass(srepType, newNodeName));
  if (!srepNode) {
    vtkErrorMacro("vtkMRMLMarkupsJsonStorageNode::CreateSRepNode failed: unable to make class by name: " << srepType);
    return nullptr;
  }
  this->ReadData(srepNode);
  if (!srepNode->GetDisplayNode()) {
    srepNode->CreateDefaultDisplayNodes();
  }
  return srepNode;
}
