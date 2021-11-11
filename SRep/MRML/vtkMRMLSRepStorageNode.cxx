#include "vtkMRMLSRepStorageNode.h"
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
  const char * const Rows = "Rows";
  const char * const Cols = "Cols";
  const char * const Skeleton = "Skeleton";
  const char * const UpSpoke = "UpSpoke";
  const char * const DownSpoke = "DownSpoke";
  const char * const CrestSpoke = "CrestSpoke";
  const char * const Direction = "Direction";
  const char * const SkeletalPoint = "SkeletalPoint";

  const char * const Display = "Display";
  const char * const Visibility = "Visibility";
  const char * const Opacity = "Opacity";
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

void write(rapidjson::PrettyWriter<rapidjson::FileWriteStream>& writer, const srep::Point3d& point) {
  writeSingleLineArray(writer, point.AsArray());
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
  return srep::Point3d(read3dArray(json));
}

void write(rapidjson::PrettyWriter<rapidjson::FileWriteStream>& writer, const srep::Vector3d& vector) {
  writeSingleLineArray(writer, vector.AsArray());
}

srep::Vector3d readVector3d(rapidjson::Value& json) {
  return srep::Vector3d(read3dArray(json));
}

void write(rapidjson::PrettyWriter<rapidjson::FileWriteStream>& writer, const srep::Spoke& spoke) {
  writer.StartObject();
  writer.Key(keys::SkeletalPoint);
  write(writer, spoke.GetSkeletalPoint());
  writer.Key(keys::Direction);
  write(writer, spoke.GetDirection());
  writer.EndObject();
}

srep::Spoke readSpoke(rapidjson::Value& json) {
  auto skeletalIter = SafeFindMember(json, keys::SkeletalPoint);
  auto directionIter = SafeFindMember(json, keys::Direction);
  return srep::Spoke(readPoint3d(skeletalIter->value), readVector3d(directionIter->value));
}

void write(rapidjson::PrettyWriter<rapidjson::FileWriteStream>& writer, vtkMRMLRectangularGridSRepNode& mrmlSRep)
{
  const srep::RectangularGridSRep* srep = mrmlSRep.GetRectangularGridSRep();
  //writing out rows and cols more for people who look at the JSON
  const auto& rows = srep ? srep->GetNumRows() : 0;
  const auto& cols = srep ? srep->GetNumCols() : 0;
  const auto& grid = srep ? srep->GetSkeletalPoints() : srep::RectangularGridSRep::SkeletalGrid{};

  writer.Key(keys::RectangularGridSRep);
  writer.StartObject();
  writer.Key(keys::Rows); writer.Uint(rows);
  writer.Key(keys::Cols); writer.Uint(cols);

  //write skeletal points as an array of arrays
  writer.Key(keys::Skeleton);
  writer.StartArray();
  for (size_t row  = 0; row < grid.size(); ++row) {
    writer.StartArray();
    for (size_t col = 0; col < grid[row].size(); ++col) {
      const auto& skeletalPoint = grid[row][col];

      writer.StartObject();
      writer.Key(keys::UpSpoke);
      write(writer, skeletalPoint.GetUpSpoke());
      writer.Key(keys::DownSpoke);
      write(writer, skeletalPoint.GetDownSpoke());
      if (skeletalPoint.IsCrest()) {
        writer.Key(keys::CrestSpoke);
        write(writer, skeletalPoint.GetCrestSpoke());
      }
      writer.EndObject();
    }
    writer.EndArray();
  }
  writer.EndArray();

  writer.EndObject();
}

void read(rapidjson::Value& json, vtkMRMLRectangularGridSRepNode* rectangularGridSRep) {
  if (!rectangularGridSRep) {
    throw std::invalid_argument("Node is not a vtkMRMLRectangularGridSRepNode");
  }

  auto skeletonIter = SafeFindMember(json, keys::Skeleton);
  auto& skeleton = skeletonIter->value;

  srep::RectangularGridSRep::SkeletalGrid grid;
  for (auto& row : skeleton.GetArray()) {
    grid.push_back(std::vector<srep::SkeletalPoint>{});
    if (!row.IsArray()) {
      throw std::runtime_error("Error parsing vtkMRMLRectangularGridSRepNode JSON. Row is not array.");
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

  // RectangularGridSRep constructor will throw if there are bad things like crest spokes not on the crest
  rectangularGridSRep->SetRectangularGridSRep(std::unique_ptr<srep::RectangularGridSRep>(new srep::RectangularGridSRep(grid)));
}

void write(rapidjson::PrettyWriter<rapidjson::FileWriteStream>& writer, vtkMRMLSRepDisplayNode& displayNode) {
  writer.Key(keys::Display);
  writer.StartObject();
  writer.Key(keys::Visibility);
  writer.Bool(static_cast<bool>(displayNode.GetVisibility()));
  writer.Key(keys::Opacity);
  writer.Double(displayNode.GetOpacity());
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
{
  this->DefaultWriteFileExtension = "srep.json";
}

//----------------------------------------------------------------------------
vtkMRMLSRepStorageNode::~vtkMRMLSRepStorageNode() = default;

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
    write(writer, *rectangularGridNode);
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
