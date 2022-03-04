/*==============================================================================

  Program: 3D Slicer

  Portions (c) Copyright Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

==============================================================================*/

// SRep Logic includes
#include "vtkSlicerSRepLogic.h"

// MRML includes
#include <vtkMRMLScene.h>
#include <vtkMRMLSRepDisplayNode.h>
#include <vtkMRMLSRepStorageNode.h>

// VTK includes
#include <vtkCellData.h>
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>

// STD includes
#include <cassert>

#include "vtkMRMLSRepNode.h"
#include "SRepInterpolation.h"

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkSlicerSRepLogic);

//----------------------------------------------------------------------------
vtkSlicerSRepLogic::vtkSlicerSRepLogic()
{
}

//----------------------------------------------------------------------------
vtkSlicerSRepLogic::~vtkSlicerSRepLogic()
{
}

//----------------------------------------------------------------------------
void vtkSlicerSRepLogic::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

//---------------------------------------------------------------------------
void vtkSlicerSRepLogic::SetMRMLSceneInternal(vtkMRMLScene * newScene)
{
  vtkNew<vtkIntArray> events;
  events->InsertNextValue(vtkMRMLScene::NodeAddedEvent);
  events->InsertNextValue(vtkMRMLScene::NodeRemovedEvent);
  events->InsertNextValue(vtkMRMLScene::EndBatchProcessEvent);
  this->SetAndObserveMRMLSceneEventsInternal(newScene, events.GetPointer());
}

//-----------------------------------------------------------------------------
void vtkSlicerSRepLogic::RegisterNodes()
{
  assert(this->GetMRMLScene() != 0);

  vtkMRMLScene *scene = this->GetMRMLScene();

  scene->RegisterNodeClass(vtkSmartPointer<vtkMRMLEllipticalSRepNode>::New());
  scene->RegisterNodeClass(vtkSmartPointer<vtkMRMLSRepDisplayNode>::New());
  scene->RegisterNodeClass(vtkSmartPointer<vtkMRMLSRepStorageNode>::New());
}

//---------------------------------------------------------------------------
void vtkSlicerSRepLogic::UpdateFromMRMLScene()
{
  assert(this->GetMRMLScene() != 0);
}

//---------------------------------------------------------------------------
void vtkSlicerSRepLogic
::OnMRMLSceneNodeAdded(vtkMRMLNode* vtkNotUsed(node))
{
}

//---------------------------------------------------------------------------
void vtkSlicerSRepLogic
::OnMRMLSceneNodeRemoved(vtkMRMLNode* vtkNotUsed(node))
{
}

//----------------------------------------------------------------------------
std::string vtkSlicerSRepLogic::AddNewEllipticalSRepNode(const std::string& name, vtkMRMLScene* scene) {
  std::string id;
  if (!scene && !this->GetMRMLScene()) {
    vtkErrorMacro("AddNewEllipticalSRepNode: no scene to add a srep node to!");
    return id;
  }

  vtkMRMLScene *addToThisScene = scene ? scene : this->GetMRMLScene();

  // create and add the node
  auto mnode = vtkSmartPointer<vtkMRMLEllipticalSRepNode>::New();
  addToThisScene->AddNode(mnode);

  // add a display node
  std::string displayID = this->AddFirstDisplayNodeForSRepNode(mnode);

  if (displayID.compare("") != 0) {
    // get the node id to return
    id = std::string(mnode->GetID());
    if (!name.empty()) {
      mnode->SetName(name.c_str());
    }
  }

  return id;
}

//----------------------------------------------------------------------------
std::string vtkSlicerSRepLogic::AddFirstDisplayNodeForSRepNode(vtkMRMLSRepNode *srepNode) {
  const std::string emptyId;
  if (!srepNode || !srepNode->GetScene()) {
    vtkErrorMacro("AddNewDisplayNodeForSRepNode: unable to add a srep display node!");
    return emptyId;
  }

  if (srepNode->GetDisplayNode()) {
    return srepNode->GetDisplayNodeID();
  }

  srepNode->CreateDefaultDisplayNodes();
  auto* displayNode = vtkMRMLSRepDisplayNode::SafeDownCast(srepNode->GetDisplayNode());
  if (!displayNode) {
    vtkErrorMacro("AddNewDisplayNodeForSRepNode: error creating new display node");
    return emptyId;
  }

  return displayNode->GetID();
}

//----------------------------------------------------------------------------
const char* vtkSlicerSRepLogic::LoadSRep(const char* fileName, const char* nodeName) {
  if (!fileName) {
    vtkErrorMacro("LoadSRep: null file, cannot load");
    return nullptr;
  }

  vtkDebugMacro("LoadSRep, file name = " << fileName << ", nodeName = " << (nodeName ? nodeName : "null"));

  vtkMRMLSRepStorageNode* storageNode = vtkMRMLSRepStorageNode::SafeDownCast(
    this->GetMRMLScene()->AddNewNodeByClass("vtkMRMLSRepStorageNode"));

  if (!storageNode) {
    vtkErrorMacro("LoadSRep: failed to instantiate srep storage node by class vtkMRMLSRepStorageNode");
    return nullptr;
  }

  storageNode->SetFileName(fileName);
  vtkMRMLSRepNode* srepNode = storageNode->CreateSRepNode(nodeName);
  if (!srepNode) {
    return nullptr;
  }
  return srepNode->GetID();

}

//----------------------------------------------------------------------------
std::string vtkSlicerSRepLogic::InterpolateSRep(vtkMRMLEllipticalSRepNode* srepNode, size_t interpolationlevel, const std::string& newNodeName) {
  auto scene = this->GetMRMLScene();
  if (!scene) {
    vtkErrorMacro("InterpolateSRep: no scene to add a srep node to!");
    return "";
  }

  const auto nodeID = AddNewEllipticalSRepNode(newNodeName, scene);
  if (nodeID.empty()) {
    vtkErrorMacro("InterpolateSRep: Error making Elliptical SRep node");
    return "";
  }

  auto interpolatedSRepNode = vtkMRMLEllipticalSRepNode::SafeDownCast(scene->GetNodeByID(nodeID));
  if (!interpolatedSRepNode) {
    vtkErrorMacro("InterpolateSRep: Unable to find newly created SRep node: " << nodeID);
    return "";
  }

  const bool success = this->InterpolateSRep(srepNode, interpolationlevel, interpolatedSRepNode);
  if (!success) {
    scene->RemoveNode(interpolatedSRepNode);
    return "";
  }
  return nodeID;
}

//----------------------------------------------------------------------------
bool vtkSlicerSRepLogic::InterpolateSRep(vtkMRMLEllipticalSRepNode* srepNode, size_t interpolationlevel, vtkMRMLEllipticalSRepNode* destination) {
  if (!destination) {
    vtkErrorMacro("InterpolateSRep: no destination");
    return false;
  }
  
  if (!srepNode) {
    vtkErrorMacro("InterpolateSRep: input node is nullptr");
    return false;
  }

  auto srep = srepNode->GetEllipticalSRep();
  if (!srep) {
    vtkErrorMacro("InterpolateSRep: input node does not have an SRep");
    return false;
  }

  if (interpolationlevel == 0) {
    destination->SetEllipticalSRep(srep->SmartClone());
    return true;
  } else {

    auto interpolatedSRep = this->SmartInterpolateSRep(*srep, interpolationlevel);
    if (!interpolatedSRep) {
      vtkErrorMacro("InterpolateSRep: Unable to interpolate SRep");
      return false;
    }

    destination->SetEllipticalSRep(interpolatedSRep);
    return true;
  }
}

//----------------------------------------------------------------------------
vtkEllipticalSRep* vtkSlicerSRepLogic::InterpolateSRep(const vtkEllipticalSRep& srep, size_t interpolationlevel) {
  auto ret = SmartInterpolateSRep(srep, interpolationlevel);
  if (ret) {
    ret->Register(nullptr);
  }
  return ret;
}

//----------------------------------------------------------------------------
vtkSmartPointer<vtkEllipticalSRep> vtkSlicerSRepLogic::SmartInterpolateSRep(const vtkEllipticalSRep& srep, size_t interpolationlevel) {
  return sreplogic::SmartInterpolateSRep(interpolationlevel, srep);
}

namespace {
struct vtkSpokeIds {
  vtkIdType boundaryId;
  vtkIdType skeletonId;
};
}

//----------------------------------------------------------------------------
vtkSmartPointer<vtkPolyData> vtkSlicerSRepLogic::SmartExportSRepToPolyData(const vtkMeshSRepInterface& srep, const vtkSRepExportPolyDataProperties& properties) {
  auto polyData = vtkSmartPointer<vtkPolyData>::New();

  vtkNew<vtkPoints> points;
  polyData->SetPoints(points);

  vtkNew<vtkCellArray> lines;
  polyData->SetLines(lines);

  auto srepArray = properties.GetSRepDataArray();
  vtkSmartPointer<vtkDataArray> pointDataArray;
  vtkSmartPointer<vtkDataArray> lineDataArray;
  if (srepArray) {
    pointDataArray = vtkSmartPointer<vtkDataArray>::Take(srepArray->NewInstance());
    pointDataArray->SetNumberOfComponents(srepArray->GetNumberOfComponents());
    pointDataArray->SetName(properties.GetPointTypeArrayName().c_str());
    polyData->GetPointData()->SetScalars(pointDataArray);

    lineDataArray = vtkSmartPointer<vtkDataArray>::Take(srepArray->NewInstance());
    lineDataArray->SetNumberOfComponents(srepArray->GetNumberOfComponents());
    lineDataArray->SetName(properties.GetLineTypeArrayName().c_str());
    polyData->GetCellData()->SetScalars(lineDataArray);
  }

  //-------------------------------
  const auto insertNextScalarData = [&srepArray](vtkDataArray* dest, int srepDataType) {
    if (srepArray) {
      dest->InsertNextTuple(srepDataType, srepArray);
    }
  };

  //-------------------------------
  const auto insertNextPoint = [&points, &pointDataArray, insertNextScalarData](const srep::Point3d& point, int pointType) {
    const auto id = points->InsertNextPoint(point.AsArray().data());
    insertNextScalarData(pointDataArray, pointType);
    return id;
  };

  //-------------------------------
  const auto insertNextLine = [&lines, &lineDataArray, insertNextScalarData](const vtkIdType start, const vtkIdType end, int lineType) {
    lines->InsertNextCell(2);
    lines->InsertCellPoint(start);
    lines->InsertCellPoint(end);
    insertNextScalarData(lineDataArray, lineType);
  };

  //-------------------------------
  const auto addSpokeMesh = [insertNextPoint, insertNextLine]
    (const vtkSRepSpokeMesh& mesh,
     int skeletonPointType,
     int boundaryPointType,
     bool addSpokes,
     int spokeType,
     bool addConnections,
     int connectionType,
     const std::vector<vtkMeshSRepInterface::IndexType> spine,
     std::function<bool(long i)> forceAddSkeletalPoint) -> std::vector<vtkSpokeIds>
  {
    std::vector<vtkSpokeIds> spokesToVTKPointIds;

    const auto isSpine = [&spine](long i) {
      return std::find(spine.begin(), spine.end(), i) != spine.end();
    };

    // add all the points and the spoke lines
    for (long i = 0; i < mesh.GetNumberOfSpokes(); ++i) {
      vtkSpokeIds ids;
      if (addSpokes || addConnections || isSpine(i) || forceAddSkeletalPoint(i)) {
        ids.skeletonId = insertNextPoint(mesh[i]->GetSkeletalPoint(), skeletonPointType);
      } else {
        ids.skeletonId = -1;
      }
      if (addSpokes) {
        ids.boundaryId = insertNextPoint(mesh[i]->GetBoundaryPoint(), boundaryPointType);
        insertNextLine(ids.skeletonId, ids.boundaryId, spokeType);
      } else {
        ids.boundaryId = -1;
      }
      spokesToVTKPointIds.push_back(ids);
    }

    // add the connection lines. It is essentially a bidirectional graph, so only add one line between two points, even if
    // it shows up twice "once in each direction"
    std::set<std::pair<vtkIdType, vtkIdType>> spineConnections;
    for (size_t i = 1; i < spine.size(); ++i) {
      //insert a sorted pair
      vtkIdType point1 = spokesToVTKPointIds[spine[i-1]].skeletonId;
      vtkIdType point2 = spokesToVTKPointIds[spine[i]].skeletonId;
      //sort the points
      if (point1 > point2) {
        std::swap(point1, point2);
      }
      spineConnections.insert(std::make_pair(point1, point2));
    }

    for (const auto& spineConnection : spineConnections) {
      insertNextLine(
        spineConnection.first,
        spineConnection.second,
        vtkSRepExportPolyDataProperties::SpineLineType);
    }

    if (addConnections) {
      std::set<std::pair<vtkIdType, vtkIdType>> connections;
      for (long i = 0; i < mesh.GetNumberOfSpokes(); ++i) {
        const auto neighbors = mesh.GetNeighbors(i);
        for (size_t neighbor : neighbors) {
          vtkIdType point1 = spokesToVTKPointIds[i].skeletonId;
          vtkIdType point2 = spokesToVTKPointIds[neighbor].skeletonId;
          //sort the points
          if (point1 > point2) {
            std::swap(point1, point2);
          }
          connections.insert(std::make_pair(point1, point2));
        }
      }

      std::set<std::pair<vtkIdType, vtkIdType>> spinelessConnections;
      std::set_difference(
        connections.begin(), connections.end(),
        spineConnections.begin(), spineConnections.end(),
        std::inserter(spinelessConnections, spinelessConnections.begin()));

      for (const auto& connection : spinelessConnections) {
        insertNextLine(connection.first, connection.second, connectionType);
      }
    }

    return spokesToVTKPointIds;
  };

  ///////////////////////////////////////
  // Start
  ///////////////////////////////////////
  const auto visibleUpSpineIndexes = properties.GetIncludeSpine() ? srep.GetUpSpine() : vtkSRepSpokeMesh::NeighborList{};
  const auto visibleDownSpineIndexes = properties.GetIncludeSpine() ? srep.GetDownSpine() : vtkSRepSpokeMesh::NeighborList{};

  const auto upSpokeToPointIds = addSpokeMesh(
    *srep.GetUpSpokes(),
    vtkSRepExportPolyDataProperties::UpSkeletalPointType, vtkSRepExportPolyDataProperties::UpBoundaryPointType,
    properties.GetIncludeUpSpokes(), vtkSRepExportPolyDataProperties::UpSpokeLineType,
    properties.GetIncludeSkeletalSheet(), vtkSRepExportPolyDataProperties::SkeletalSheetLineType,
    visibleUpSpineIndexes,
    [&](long i) {
      const auto crestSkeletonConnections = srep.GetCrestToUpSpokeConnections();
      return properties.GetIncludeSkeletonToCrestConnection()
        && crestSkeletonConnections.end() != std::find(crestSkeletonConnections.begin(), crestSkeletonConnections.end(), i);
    }
  );

  const auto downSpokeToPointIds = addSpokeMesh(
    *srep.GetDownSpokes(),
    vtkSRepExportPolyDataProperties::DownSkeletalPointType, vtkSRepExportPolyDataProperties::DownBoundaryPointType,
    properties.GetIncludeDownSpokes(), vtkSRepExportPolyDataProperties::DownSpokeLineType,
    properties.GetIncludeSkeletalSheet(), vtkSRepExportPolyDataProperties::SkeletalSheetLineType,
    visibleDownSpineIndexes,
    [&](long i) {
      const auto crestSkeletonConnections = srep.GetCrestToDownSpokeConnections();
      return properties.GetIncludeSkeletonToCrestConnection()
        && crestSkeletonConnections.end() != std::find(crestSkeletonConnections.begin(), crestSkeletonConnections.end(), i);
    }
  );

  const auto crestSpokeToPointIds = addSpokeMesh(
    *srep.GetCrestSpokes(),
    vtkSRepExportPolyDataProperties::CrestSkeletalPointType, vtkSRepExportPolyDataProperties::CrestBoundaryPointType,
    properties.GetIncludeCrestSpokes(), vtkSRepExportPolyDataProperties::CrestSpokeLineType,
    properties.GetIncludeCrestCurve(), vtkSRepExportPolyDataProperties::CrestCurveLineType,
    {},
    [&](long){ return properties.GetIncludeSkeletonToCrestConnection(); }
  );

  // connect the crest to skeleton
  if (properties.GetIncludeSkeletonToCrestConnection()) {
    for (size_t crestIndex = 0; crestIndex < srep.GetCrestToUpSpokeConnections().size(); ++crestIndex) {
      const auto skeletonIndex = srep.GetCrestToUpSpokeConnections()[crestIndex];
      insertNextLine(crestSpokeToPointIds[crestIndex].skeletonId, upSpokeToPointIds[skeletonIndex].skeletonId,
        vtkSRepExportPolyDataProperties::SkeletonToCrestConnectionLineType);
    }
    for (size_t crestIndex = 0; crestIndex < srep.GetCrestToDownSpokeConnections().size(); ++crestIndex) {
      const auto skeletonIndex = srep.GetCrestToDownSpokeConnections()[crestIndex];
      insertNextLine(crestSpokeToPointIds[crestIndex].skeletonId, downSpokeToPointIds[skeletonIndex].skeletonId,
        vtkSRepExportPolyDataProperties::SkeletonToCrestConnectionLineType);
    }
  }

  return polyData;
}

//----------------------------------------------------------------------------
vtkPolyData* vtkSlicerSRepLogic::ExportSRepToPolyData(const vtkMeshSRepInterface& srep, const vtkSRepExportPolyDataProperties& properties) {
  auto ret = SmartExportSRepToPolyData(srep, properties);
  if (ret) {
    ret->Register(nullptr);
  }
  return ret;
}
