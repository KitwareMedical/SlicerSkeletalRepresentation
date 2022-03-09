/*==============================================================================

  Copyright (c) Laboratory for Percutaneous Surgery (PerkLab)
  Queen's University, Kingston, ON, Canada. All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

  This file was originally developed by Kyle Sunderland, PerkLab, Queen's University
  and was supported through CANARIE's Research Software Program, Cancer
  Care Ontario, OpenAnatomy, and Brigham and Women�s Hospital through NIH grant R01MH112748.

==============================================================================*/

// Qt includes
#include <QDebug>

#include "qSlicerSRepWriter.h"

// QTCore includes
#include "qSlicerCoreApplication.h"
#include "qSlicerCoreIOManager.h"

// MRML includes
#include <vtkMRMLScene.h>
#include <vtkMRMLSceneViewNode.h>
#include <vtkMRMLSRepStorageNode.h>
#include <vtkMRMLStorableNode.h>

// VTK includes
#include <vtkStdString.h>
#include <vtkStringArray.h>

//----------------------------------------------------------------------------
qSlicerSRepWriter::qSlicerSRepWriter(QObject* parentObject)
  : qSlicerNodeWriter("SRep", QString("SRepFile"), QStringList() << "vtkMRMLSRepNode", true, parentObject)
{
}

//----------------------------------------------------------------------------
qSlicerSRepWriter::~qSlicerSRepWriter() = default;

//----------------------------------------------------------------------------
QStringList qSlicerSRepWriter::extensions(vtkObject* vtkNotUsed(object))const
{
  return QStringList() << "SRep (*.srep.json)";
}

//----------------------------------------------------------------------------
bool qSlicerSRepWriter::write(const qSlicerIO::IOProperties& properties)
{
  Q_ASSERT(!properties["nodeID"].toString().isEmpty());

  vtkMRMLStorableNode* node = vtkMRMLStorableNode::SafeDownCast(
    this->getNodeByID(properties["nodeID"].toString().toUtf8().data()));
  if (!this->canWriteObject(node)) {
    return false;
  }

  vtkMRMLSRepNode* srepNode = vtkMRMLSRepNode::SafeDownCast(node);
  if (!srepNode) {
    qDebug() << "SRep writer can only write vtkMRMLSRepNodes";
    return false;
  }
  auto storageNode = vtkSmartPointer<vtkMRMLStorageNode>::Take(srepNode->CreateDefaultStorageNode());
  if (!storageNode) {
    qDebug() << "Error creating srep storage node";
    return false;
  }
  srepNode->SetAndObserveStorageNodeID(storageNode->GetID());

  return Superclass::write(properties);
}
