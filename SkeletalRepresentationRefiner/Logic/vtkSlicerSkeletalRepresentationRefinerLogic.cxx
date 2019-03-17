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

// SkeletalRepresentationRefiner Logic includes
#include "vtkSlicerSkeletalRepresentationRefinerLogic.h"

// MRML includes
#include <vtkMRMLScene.h>

// VTK includes
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>

// STD includes
#include <cassert>

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkSlicerSkeletalRepresentationRefinerLogic);

//----------------------------------------------------------------------------
vtkSlicerSkeletalRepresentationRefinerLogic::vtkSlicerSkeletalRepresentationRefinerLogic()
{
}

//----------------------------------------------------------------------------
vtkSlicerSkeletalRepresentationRefinerLogic::~vtkSlicerSkeletalRepresentationRefinerLogic()
{
}

//----------------------------------------------------------------------------
void vtkSlicerSkeletalRepresentationRefinerLogic::PrintSelf(ostream& os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os, indent);
}

void vtkSlicerSkeletalRepresentationRefinerLogic::SetImageFileName(const std::__cxx11::string &imageFilePath)
{
    mImageFilePath = imageFilePath;
}

void vtkSlicerSkeletalRepresentationRefinerLogic::SetSrepFileName(const std::__cxx11::string &srepFilePath)
{
    mSrepFilePath = srepFilePath;
}

void vtkSlicerSkeletalRepresentationRefinerLogic::Refine()
{

}

//---------------------------------------------------------------------------
void vtkSlicerSkeletalRepresentationRefinerLogic::SetMRMLSceneInternal(vtkMRMLScene * newScene)
{
  vtkNew<vtkIntArray> events;
  events->InsertNextValue(vtkMRMLScene::NodeAddedEvent);
  events->InsertNextValue(vtkMRMLScene::NodeRemovedEvent);
  events->InsertNextValue(vtkMRMLScene::EndBatchProcessEvent);
  this->SetAndObserveMRMLSceneEventsInternal(newScene, events.GetPointer());
}

//-----------------------------------------------------------------------------
void vtkSlicerSkeletalRepresentationRefinerLogic::RegisterNodes()
{
  assert(this->GetMRMLScene() != 0);
}

//---------------------------------------------------------------------------
void vtkSlicerSkeletalRepresentationRefinerLogic::UpdateFromMRMLScene()
{
  assert(this->GetMRMLScene() != 0);
}

//---------------------------------------------------------------------------
void vtkSlicerSkeletalRepresentationRefinerLogic
::OnMRMLSceneNodeAdded(vtkMRMLNode* vtkNotUsed(node))
{
}

//---------------------------------------------------------------------------
void vtkSlicerSkeletalRepresentationRefinerLogic
::OnMRMLSceneNodeRemoved(vtkMRMLNode* vtkNotUsed(node))
{
}
