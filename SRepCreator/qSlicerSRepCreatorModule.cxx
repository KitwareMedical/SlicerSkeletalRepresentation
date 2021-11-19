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

// SRepCreator Logic includes
#include <vtkSlicerSRepCreatorLogic.h>

// SRepCreator includes
#include "qSlicerSRepCreatorModule.h"
#include "qSlicerSRepCreatorModuleWidget.h"

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerSRepCreatorModulePrivate
{
public:
  qSlicerSRepCreatorModulePrivate();
};

//-----------------------------------------------------------------------------
// qSlicerSRepCreatorModulePrivate methods

//-----------------------------------------------------------------------------
qSlicerSRepCreatorModulePrivate::qSlicerSRepCreatorModulePrivate()
{
}

//-----------------------------------------------------------------------------
// qSlicerSRepCreatorModule methods

//-----------------------------------------------------------------------------
qSlicerSRepCreatorModule::qSlicerSRepCreatorModule(QObject* _parent)
  : Superclass(_parent)
  , d_ptr(new qSlicerSRepCreatorModulePrivate)
{
}

//-----------------------------------------------------------------------------
qSlicerSRepCreatorModule::~qSlicerSRepCreatorModule()
{
}

//-----------------------------------------------------------------------------
QString qSlicerSRepCreatorModule::helpText() const
{
  return "This module is used to create SReps from existing models";
}

//-----------------------------------------------------------------------------
QString qSlicerSRepCreatorModule::acknowledgementText() const
{
  return "This work was partially funded by NIH grant NXNNXXNNNNNN-NNXN";
}

//-----------------------------------------------------------------------------
QStringList qSlicerSRepCreatorModule::contributors() const
{
  QStringList moduleContributors;
  moduleContributors << QString("John Doe (AnyWare Corp.)");
  return moduleContributors;
}

//-----------------------------------------------------------------------------
QIcon qSlicerSRepCreatorModule::icon() const
{
  return QIcon(":/Icons/SRepCreator.png");
}

//-----------------------------------------------------------------------------
QStringList qSlicerSRepCreatorModule::categories() const
{
  return QStringList() << "SRep2";
}

//-----------------------------------------------------------------------------
QStringList qSlicerSRepCreatorModule::dependencies() const
{
  return QStringList();
}

//-----------------------------------------------------------------------------
void qSlicerSRepCreatorModule::setup()
{
  this->Superclass::setup();
}

//-----------------------------------------------------------------------------
qSlicerAbstractModuleRepresentation* qSlicerSRepCreatorModule
::createWidgetRepresentation()
{
  return new qSlicerSRepCreatorModuleWidget;
}

//-----------------------------------------------------------------------------
vtkMRMLAbstractLogic* qSlicerSRepCreatorModule::createLogic()
{
  return vtkSlicerSRepCreatorLogic::New();
}
