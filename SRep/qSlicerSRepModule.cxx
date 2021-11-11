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

// SRep includes
#include "qSlicerSRepModule.h"
#include "qSlicerSRepModuleWidget.h"
#include "qSlicerSRepReader.h"
#include "qSlicerSRepWriter.h"
#include <vtkMRMLSRepDisplayableManager.h>
#include <vtkSlicerSRepLogic.h>

#include <qSlicerCoreApplication.h>
#include <qSlicerModuleManager.h>
#include <qSlicerIOManager.h>

#include <vtkMRMLSliceViewDisplayableManagerFactory.h>
#include <vtkMRMLThreeDViewDisplayableManagerFactory.h>

// SubjectHierarchy Plugins includes
#include "qSlicerSubjectHierarchyPluginHandler.h"
#include "qSlicerSubjectHierarchySRepPlugin.h"

// DisplayableManager initialization
#include <vtkAutoInit.h>
VTK_MODULE_INIT(vtkSlicerSRepModuleMRMLDisplayableManager)

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerSRepModulePrivate
{
public:
  qSlicerSRepModulePrivate();
};

//-----------------------------------------------------------------------------
// qSlicerSRepModulePrivate methods

//-----------------------------------------------------------------------------
qSlicerSRepModulePrivate::qSlicerSRepModulePrivate()
{
}

//-----------------------------------------------------------------------------
// qSlicerSRepModule methods

//-----------------------------------------------------------------------------
qSlicerSRepModule::qSlicerSRepModule(QObject* _parent)
  : Superclass(_parent)
  , d_ptr(new qSlicerSRepModulePrivate)
{
}

//-----------------------------------------------------------------------------
qSlicerSRepModule::~qSlicerSRepModule()
{
}

//-----------------------------------------------------------------------------
QString qSlicerSRepModule::helpText() const
{
  return "This module is for loading and visualizing SReps";
}

//-----------------------------------------------------------------------------
QString qSlicerSRepModule::acknowledgementText() const
{
  return "This work was partially funded by NIH grant TBD";
}

//-----------------------------------------------------------------------------
QStringList qSlicerSRepModule::contributors() const
{
  QStringList moduleContributors;
  moduleContributors << QString("Connor Bowley (Kitware, Inc.)");
  return moduleContributors;
}

//-----------------------------------------------------------------------------
QIcon qSlicerSRepModule::icon() const
{
  return QIcon(":/Icons/SRep.png");
}

//-----------------------------------------------------------------------------
QStringList qSlicerSRepModule::categories() const
{
  return QStringList() << "SRep2";
}

//-----------------------------------------------------------------------------
QStringList qSlicerSRepModule::dependencies() const
{
  return QStringList();
}

//-----------------------------------------------------------------------------
void qSlicerSRepModule::setup()
{
  this->Superclass::setup();

  // Register displayable manager with 3d view
  vtkMRMLThreeDViewDisplayableManagerFactory::GetInstance()->RegisterDisplayableManager("vtkMRMLSRepDisplayableManager");

  // Register Subject Hierarchy core plugins
  qSlicerSubjectHierarchyPluginHandler::instance()->registerPlugin(new qSlicerSubjectHierarchySRepPlugin());

  // Register Reader/Writer
  qSlicerIOManager* ioManager = qSlicerApplication::application()->ioManager();
  // ioManager takes ownership of pointer
  ioManager->registerIO(new qSlicerSRepReader(vtkSlicerSRepLogic::SafeDownCast(this->logic()), this));
  ioManager->registerIO(new qSlicerSRepWriter(this));
}

//-----------------------------------------------------------------------------
qSlicerAbstractModuleRepresentation* qSlicerSRepModule
::createWidgetRepresentation()
{
  return new qSlicerSRepModuleWidget;
}

//-----------------------------------------------------------------------------
vtkMRMLAbstractLogic* qSlicerSRepModule::createLogic()
{
  return vtkSlicerSRepLogic::New();
}

//-----------------------------------------------------------------------------
QStringList qSlicerSRepModule::associatedNodeTypes() const {
  return QStringList()
    << "vtkMRMLRectangularGridSRepNode"
    << "vtkMRMLSRepDisplayNode"
    << "vtkMRMLSRepStorageNode"
  ;
}
