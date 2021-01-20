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

// SRepRefiner Logic includes
#include <vtkSlicerSRepRefinerLogic.h>

// SRepRefiner includes
#include "qSlicerSRepRefinerModule.h"
#include "qSlicerSRepRefinerModuleWidget.h"

//-----------------------------------------------------------------------------
#if (QT_VERSION < QT_VERSION_CHECK(5, 0, 0))
#include <QtPlugin>
Q_EXPORT_PLUGIN2(qSlicerSRepRefinerModule, qSlicerSRepRefinerModule);
#endif

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerSRepRefinerModulePrivate
{
public:
  qSlicerSRepRefinerModulePrivate();
};

//-----------------------------------------------------------------------------
// qSlicerSRepRefinerModulePrivate methods

//-----------------------------------------------------------------------------
qSlicerSRepRefinerModulePrivate::qSlicerSRepRefinerModulePrivate()
{
}

//-----------------------------------------------------------------------------
// qSlicerSRepRefinerModule methods

//-----------------------------------------------------------------------------
qSlicerSRepRefinerModule::qSlicerSRepRefinerModule(QObject* _parent)
  : Superclass(_parent)
  , d_ptr(new qSlicerSRepRefinerModulePrivate)
{
}

//-----------------------------------------------------------------------------
qSlicerSRepRefinerModule::~qSlicerSRepRefinerModule()
{
}

//-----------------------------------------------------------------------------
QString qSlicerSRepRefinerModule::helpText() const
{
  return "This is a loadable module that can be bundled in an extension";
}

//-----------------------------------------------------------------------------
QString qSlicerSRepRefinerModule::acknowledgementText() const
{
  return "This file was originally developed by Zhiyuan Liu, and currently maintained by the SlicerSALT team.";
}

//-----------------------------------------------------------------------------
QStringList qSlicerSRepRefinerModule::contributors() const
{
  QStringList moduleContributors;
  moduleContributors << QString("Zhiyuan Liu, Megan Stuart, Jiyao Wang, Pablo Hernandez-Cerdan, Jean-Christophe Fillion-Robin");
  return moduleContributors;
}

//-----------------------------------------------------------------------------
QIcon qSlicerSRepRefinerModule::icon() const
{
  return QIcon(":/Icons/SRepRefiner.png");
}

//-----------------------------------------------------------------------------
QStringList qSlicerSRepRefinerModule::categories() const
{
  return QStringList() << "Skeleton, topology";
}

//-----------------------------------------------------------------------------
QStringList qSlicerSRepRefinerModule::dependencies() const
{
  return QStringList();
}

//-----------------------------------------------------------------------------
void qSlicerSRepRefinerModule::setup()
{
  this->Superclass::setup();
}

//-----------------------------------------------------------------------------
qSlicerAbstractModuleRepresentation* qSlicerSRepRefinerModule
::createWidgetRepresentation()
{
  return new qSlicerSRepRefinerModuleWidget;
}

//-----------------------------------------------------------------------------
vtkMRMLAbstractLogic* qSlicerSRepRefinerModule::createLogic()
{
  return vtkSlicerSRepRefinerLogic::New();
}
