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

// SRepInitializer Logic includes
#include <vtkSlicerSRepInitializerLogic.h>

// SRepInitializer includes
#include "qSlicerSRepInitializerModule.h"
#include "qSlicerSRepInitializerModuleWidget.h"

//-----------------------------------------------------------------------------
#if (QT_VERSION < QT_VERSION_CHECK(5, 0, 0))
#include <QtPlugin>
Q_EXPORT_PLUGIN2(qSlicerSRepInitializerModule, qSlicerSRepInitializerModule);
#endif

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerSRepInitializerModulePrivate
{
public:
  qSlicerSRepInitializerModulePrivate();
};

//-----------------------------------------------------------------------------
// qSlicerSRepInitializerModulePrivate methods

//-----------------------------------------------------------------------------
qSlicerSRepInitializerModulePrivate::qSlicerSRepInitializerModulePrivate()
{
}

//-----------------------------------------------------------------------------
// qSlicerSRepInitializerModule methods

//-----------------------------------------------------------------------------
qSlicerSRepInitializerModule::qSlicerSRepInitializerModule(QObject* _parent)
  : Superclass(_parent)
  , d_ptr(new qSlicerSRepInitializerModulePrivate)
{
}

//-----------------------------------------------------------------------------
qSlicerSRepInitializerModule::~qSlicerSRepInitializerModule()
{
}

//-----------------------------------------------------------------------------
QString qSlicerSRepInitializerModule::helpText() const
{
    return "Given a surface, use geometric flow to initialize its skeletal model.";
}

//-----------------------------------------------------------------------------
QString qSlicerSRepInitializerModule::acknowledgementText() const
{
    return "This file was originally developed by Zhiyuan Liu and Junpyo Hong, and currently maintained by the SlicerSALT team.";
}

//-----------------------------------------------------------------------------
QStringList qSlicerSRepInitializerModule::contributors() const
{
  QStringList moduleContributors;
  moduleContributors << QString("Zhiyuan Liu, Pablo Hernandez-Cerdan, Megan Stuart, Jiyao Wang");
  return moduleContributors;
}

//-----------------------------------------------------------------------------
QIcon qSlicerSRepInitializerModule::icon() const
{
  return QIcon(":/Icons/SRepInitializer.png");
}

//-----------------------------------------------------------------------------
QStringList qSlicerSRepInitializerModule::categories() const
{
  return QStringList() << "Skeleton, topology";
}

//-----------------------------------------------------------------------------
QStringList qSlicerSRepInitializerModule::dependencies() const
{
  return QStringList();
}

//-----------------------------------------------------------------------------
void qSlicerSRepInitializerModule::setup()
{
  this->Superclass::setup();
}

//-----------------------------------------------------------------------------
qSlicerAbstractModuleRepresentation* qSlicerSRepInitializerModule
::createWidgetRepresentation()
{
  return new qSlicerSRepInitializerModuleWidget;
}

//-----------------------------------------------------------------------------
vtkMRMLAbstractLogic* qSlicerSRepInitializerModule::createLogic()
{
  return vtkSlicerSRepInitializerLogic::New();
}
