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

// SkeletalRepresentationInitializer Logic includes
#include <vtkSlicerSkeletalRepresentationInitializerLogic.h>

// SkeletalRepresentationInitializer includes
#include "qSlicerSkeletalRepresentationInitializerModule.h"
#include "qSlicerSkeletalRepresentationInitializerModuleWidget.h"

//-----------------------------------------------------------------------------
#if (QT_VERSION < QT_VERSION_CHECK(5, 0, 0))
#include <QtPlugin>
Q_EXPORT_PLUGIN2(qSlicerSkeletalRepresentationInitializerModule, qSlicerSkeletalRepresentationInitializerModule);
#endif

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerSkeletalRepresentationInitializerModulePrivate
{
public:
  qSlicerSkeletalRepresentationInitializerModulePrivate();
};

//-----------------------------------------------------------------------------
// qSlicerSkeletalRepresentationInitializerModulePrivate methods

//-----------------------------------------------------------------------------
qSlicerSkeletalRepresentationInitializerModulePrivate::qSlicerSkeletalRepresentationInitializerModulePrivate()
{
}

//-----------------------------------------------------------------------------
// qSlicerSkeletalRepresentationInitializerModule methods

//-----------------------------------------------------------------------------
qSlicerSkeletalRepresentationInitializerModule::qSlicerSkeletalRepresentationInitializerModule(QObject* _parent)
  : Superclass(_parent)
  , d_ptr(new qSlicerSkeletalRepresentationInitializerModulePrivate)
{
}

//-----------------------------------------------------------------------------
qSlicerSkeletalRepresentationInitializerModule::~qSlicerSkeletalRepresentationInitializerModule()
{
}

//-----------------------------------------------------------------------------
QString qSlicerSkeletalRepresentationInitializerModule::helpText() const
{
    return "Given a surface, use geometric flow to initialize its skeletal model.";
}

//-----------------------------------------------------------------------------
QString qSlicerSkeletalRepresentationInitializerModule::acknowledgementText() const
{
    return "This file was originally developed by Junpyo Hong, Zhiyuan Liu, and currently maintained by the SlicerSALT team.";
}

//-----------------------------------------------------------------------------
QStringList qSlicerSkeletalRepresentationInitializerModule::contributors() const
{
  QStringList moduleContributors;
  moduleContributors << QString("Junpyo Hong, Zhiyuan Liu, Pablo Hernandez-Cerdan");
  return moduleContributors;
}

//-----------------------------------------------------------------------------
QIcon qSlicerSkeletalRepresentationInitializerModule::icon() const
{
  return QIcon(":/Icons/SkeletalRepresentationInitializer.png");
}

//-----------------------------------------------------------------------------
QStringList qSlicerSkeletalRepresentationInitializerModule::categories() const
{
  return QStringList() << "Skeleton, topology";
}

//-----------------------------------------------------------------------------
QStringList qSlicerSkeletalRepresentationInitializerModule::dependencies() const
{
  return QStringList();
}

//-----------------------------------------------------------------------------
void qSlicerSkeletalRepresentationInitializerModule::setup()
{
  this->Superclass::setup();
}

//-----------------------------------------------------------------------------
qSlicerAbstractModuleRepresentation* qSlicerSkeletalRepresentationInitializerModule
::createWidgetRepresentation()
{
  return new qSlicerSkeletalRepresentationInitializerModuleWidget;
}

//-----------------------------------------------------------------------------
vtkMRMLAbstractLogic* qSlicerSkeletalRepresentationInitializerModule::createLogic()
{
  return vtkSlicerSkeletalRepresentationInitializerLogic::New();
}
