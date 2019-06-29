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
#include <vtkSlicerSkeletalRepresentationRefinerLogic.h>

// SkeletalRepresentationRefiner includes
#include "qSlicerSkeletalRepresentationRefinerModule.h"
#include "qSlicerSkeletalRepresentationRefinerModuleWidget.h"

//-----------------------------------------------------------------------------
#if (QT_VERSION < QT_VERSION_CHECK(5, 0, 0))
#include <QtPlugin>
Q_EXPORT_PLUGIN2(qSlicerSkeletalRepresentationRefinerModule, qSlicerSkeletalRepresentationRefinerModule);
#endif

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerSkeletalRepresentationRefinerModulePrivate
{
public:
  qSlicerSkeletalRepresentationRefinerModulePrivate();
};

//-----------------------------------------------------------------------------
// qSlicerSkeletalRepresentationRefinerModulePrivate methods

//-----------------------------------------------------------------------------
qSlicerSkeletalRepresentationRefinerModulePrivate::qSlicerSkeletalRepresentationRefinerModulePrivate()
{
}

//-----------------------------------------------------------------------------
// qSlicerSkeletalRepresentationRefinerModule methods

//-----------------------------------------------------------------------------
qSlicerSkeletalRepresentationRefinerModule::qSlicerSkeletalRepresentationRefinerModule(QObject* _parent)
  : Superclass(_parent)
  , d_ptr(new qSlicerSkeletalRepresentationRefinerModulePrivate)
{
}

//-----------------------------------------------------------------------------
qSlicerSkeletalRepresentationRefinerModule::~qSlicerSkeletalRepresentationRefinerModule()
{
}

//-----------------------------------------------------------------------------
QString qSlicerSkeletalRepresentationRefinerModule::helpText() const
{
  return "This is a loadable module that can be bundled in an extension";
}

//-----------------------------------------------------------------------------
QString qSlicerSkeletalRepresentationRefinerModule::acknowledgementText() const
{
  return "This file was originally developed by Zhiyuan Liu, and currently maintained by the SlicerSALT team.";
}

//-----------------------------------------------------------------------------
QStringList qSlicerSkeletalRepresentationRefinerModule::contributors() const
{
  QStringList moduleContributors;
  moduleContributors << QString("Zhiyuan Liu, Pablo Hernandez-Cerdan, Megan Stuart, Jiyao Wang");
  return moduleContributors;
}

//-----------------------------------------------------------------------------
QIcon qSlicerSkeletalRepresentationRefinerModule::icon() const
{
  return QIcon(":/Icons/SkeletalRepresentationRefiner.png");
}

//-----------------------------------------------------------------------------
QStringList qSlicerSkeletalRepresentationRefinerModule::categories() const
{
  return QStringList() << "Skeleton, topology";
}

//-----------------------------------------------------------------------------
QStringList qSlicerSkeletalRepresentationRefinerModule::dependencies() const
{
  return QStringList();
}

//-----------------------------------------------------------------------------
void qSlicerSkeletalRepresentationRefinerModule::setup()
{
  this->Superclass::setup();
}

//-----------------------------------------------------------------------------
qSlicerAbstractModuleRepresentation* qSlicerSkeletalRepresentationRefinerModule
::createWidgetRepresentation()
{
  return new qSlicerSkeletalRepresentationRefinerModuleWidget;
}

//-----------------------------------------------------------------------------
vtkMRMLAbstractLogic* qSlicerSkeletalRepresentationRefinerModule::createLogic()
{
  return vtkSlicerSkeletalRepresentationRefinerLogic::New();
}
