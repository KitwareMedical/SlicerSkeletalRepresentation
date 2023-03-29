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

// SRepWarper Logic includes
#include <vtkSlicerSRepWarperLogic.h>

// SRepWarper includes
#include "qSlicerSRepWarperModule.h"
#include "qSlicerSRepWarperModuleWidget.h"

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerSRepWarperModulePrivate
{
public:
  qSlicerSRepWarperModulePrivate();
};

//-----------------------------------------------------------------------------
// qSlicerSRepWarperModulePrivate methods

//-----------------------------------------------------------------------------
qSlicerSRepWarperModulePrivate::qSlicerSRepWarperModulePrivate()
{
}

//-----------------------------------------------------------------------------
// qSlicerSRepWarperModule methods

//-----------------------------------------------------------------------------
qSlicerSRepWarperModule::qSlicerSRepWarperModule(QObject* _parent)
  : Superclass(_parent)
  , d_ptr(new qSlicerSRepWarperModulePrivate)
{
}

//-----------------------------------------------------------------------------
qSlicerSRepWarperModule::~qSlicerSRepWarperModule()
{
}

//-----------------------------------------------------------------------------
QString qSlicerSRepWarperModule::helpText() const
{
  return "<p>This module is used for warping an SRep fit to a model to another corresponding model using TPS.</p>"
        "<p>Parameters</p>"
        "<ul>"
        "  <li>Source Model: the model the SRep is fit to</li>"
        "  <li>Source SRep: the SRep to be warped</li>"
        "  <li>Target Model: the model to warp the SRep to</li>"
        "  <li>Output SRep: the warped SRep </li>"
        "</ul>"
  ;
}

//-----------------------------------------------------------------------------
QString qSlicerSRepWarperModule::acknowledgementText() const
{
  return "This work was heavily based off the original SRep work by Zhiyuan Liu, Megan Stuart, Jiyao Wang, Pablo Hernandez-Cerdan, Jean-Christophe Fillion-Robin";
}

//-----------------------------------------------------------------------------
QStringList qSlicerSRepWarperModule::contributors() const
{
  QStringList moduleContributors;
  moduleContributors << QString("Jared Vicory (Kitware, Inc.)");
  return moduleContributors;
}

//-----------------------------------------------------------------------------
QIcon qSlicerSRepWarperModule::icon() const
{
  return QIcon(":/Icons/SRepWarper.png");
}

//-----------------------------------------------------------------------------
QStringList qSlicerSRepWarperModule::categories() const
{
  return QStringList() << "Skeleton, topology";
}

//-----------------------------------------------------------------------------
QStringList qSlicerSRepWarperModule::dependencies() const
{
  return QStringList();
}

//-----------------------------------------------------------------------------
void qSlicerSRepWarperModule::setup()
{
  this->Superclass::setup();
}

//-----------------------------------------------------------------------------
qSlicerAbstractModuleRepresentation* qSlicerSRepWarperModule
::createWidgetRepresentation()
{
  return new qSlicerSRepWarperModuleWidget;
}

//-----------------------------------------------------------------------------
vtkMRMLAbstractLogic* qSlicerSRepWarperModule::createLogic()
{
  return vtkSlicerSRepWarperLogic::New();
}
