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

// SRepRefinement Logic includes
#include <vtkSlicerSRepRefinementLogic.h>

// SRepRefinement includes
#include "qSlicerSRepRefinementModule.h"
#include "qSlicerSRepRefinementModuleWidget.h"

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerSRepRefinementModulePrivate
{
public:
  qSlicerSRepRefinementModulePrivate();
};

//-----------------------------------------------------------------------------
// qSlicerSRepRefinementModulePrivate methods

//-----------------------------------------------------------------------------
qSlicerSRepRefinementModulePrivate::qSlicerSRepRefinementModulePrivate()
{
}

//-----------------------------------------------------------------------------
// qSlicerSRepRefinementModule methods

//-----------------------------------------------------------------------------
qSlicerSRepRefinementModule::qSlicerSRepRefinementModule(QObject* _parent)
  : Superclass(_parent)
  , d_ptr(new qSlicerSRepRefinementModulePrivate)
{
}

//-----------------------------------------------------------------------------
qSlicerSRepRefinementModule::~qSlicerSRepRefinementModule()
{
}

//-----------------------------------------------------------------------------
QString qSlicerSRepRefinementModule::helpText() const
{
  return "<p>This module is used for refining SReps using the NEWUOA algorithm to better fit to a model"
        " by minimizing an objective function consisting of</p>"
        "<ul>"
        "  <li>L<sub>0</sub>: The squared distance between the SRep implied boundary and the model boundary</li>"
        "  <li>L<sub>1</sub>: The deviation of spokes from being perpendicular to the boundary</li>"
        "  <li>L<sub>2</sub>: The geometric illegality of spokes (spokes cannot cross each other)</li>"
        "</ul>"
        "<p>Parameters</p>"
        "<ul>"
        "  <li>Input Model: the model to refine the SRep to</li>"
        "  <li>Input SRep: the SRep to be refined</li>"
        "  <li>Output SRep: the SRep object to put the refined SRep into</li>"
        "  <li>Interpolation level: how much to interpolate between spokes. The interpolated"
             " spokes are used to define the implied boundary used by the objective function."
             " Interpolated spokes are produced at 2^level times the original density.</li>"
        "  <li>Initial region size: the initial value of the newuoa trust region radius.</li>"
        "  <li>Final region size: the final value of the newuoa trust region radius."
             " Typically this is around one tenth the greatest expected change to a variable.</li>"
        "  <li>Max iterations: the maximum amount of iterations to run</li>"
        "  <li>Image match weight: the amount of weight to give to the L<sub>0</sub> piece of"
             " the objective function</li>"
        "  <li>Normal match weight: the amount of weight to give to the L<sub>1</sub> piece of"
             " the objective function</li>"
        "  <li>Geometric illegality weight: the amount of weight to give to the L<sub>2</sub> piece of"
             " the objective function</li>"
        "</ul>"
        "<p>Operations</p>"
        "<ul>"
        "  <li>Visualize interpolation: visualize the interpolated spokes</li>"
        "</ul>"
  ;
}

//-----------------------------------------------------------------------------
QString qSlicerSRepRefinementModule::acknowledgementText() const
{
  return "This work was heavily based off the original SRep work by Zhiyuan Liu, Megan Stuart, Jiyao Wang, Pablo Hernandez-Cerdan, Jean-Christophe Fillion-Robin";
}

//-----------------------------------------------------------------------------
QStringList qSlicerSRepRefinementModule::contributors() const
{
  QStringList moduleContributors;
  moduleContributors << QString("Connor Bowley (Kitware, Inc.)");
  return moduleContributors;
}

//-----------------------------------------------------------------------------
QIcon qSlicerSRepRefinementModule::icon() const
{
  return QIcon(":/Icons/SRepRefinement.png");
}

//-----------------------------------------------------------------------------
QStringList qSlicerSRepRefinementModule::categories() const
{
  return QStringList() << "SRep2";
}

//-----------------------------------------------------------------------------
QStringList qSlicerSRepRefinementModule::dependencies() const
{
  return QStringList();
}

//-----------------------------------------------------------------------------
void qSlicerSRepRefinementModule::setup()
{
  this->Superclass::setup();
}

//-----------------------------------------------------------------------------
qSlicerAbstractModuleRepresentation* qSlicerSRepRefinementModule
::createWidgetRepresentation()
{
  return new qSlicerSRepRefinementModuleWidget;
}

//-----------------------------------------------------------------------------
vtkMRMLAbstractLogic* qSlicerSRepRefinementModule::createLogic()
{
  return vtkSlicerSRepRefinementLogic::New();
}
