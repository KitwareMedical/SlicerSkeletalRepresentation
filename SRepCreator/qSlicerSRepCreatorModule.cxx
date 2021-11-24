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
  return "<p>This module initializes an SRep from a model in two steps.</p>\n"
         "<ol>"
         "<li>Forward Flow: creates an ellipsoid that best fits the input model"
            " and creates an SRep to fit that ellipsoid. The model is \"flowed\""
            " toward an ellipsoidal shape for a number of iterations, then a best"
            " fit ellipsoid is made on the flowed shaped.</li>"
         "<li>Backward Flow: adjusts the SRep created from the ellipsoid to fit"
            " the original model by reversing the flow transformations made during"
            " the forward flow.</li>"
         "</ol>"
         "<p>Parameters</p>"
         "<ul>"
         "  <li>Input Mesh: model to create the SRep of.</li>"
         "  <li>Max iterations: the number of iterations to run in the forward flow.</li>"
         "  <li>Step size: the size of step to take in during the forward flow.</li>"
         "  <li>Smooth amount: the amount of smoothing that should be applied to the model"
              " while flowing toward the ellipsoidal shape. The algorithm doesn't work well"
              " with sharp edges or points, so smoothing can help with that.</li>"
         "  <li># Fold Points: The number of fold (aka crest) points in the generated SRep.</li>"
         "  <li># Steps to Fold: The number of steps from the spine to outer boundary of the"
              " skeletal sheet. The point on the spine is not included in this number.</li>"
         "</ul>"
         "<p>Operations</p>"
         "<ul>"
         "  <li>Run: Run both forward and backward flows to generate an SRep.</li>"
         "</ul>"
         "<p>Advanced</p>"
         "<ul>"
         "  <li>Forward output every # iterations: Outputs the flowed mesh as a new model node"
              " every given iterations.</li>"
         "  <li>Backward output every # iterations: Outputs the backflowed SRep as a new SRep"
              " node every given iterations.</li>"
         "  <li>Output fitted ellipsoid: Outputs a model node and an SRep node for the best"
              " fitting ellipse at the end of the forward flow</li>"
         "  <li>Run forward only: Only run the forward flow.</li>"
         "  <li>Run backward only: Only run the backward flow. Will fail if the forward flow"
              " has not been run.</li>"
         "</ul>"
  ;
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
