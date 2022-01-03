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

// Qt includes
#include <QDebug>
#include <QMessageBox>

// Slicer includes
#include "qSlicerSRepRefinementModuleWidget.h"
#include "ui_qSlicerSRepRefinementModuleWidget.h"
#include "vtkSlicerSRepRefinementLogic.h"

#include "SRepProgressHelper.h"

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerSRepRefinementModuleWidgetPrivate: public Ui_qSlicerSRepRefinementModuleWidget
{
  Q_DECLARE_PUBLIC(qSlicerSRepRefinementModuleWidget)
public:
  qSlicerSRepRefinementModuleWidgetPrivate(qSlicerSRepRefinementModuleWidget* object);
  vtkSlicerSRepRefinementLogic* logic() const;
private:
  qSlicerSRepRefinementModuleWidget* const q_ptr;
};

//-----------------------------------------------------------------------------
// qSlicerSRepRefinementModuleWidgetPrivate methods

//-----------------------------------------------------------------------------
qSlicerSRepRefinementModuleWidgetPrivate::qSlicerSRepRefinementModuleWidgetPrivate(qSlicerSRepRefinementModuleWidget* object)
  : q_ptr(object)
{}

//-----------------------------------------------------------------------------
vtkSlicerSRepRefinementLogic* qSlicerSRepRefinementModuleWidgetPrivate::logic() const
{
    Q_Q(const qSlicerSRepRefinementModuleWidget);
    return vtkSlicerSRepRefinementLogic::SafeDownCast(q->logic());
}

//-----------------------------------------------------------------------------
// qSlicerSRepRefinementModuleWidget methods

//-----------------------------------------------------------------------------
qSlicerSRepRefinementModuleWidget::qSlicerSRepRefinementModuleWidget(QWidget* _parent)
  : Superclass( _parent )
  , d_ptr( new qSlicerSRepRefinementModuleWidgetPrivate(this) )
{}

//-----------------------------------------------------------------------------
qSlicerSRepRefinementModuleWidget::~qSlicerSRepRefinementModuleWidget() = default;

//-----------------------------------------------------------------------------
void qSlicerSRepRefinementModuleWidget::setup()
{
  Q_D(qSlicerSRepRefinementModuleWidget);
  d->setupUi(this);
  this->Superclass::setup();
  d->progressBar->hide();

  QObject::connect(d->refineButton, SIGNAL(clicked()), this, SLOT(refine()));
}

//-----------------------------------------------------------------------------
void qSlicerSRepRefinementModuleWidget::setMRMLScene(vtkMRMLScene* scene) {
  Q_D(qSlicerSRepRefinementModuleWidget);
  d->inputModelCbox->setMRMLScene(scene);
  d->inputSRepCbox->setMRMLScene(scene);
  d->outputSRepCbox->setMRMLScene(scene);
}

//-----------------------------------------------------------------------------
void qSlicerSRepRefinementModuleWidget::refine()
{
  Q_D(qSlicerSRepRefinementModuleWidget);

  auto model = vtkMRMLModelNode::SafeDownCast(d->inputModelCbox->currentNode());
  auto inputSRep = vtkMRMLEllipticalSRepNode::SafeDownCast(d->inputSRepCbox->currentNode());
  auto outputSRep = vtkMRMLEllipticalSRepNode::SafeDownCast(d->outputSRepCbox->currentNode());

  const auto interpolationLevel = std::lround(d->interpolationLevelCTKSlider->value());
  const auto initialRegionSize = d->initialRegionSizeCTKSlider->value();
  const auto finalRegionSize = d->finalRegionSizeCTKSlider->value();
  const auto maxIterations = std::lround(d->maxIterationsCTKSlider->value());
  const auto imageMatchWeight = d->imageMatchWeightCTKSlider->value();
  const auto normalMatchWeight = d->normalMatchWeightCTKSlider->value();
  const auto geometricIllegalityWeight = d->geometricIllegalityWeightCTKSlider->value();

  try {
  d->progressBar->show();
  const auto fin = srep::util::finally([&d](){ d->progressBar->hide(); });
  SRepProgressHelper<QProgressBar> progressManager(*(d->logic()), d->progressBar);
  d->logic()->Run(
    model,
    inputSRep,
    initialRegionSize,
    finalRegionSize,
    maxIterations,
    interpolationLevel,
    imageMatchWeight,
    normalMatchWeight,
    geometricIllegalityWeight,
    outputSRep);
  } catch (const std::exception& e) {
    QMessageBox::warning(this, "Error refining SRep", e.what());
  } 
}

