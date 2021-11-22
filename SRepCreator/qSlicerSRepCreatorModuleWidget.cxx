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
#include <QtWidgets/QLabel>
#include "qMRMLNodeComboBox.h"


// Slicer includes
#include "qSlicerSRepCreatorModuleWidget.h"
#include "ui_qSlicerSRepCreatorModuleWidget.h"
#include "vtkSlicerSRepCreatorLogic.h"

#include <srep/RectangularGridSRep.h>

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerSRepCreatorModuleWidgetPrivate: public Ui_qSlicerSRepCreatorModuleWidget
{
  Q_DECLARE_PUBLIC(qSlicerSRepCreatorModuleWidget)
public:
  qSlicerSRepCreatorModuleWidgetPrivate(qSlicerSRepCreatorModuleWidget* object);
  vtkSlicerSRepCreatorLogic* logic() const;
private:
  qSlicerSRepCreatorModuleWidget* const q_ptr;
};

//-----------------------------------------------------------------------------
// qSlicerSRepCreatorModuleWidgetPrivate methods

//-----------------------------------------------------------------------------
qSlicerSRepCreatorModuleWidgetPrivate::qSlicerSRepCreatorModuleWidgetPrivate(qSlicerSRepCreatorModuleWidget* object)
  : q_ptr(object)
{}

//-----------------------------------------------------------------------------
vtkSlicerSRepCreatorLogic* qSlicerSRepCreatorModuleWidgetPrivate::logic() const
{
    Q_Q(const qSlicerSRepCreatorModuleWidget);
    return vtkSlicerSRepCreatorLogic::SafeDownCast(q->logic());
}

//-----------------------------------------------------------------------------
// qSlicerSRepCreatorModuleWidget methods

//-----------------------------------------------------------------------------
qSlicerSRepCreatorModuleWidget::qSlicerSRepCreatorModuleWidget(QWidget* _parent)
  : Superclass( _parent )
  , d_ptr( new qSlicerSRepCreatorModuleWidgetPrivate(this) )
{
}

//-----------------------------------------------------------------------------
qSlicerSRepCreatorModuleWidget::~qSlicerSRepCreatorModuleWidget()
{
}

//-----------------------------------------------------------------------------
void qSlicerSRepCreatorModuleWidget::setMRMLScene(vtkMRMLScene* scene) {
  Q_D(qSlicerSRepCreatorModuleWidget);
  d->inputModelComboBox->setMRMLScene(scene);
}

//-----------------------------------------------------------------------------
void qSlicerSRepCreatorModuleWidget::setup()
{
  Q_D(qSlicerSRepCreatorModuleWidget);
  d->setupUi(this);
  this->Superclass::setup();

  QObject::connect(d->runForwardButton, &QPushButton::clicked,
    [this](){this->onRunForward();});
  QObject::connect(d->runBackwardButton, &QPushButton::clicked,
    [this](){this->onRunBackward();});
  QObject::connect(d->runButton, &QPushButton::clicked,
    [this](){this->onRun();});
  QObject::connect(d->numFoldPointsCTKSlider, &ctkSliderWidget::valueChanged,
    [this](){this->onNumFoldPointsValueChanged();});
}

//-----------------------------------------------------------------------------
void qSlicerSRepCreatorModuleWidget::onNumFoldPointsValueChanged() {
  Q_D(qSlicerSRepCreatorModuleWidget);
  //make it even
  d->numFoldPointsCTKSlider->setValue(d->numFoldPointsCTKSlider->value() - (std::lround(d->numFoldPointsCTKSlider->value()) % 2));
}

//-----------------------------------------------------------------------------
void qSlicerSRepCreatorModuleWidget::onRunForward() {
  Q_D(qSlicerSRepCreatorModuleWidget);
  auto model = vtkMRMLModelNode::SafeDownCast(d->inputModelComboBox->currentNode());
  if (!model) {
    QMessageBox::warning(this, "Error", "No model selected to make an SRep from.");
    return;
  }
  const auto numFoldPoints = std::lround(d->numFoldPointsCTKSlider->value());
  const auto numStepsToFold = std::lround(d->numStepsToFoldCTKSlider->value());
  const auto dt = d->stepSizeCTKSlider->value();
  const auto smoothAmount = d->smoothAmountCTKSlider->value();
  const auto maxIterations = std::lround(d->maxIterationsCTKSlider->value());
  const auto outputEllipsoidModel = d->outputFittedEllipsoidCheckbox->isChecked();
  const auto outputEveryNumIterations = std::lround(d->forwardOutputCTKSlider->value());

  auto srepNode = d->logic()->RunForward(model, numFoldPoints, numStepsToFold, dt, smoothAmount, maxIterations,
    outputEllipsoidModel, outputEveryNumIterations);
  if (!srepNode) {
    QMessageBox::warning(this, "Error creating SRep", "SRep node from forward flow unable to be created. Try checking the error log for more details.");
  }
}

//-----------------------------------------------------------------------------
void qSlicerSRepCreatorModuleWidget::onRunBackward() {
  Q_D(qSlicerSRepCreatorModuleWidget);
  const auto outputEveryNumIterations = std::lround(d->backwardOutputCTKSlider->value());
  auto srepNode = d->logic()->RunBackward(outputEveryNumIterations);
  if (!srepNode) {
    QMessageBox::warning(this, "Error creating SRep", "SRep node from backward flow unable to be created. Try checking the error log for more details.");
  }
}

//-----------------------------------------------------------------------------
void qSlicerSRepCreatorModuleWidget::onRun() {
  Q_D(qSlicerSRepCreatorModuleWidget);
  auto model = vtkMRMLModelNode::SafeDownCast(d->inputModelComboBox->currentNode());
  if (!model) {
    QMessageBox::warning(this, "Error", "No model selected to make an SRep from.");
    return;
  }
  const auto numFoldPoints = std::lround(d->numFoldPointsCTKSlider->value());
  const auto numStepsToFold = std::lround(d->numStepsToFoldCTKSlider->value());
  const auto dt = d->stepSizeCTKSlider->value();
  const auto smoothAmount = d->smoothAmountCTKSlider->value();
  const auto maxIterations = std::lround(d->maxIterationsCTKSlider->value());
  const auto outputEllipsoidModel = d->outputFittedEllipsoidCheckbox->isChecked();
  const auto forwardOutputEveryNumIterations = std::lround(d->forwardOutputCTKSlider->value());
  const auto backwardOutputEveryNumIterations = std::lround(d->backwardOutputCTKSlider->value());

  auto srepNode = d->logic()->Run(model, numFoldPoints, numStepsToFold, dt, smoothAmount, maxIterations,
    outputEllipsoidModel, forwardOutputEveryNumIterations, backwardOutputEveryNumIterations);
  if (!srepNode) {
    QMessageBox::warning(this, "Error creating SRep", "SRep node unable to be created. Try checking the error log for more details.");
  }
}
