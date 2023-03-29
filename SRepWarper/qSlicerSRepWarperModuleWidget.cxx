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
#include "qSlicerSRepWarperModuleWidget.h"
#include "ui_qSlicerSRepWarperModuleWidget.h"
#include "vtkSlicerSRepWarperLogic.h"

#include "SRepProgressHelper.h"
#include <srepUtil.h>

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerSRepWarperModuleWidgetPrivate: public Ui_qSlicerSRepWarperModuleWidget
{
  Q_DECLARE_PUBLIC(qSlicerSRepWarperModuleWidget)
public:
  qSlicerSRepWarperModuleWidgetPrivate(qSlicerSRepWarperModuleWidget* object);
  vtkSlicerSRepWarperLogic* logic() const;
private:
  qSlicerSRepWarperModuleWidget* const q_ptr;
};

//-----------------------------------------------------------------------------
// qSlicerSRepWarperModuleWidgetPrivate methods

//-----------------------------------------------------------------------------
qSlicerSRepWarperModuleWidgetPrivate::qSlicerSRepWarperModuleWidgetPrivate(qSlicerSRepWarperModuleWidget* object)
  : q_ptr(object)
{}

//-----------------------------------------------------------------------------
vtkSlicerSRepWarperLogic* qSlicerSRepWarperModuleWidgetPrivate::logic() const
{
    Q_Q(const qSlicerSRepWarperModuleWidget);
    return vtkSlicerSRepWarperLogic::SafeDownCast(q->logic());
}

//-----------------------------------------------------------------------------
// qSlicerSRepWarperModuleWidget methods

//-----------------------------------------------------------------------------
qSlicerSRepWarperModuleWidget::qSlicerSRepWarperModuleWidget(QWidget* _parent)
  : Superclass( _parent )
  , d_ptr( new qSlicerSRepWarperModuleWidgetPrivate(this) )
{}

//-----------------------------------------------------------------------------
qSlicerSRepWarperModuleWidget::~qSlicerSRepWarperModuleWidget() = default;

//-----------------------------------------------------------------------------
void qSlicerSRepWarperModuleWidget::setup()
{
  Q_D(qSlicerSRepWarperModuleWidget);
  d->setupUi(this);
  this->Superclass::setup();
  d->progressBar->hide();

  QObject::connect(d->runButton, SIGNAL(clicked()), this, SLOT(run()));
}

//-----------------------------------------------------------------------------
void qSlicerSRepWarperModuleWidget::setMRMLScene(vtkMRMLScene* scene) {
  Q_D(qSlicerSRepWarperModuleWidget);
  d->sourceModelCbox->setMRMLScene(scene);
  d->sourceSRepCbox->setMRMLScene(scene);
  d->targetModelCbox->setMRMLScene(scene);
  d->outputSRepCbox->setMRMLScene(scene);
  Superclass::setMRMLScene(scene);
}

//-----------------------------------------------------------------------------
void qSlicerSRepWarperModuleWidget::run()
{
  Q_D(qSlicerSRepWarperModuleWidget);

  auto sourceModelNode = vtkMRMLModelNode::SafeDownCast(d->sourceModelCbox->currentNode());
  auto sourceSRepNode = vtkMRMLEllipticalSRepNode::SafeDownCast(d->sourceSRepCbox->currentNode());
  auto targetModelNode = vtkMRMLModelNode::SafeDownCast(d->targetModelCbox->currentNode());

  auto outputSRepNode = vtkMRMLEllipticalSRepNode::SafeDownCast(d->outputSRepCbox->currentNode());

  try {
  d->progressBar->show();
  const auto fin = srep::util::finally([&d](){ d->progressBar->hide(); });
  SRepProgressHelper<QProgressBar> progressManager(*(d->logic()), d->progressBar);
  d->logic()->Run(
    sourceModelNode,
    sourceSRepNode,
    targetModelNode,
    outputSRepNode);
  
  // Make the result visible
  sourceModelNode->SetDisplayVisibility(false);
  sourceSRepNode->SetDisplayVisibility(false);
  targetModelNode->SetDisplayVisibility(false);

  outputSRepNode->CreateDefaultDisplayNodes();
  outputSRepNode->SetDisplayVisibility(true);
  } catch (const std::exception& e) {
    QMessageBox::warning(this, "Error warping SRep", e.what());
  } 
}

