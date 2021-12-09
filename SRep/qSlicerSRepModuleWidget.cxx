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
#include <QFileDialog>
#include <QMessageBox>
#include <QSignalBlocker>

// Slicer includes
#include "qSlicerSRepModuleWidget.h"
#include "ui_qSlicerSRepModuleWidget.h"
#include "vtkSlicerSRepLogic.h"
#include <qMRMLSubjectHierarchyModel.h>
#include <vtkMRMLScene.h>

#include <functional>

namespace {

QColor toQColor(const vtkColor3ub& color) {
  return QColor(color[0], color[1], color[2]);
}

vtkColor3ub toVTKColor3ub(const QColor& color) {
  return vtkColor3ub(color.red(), color.green(), color.blue());
}

} // namespace {}

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerSRepModuleWidgetPrivate: public Ui_qSlicerSRepModuleWidget
{
  Q_DECLARE_PUBLIC(qSlicerSRepModuleWidget)
public:
  qSlicerSRepModuleWidgetPrivate(qSlicerSRepModuleWidget& object);
  vtkSlicerSRepLogic* logic() const;

  void setupSRepUi(qSlicerWidget* widget);
  vtkWeakPointer<vtkMRMLSRepNode> activeSRepNode;
private:
  using SetColorFunc = void (vtkMRMLSRepDisplayNode::*)(const vtkColor3ub&);

  qSlicerSRepModuleWidget* const q_ptr;

  void onColorChanged(ctkColorPickerButton* button, SetColorFunc setFunc);
  void setupThicknessSlider(bool isAbsolute);
};

//-----------------------------------------------------------------------------
// qSlicerSRepModuleWidgetPrivate methods

//-----------------------------------------------------------------------------
qSlicerSRepModuleWidgetPrivate::qSlicerSRepModuleWidgetPrivate(qSlicerSRepModuleWidget& object)
  : activeSRepNode(nullptr)
  , q_ptr(&object)
{
}

//-----------------------------------------------------------------------------
vtkSlicerSRepLogic* qSlicerSRepModuleWidgetPrivate::logic() const
{
    Q_Q(const qSlicerSRepModuleWidget);
    return vtkSlicerSRepLogic::SafeDownCast(q->logic());
}

//-----------------------------------------------------------------------------
void qSlicerSRepModuleWidgetPrivate::setupThicknessSlider(bool isAbsolute) {
  if (isAbsolute) {
    this->useAbsoluteThicknessButton->setText("absolute");
    this->thicknessCTKSlider->setDecimals(2);
    this->thicknessCTKSlider->setMinimum(0.01);
    this->thicknessCTKSlider->setMaximum(10.0);
    this->thicknessCTKSlider->setSingleStep(0.01);
    this->thicknessCTKSlider->setPageStep(1);
    this->thicknessCTKSlider->setSuffix("mm");
  } else {
    this->useAbsoluteThicknessButton->setText("relative");
    this->thicknessCTKSlider->setDecimals(1);
    this->thicknessCTKSlider->setMinimum(0.1);
    this->thicknessCTKSlider->setMaximum(10);
    this->thicknessCTKSlider->setSingleStep(0.1);
    this->thicknessCTKSlider->setPageStep(2);
    this->thicknessCTKSlider->setSuffix("%");
  }
}

//-----------------------------------------------------------------------------
void qSlicerSRepModuleWidgetPrivate::setupSRepUi(qSlicerWidget* widget) {
  Q_Q(qSlicerSRepModuleWidget);
  this->setupUi(widget);

  this->activeSRepTreeView->setNodeTypes(QStringList()
    << "vtkMRMLEllipticalSRepNode"
    << "vtkMRMLRectangularGridSRepNode");
  this->activeSRepTreeView->setColumnHidden(this->activeSRepTreeView->model()->idColumn(), true);
  this->activeSRepTreeView->setColumnHidden(this->activeSRepTreeView->model()->transformColumn(), true);
  this->activeSRepTreeView->setColumnHidden(this->activeSRepTreeView->model()->descriptionColumn(), false);
  QObject::connect(this->activeSRepTreeView, SIGNAL(currentItemChanged(vtkIdType)),
    q, SLOT(onActiveSRepItemChanged(vtkIdType)));
  QObject::connect(widget, SIGNAL(mrmlSceneChanged(vtkMRMLScene*)),
    this->activeSRepTreeView, SLOT(setMRMLScene(vtkMRMLScene*)));

  QObject::connect(this->inputFileBrowseButton, SIGNAL(clicked()), q, SLOT(onInputFileBrowse()));
  QObject::connect(this->importButton, SIGNAL(clicked()), q, SLOT(onImport()));

  QObject::connect(this->outputDirectoryBrowseButton, SIGNAL(clicked()), q, SLOT(onExportDirectoryBrowse()));
  QObject::connect(this->exportButton, SIGNAL(clicked()), q, SLOT(onExport()));

  //visibility
  QObject::connect(this->visibilityCheckbox, SIGNAL(clicked()), q, SLOT(onVisibilityChanged()));

  //opacity
  QObject::connect(this->opacitySlider, SIGNAL(valueChanged(int)), q, SLOT(onOpacitySliderChanged()));
  QObject::connect(this->opacitySpinbox, SIGNAL(valueChanged(double)), q, SLOT(onOpacitySpinboxChanged()));

  //colors
  QObject::connect(this->upSpokeColorButton, &ctkColorPickerButton::colorChanged,
    [this](){ this->onColorChanged(this->upSpokeColorButton, &vtkMRMLSRepDisplayNode::SetUpSpokeColor);});
  QObject::connect(this->downSpokeColorButton, &ctkColorPickerButton::colorChanged,
    [this](){ this->onColorChanged(this->downSpokeColorButton, &vtkMRMLSRepDisplayNode::SetDownSpokeColor);});
  QObject::connect(this->crestSpokeColorButton, &ctkColorPickerButton::colorChanged,
    [this](){ this->onColorChanged(this->crestSpokeColorButton, &vtkMRMLSRepDisplayNode::SetCrestSpokeColor);});
  QObject::connect(this->skeletalSheetColorButton, &ctkColorPickerButton::colorChanged,
    [this](){ this->onColorChanged(this->skeletalSheetColorButton, &vtkMRMLSRepDisplayNode::SetSkeletalSheetColor);});
  QObject::connect(this->crestCurveColorButton, &ctkColorPickerButton::colorChanged,
    [this](){ this->onColorChanged(this->crestCurveColorButton, &vtkMRMLSRepDisplayNode::SetCrestCurveColor);});
  QObject::connect(this->skeletonToCrestConnectionColorButton, &ctkColorPickerButton::colorChanged,
    [this](){ this->onColorChanged(this->skeletonToCrestConnectionColorButton, &vtkMRMLSRepDisplayNode::SetSkeletonToCrestConnectionColor);});

  QObject::connect(widget, SIGNAL(mrmlSceneChanged(vtkMRMLScene*)),
    this->interpolationOutputNodeCbox, SLOT(setMRMLScene(vtkMRMLScene*)));
  QObject::connect(this->interpolateButton, SIGNAL(clicked()), q, SLOT(onInterpolate()));

  QObject::connect(this->thicknessCTKSlider, SIGNAL(valueChanged(double)),
    q, SLOT(onThicknessChanged()));
  QObject::connect(this->useAbsoluteThicknessButton, SIGNAL(clicked()),
    q , SLOT(onUseAbsoluteThicknessChanged()));
  this->setupThicknessSlider(false);
}

//-----------------------------------------------------------------------------
void qSlicerSRepModuleWidgetPrivate::onColorChanged(ctkColorPickerButton* button, SetColorFunc setFunc) {
  auto displayNode = this->activeSRepNode->GetSRepDisplayNode();
  if (displayNode) {
    (displayNode->*setFunc)(toVTKColor3ub(button->color()));
  }
}

//-----------------------------------------------------------------------------
// qSlicerSRepModuleWidget methods

//-----------------------------------------------------------------------------
qSlicerSRepModuleWidget::qSlicerSRepModuleWidget(QWidget* _parent)
  : Superclass( _parent )
  , d_ptr( new qSlicerSRepModuleWidgetPrivate(*this) )
{
}

//-----------------------------------------------------------------------------
qSlicerSRepModuleWidget::~qSlicerSRepModuleWidget()
{
}

//-----------------------------------------------------------------------------
void qSlicerSRepModuleWidget::setup()
{
  Q_D(qSlicerSRepModuleWidget);
  d->setupSRepUi(this);
  this->Superclass::setup();
}

//-----------------------------------------------------------------------------
void qSlicerSRepModuleWidget::onInputFileBrowse()
{
  Q_D(qSlicerSRepModuleWidget);
  QString selectedFilter = tr("XML (*.xml)");
  QString fileName = QFileDialog::getOpenFileName(
    this,
    "Select input mesh",
    QString(),
    tr("All files (*.*);;XML (*.xml)"),
    &selectedFilter
  );
  d->inputFileLineEdit->setText(fileName);
}

//-----------------------------------------------------------------------------
void qSlicerSRepModuleWidget::onImport()
{
  Q_D(qSlicerSRepModuleWidget);
  const QString inputFile = d->inputFileLineEdit->text();
  if (inputFile.isEmpty()) {
    QMessageBox::critical(this, QObject::tr("Error"), "Input file must not be empty");
    return;
  }

  d->logic()->ImportRectangularGridSRepFromXML(inputFile.toStdString());
}

//-----------------------------------------------------------------------------
void qSlicerSRepModuleWidget::onExportDirectoryBrowse() {
  Q_D(qSlicerSRepModuleWidget);
  QString dir = QFileDialog::getExistingDirectory(this, tr("Export Directory"),
                                                QString(),
                                                QFileDialog::ShowDirsOnly
                                                | QFileDialog::DontResolveSymlinks);

  d->outputDirectoryLineEdit->setText(dir);
}

//-----------------------------------------------------------------------------
void qSlicerSRepModuleWidget::onExport() {
  Q_D(qSlicerSRepModuleWidget);
  const auto exportDirectoryPath = d->outputDirectoryLineEdit->text();
  const auto exportBaseName = d->outputBaseNameLineEdit->text();
  if (exportDirectoryPath.isEmpty()) {
    QMessageBox::critical(this, QObject::tr("Error"), "Export directory must not be empty");
    return;
  }
  if (exportBaseName.isEmpty()) {
    QMessageBox::critical(this, QObject::tr("Error"), "Export base name must not be empty");
    return;
  }

  const QDir exportDirectory(exportDirectoryPath);
  const auto headerFile = exportDirectory.filePath(exportBaseName + QString("-header.xml"));
  const auto upFile = exportDirectory.filePath(exportBaseName + QString("-up-spokes.vtp"));
  const auto downFile = exportDirectory.filePath(exportBaseName + QString("-down-spokes.vtp"));
  const auto crestFile = exportDirectory.filePath(exportBaseName + QString("-crest-spokes.vtp"));

  const auto success = d->logic()->ExportRectangularGridSRepToXML(d->activeSRepNode,
                                              headerFile.toStdString(),
                                              upFile.toStdString(),
                                              downFile.toStdString(),
                                              crestFile.toStdString());
  if (success) {
    QMessageBox::information(this, QObject::tr("Export Successful"), "Exported SRep successfully");
  } else {
    QMessageBox::critical(this, QObject::tr("Export Failed"), "Failed to export SRep");
  }
}

//-----------------------------------------------------------------------------
void qSlicerSRepModuleWidget::onActiveSRepItemChanged(vtkIdType) {
  if (!this->isEntered()) {
    // ignore any changes if the GUI is not shown
    return;
  }
  Q_D(qSlicerSRepModuleWidget);
  this->onActiveSRepMRMLNodeChanged(d->activeSRepTreeView->currentNode());
}

//-----------------------------------------------------------------------------
void qSlicerSRepModuleWidget::onActiveSRepMRMLNodeChanged(vtkMRMLNode *node)
{
  if (!this->isEntered()) {
    // ignore any changes if the GUI is not shown
    return;
  }
  Q_D(qSlicerSRepModuleWidget);
  vtkMRMLSRepNode* srepNode = vtkMRMLSRepNode::SafeDownCast(node);

  this->setMRMLSRepNode(srepNode);
}

//-----------------------------------------------------------------------------
void qSlicerSRepModuleWidget::setMRMLSRepNode(vtkMRMLSRepNode* srepNode, bool forceReconnect) {
  if (!this->mrmlScene()) {
    srepNode = nullptr;
  }
  Q_D(qSlicerSRepModuleWidget);
  if (srepNode == d->activeSRepNode && !forceReconnect) {
    // no change
    return;
  }
  qvtkReconnect(d->activeSRepNode, srepNode, vtkCommand::ModifiedEvent,
    this, SLOT(onActiveSRepNodeModifiedEvent()));
  qvtkReconnect(d->activeSRepNode, srepNode, vtkMRMLDisplayableNode::DisplayModifiedEvent,
    this, SLOT(onActiveSRepNodeDisplayModifiedEvent()));
  d->activeSRepNode = srepNode;
  this->updateWidgetFromMRML();
}


//-----------------------------------------------------------------------------
void qSlicerSRepModuleWidget::onVisibilityChanged() {
  Q_D(qSlicerSRepModuleWidget);
  if (!d->activeSRepNode) {
    return;
  }

  auto displayNode = d->activeSRepNode->GetSRepDisplayNode();
  if (displayNode) {
    displayNode->SetVisibility(d->visibilityCheckbox->isChecked());
  }
}

//-----------------------------------------------------------------------------
void qSlicerSRepModuleWidget::onOpacitySliderChanged() {
  Q_D(qSlicerSRepModuleWidget);
  d->opacitySpinbox->setValue(d->opacitySlider->value() / 100.);
}

//-----------------------------------------------------------------------------
void qSlicerSRepModuleWidget::onOpacitySpinboxChanged() {
  Q_D(qSlicerSRepModuleWidget);
  d->opacitySlider->setValue(static_cast<int>(d->opacitySpinbox->value() * 100));
  auto displayNode = d->activeSRepNode->GetSRepDisplayNode();
  if (displayNode) {
    displayNode->SetOpacity(d->opacitySpinbox->value());
  }
}

//-----------------------------------------------------------------------------
void qSlicerSRepModuleWidget::updateWidgetFromMRML() {
  Q_D(qSlicerSRepModuleWidget);

  {
    QSignalBlocker block(d->activeSRepTreeView);
    d->activeSRepTreeView->setCurrentNode(d->activeSRepNode);
  }

  const bool haveActiveSRepNode = static_cast<bool>(d->activeSRepNode);
  d->displayContainer->setEnabled(haveActiveSRepNode);
  d->exportContainer->setEnabled(haveActiveSRepNode);
  d->informationContainer->setEnabled(haveActiveSRepNode);
  d->interpolationContainer->setEnabled(haveActiveSRepNode && nullptr != vtkMRMLEllipticalSRepNode::SafeDownCast(d->activeSRepNode));

  if (haveActiveSRepNode) {
    auto displayNode = d->activeSRepNode->GetSRepDisplayNode();
    if (displayNode) {
      QSignalBlocker block1(d->thicknessCTKSlider);
      QSignalBlocker block2(d->useAbsoluteThicknessButton);

      d->visibilityCheckbox->setChecked(displayNode->GetVisibility());
      d->opacitySlider->setValue(static_cast<int>(displayNode->GetOpacity() * 100));
      d->opacitySpinbox->setValue(displayNode->GetOpacity());
      if (displayNode->GetUseAbsoluteThickness()) {
        if (!d->useAbsoluteThicknessButton->isChecked()) {
          d->useAbsoluteThicknessButton->click();
        }
        d->setupThicknessSlider(true);
        d->thicknessCTKSlider->setValue(displayNode->GetAbsoluteThickness());
      } else {
        if (d->useAbsoluteThicknessButton->isChecked()) {
          d->useAbsoluteThicknessButton->click();
        }
        d->setupThicknessSlider(false);
        d->thicknessCTKSlider->setValue(displayNode->GetRelativeThickness() * 100.);
      }

      d->upSpokeColorButton->setColor(toQColor(displayNode->GetUpSpokeColor()));
      d->downSpokeColorButton->setColor(toQColor(displayNode->GetDownSpokeColor()));
      d->crestSpokeColorButton->setColor(toQColor(displayNode->GetCrestSpokeColor()));
      d->crestCurveColorButton->setColor(toQColor(displayNode->GetCrestCurveColor()));
      d->skeletalSheetColorButton->setColor(toQColor(displayNode->GetSkeletalSheetColor()));
      d->skeletonToCrestConnectionColorButton->setColor(toQColor(displayNode->GetSkeletonToCrestConnectionColor()));
    }

    const auto srep = d->activeSRepNode->GetSRep();
    if (srep) {
      d->numUpSpokesLineEdit->setText(QString::number(srep->GetUpSpokes().GetNumberOfSpokes()));
      d->numDownSpokesLineEdit->setText(QString::number(srep->GetDownSpokes().GetNumberOfSpokes()));
      d->numCrestSpokesLineEdit->setText(QString::number(srep->GetCrestSpokes().GetNumberOfSpokes()));
    } else {
      d->numUpSpokesLineEdit->setText("-");
      d->numDownSpokesLineEdit->setText("-");
      d->numCrestSpokesLineEdit->setText("-");
    }
  }
}

//-----------------------------------------------------------------------------
void qSlicerSRepModuleWidget::onMRMLSceneEndBatchProcessEvent()
{
  Q_D(qSlicerSRepModuleWidget);
  if (!this->mrmlScene()) {
    return;
  }
  // d->setMRMLSRepNodeFromSelectionNode(); //TODO??
  // force update (clear GUI if no node is selected anymore)
  this->updateWidgetFromMRML();
}

//-----------------------------------------------------------------------------
void qSlicerSRepModuleWidget::enter() {
  this->Superclass::enter();

  Q_D(qSlicerSRepModuleWidget);

  // set up mrml scene observations so that the GUI gets updated
  this->qvtkConnect(this->mrmlScene(), vtkMRMLScene::NodeAddedEvent,
                    this, SLOT(onNodeAddedEvent(vtkObject*, vtkObject*)));
  // this->qvtkConnect(this->mrmlScene(), vtkMRMLScene::NodeRemovedEvent,
  //                   this, SLOT(onNodeRemovedEvent(vtkObject*, vtkObject*)));
  this->qvtkConnect(this->mrmlScene(), vtkMRMLScene::EndImportEvent,
                    this, SLOT(onMRMLSceneEndImportEvent()));
  this->qvtkConnect(this->mrmlScene(), vtkMRMLScene::EndBatchProcessEvent,
                    this, SLOT(onMRMLSceneEndBatchProcessEvent()));
  this->qvtkConnect(this->mrmlScene(), vtkMRMLScene::EndCloseEvent,
                    this, SLOT(onMRMLSceneEndCloseEvent()));
  this->qvtkConnect(this->mrmlScene(), vtkMRMLScene::EndRestoreEvent,
                    this, SLOT(onMRMLSceneEndRestoreEvent()));

  this->setMRMLSRepNode(vtkMRMLSRepNode::SafeDownCast(d->activeSRepTreeView->currentNode()), true);
}
//-----------------------------------------------------------------------------
void qSlicerSRepModuleWidget::exit() {
  this->Superclass::exit();

  // remove mrml scene observations, don't need to update the GUI while the
  // module is not showing
  this->qvtkDisconnectAll();
}

//-----------------------------------------------------------------------------
void qSlicerSRepModuleWidget::onMRMLSceneEndCloseEvent() {
  if (!this->mrmlScene() || this->mrmlScene()->IsBatchProcessing())
    {
    return;
    }
  this->setMRMLSRepNode(nullptr);
  this->updateWidgetFromMRML();
}

//-----------------------------------------------------------------------------
void qSlicerSRepModuleWidget::onMRMLSceneEndImportEvent() {
  this->updateWidgetFromMRML();
}
//-----------------------------------------------------------------------------
void qSlicerSRepModuleWidget::onMRMLSceneEndRestoreEvent() {
  this->updateWidgetFromMRML();
}
//-----------------------------------------------------------------------------
void qSlicerSRepModuleWidget::onActiveSRepNodeModifiedEvent() {
  this->updateWidgetFromMRML();
}
//-----------------------------------------------------------------------------
void qSlicerSRepModuleWidget::onActiveSRepNodeDisplayModifiedEvent() {
  this->updateWidgetFromMRML();
}

//-----------------------------------------------------------------------------
void qSlicerSRepModuleWidget::onNodeAddedEvent(vtkObject*, vtkObject* node) {
  if (!this->mrmlScene() || this->mrmlScene()->IsBatchProcessing()) {
    return;
  }

  Q_D(qSlicerSRepModuleWidget);
  vtkMRMLSRepNode* srepNode = vtkMRMLSRepNode::SafeDownCast(node);
  if (srepNode) {
    // make it active
    d->activeSRepTreeView->setCurrentNode(srepNode);
  }
}

//-----------------------------------------------------------------------------
bool qSlicerSRepModuleWidget::setEditedNode(vtkMRMLNode* node, QString role, QString context) {
  Q_D(qSlicerSRepModuleWidget);
  Q_UNUSED(role);
  Q_UNUSED(context);

  {
    auto srepNode = vtkMRMLSRepNode::SafeDownCast(node);
    if (srepNode) {
      d->activeSRepTreeView->setCurrentNode(srepNode);
      return true;
    }
  }
  {
    auto srepDisplayNode = vtkMRMLSRepDisplayNode::SafeDownCast(node);
    if (srepDisplayNode) {
      auto srepNode = vtkMRMLSRepNode::SafeDownCast(srepDisplayNode->GetDisplayableNode());
      if (srepNode) {
        d->activeSRepTreeView->setCurrentNode(srepNode);
        return true;
      } else {
        return false;
      }
    }
  }
  return false;
}

//-----------------------------------------------------------------------------
double qSlicerSRepModuleWidget::nodeEditable(vtkMRMLNode* node) {
  if (vtkMRMLSRepNode::SafeDownCast(node)
    || vtkMRMLSRepDisplayNode::SafeDownCast(node))
  {
    return 0.5;
  }
  return 0.0;
}

//-----------------------------------------------------------------------------
void qSlicerSRepModuleWidget::onInterpolate() {
  Q_D(qSlicerSRepModuleWidget);

  const size_t interpolationLevel = std::lround(d->interpolationLevelCTKSlider->value());
  auto srepNode = vtkMRMLEllipticalSRepNode::SafeDownCast(d->activeSRepTreeView->currentNode());
  auto destination = vtkMRMLEllipticalSRepNode::SafeDownCast(d->interpolationOutputNodeCbox->currentNode());
  if (!destination) {
    QMessageBox::warning(this, "Error", "No SRep node selected to put interpolated SRep into.");
    return;
  }

  d->logic()->InterpolateSRep(srepNode, interpolationLevel, destination);
}

//-----------------------------------------------------------------------------
void qSlicerSRepModuleWidget::onThicknessChanged() {
  Q_D(qSlicerSRepModuleWidget);

  if (d->activeSRepNode) {
    auto displayNode = d->activeSRepNode->GetSRepDisplayNode();
    if (displayNode) {
      if (displayNode->GetUseAbsoluteThickness()) {
        displayNode->SetAbsoluteThickness(d->thicknessCTKSlider->value());
      } else {
        displayNode->SetRelativeThickness(d->thicknessCTKSlider->value() / 100.);
      }
    }
  }
}

//-----------------------------------------------------------------------------
void qSlicerSRepModuleWidget::onUseAbsoluteThicknessChanged() {
  Q_D(qSlicerSRepModuleWidget);

  if (d->activeSRepNode) {
    auto displayNode = d->activeSRepNode->GetSRepDisplayNode();
    if (displayNode) {
      displayNode->SetUseAbsoluteThickness(d->useAbsoluteThicknessButton->isChecked());
    }
  }
}
