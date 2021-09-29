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

// Slicer includes
#include "qSlicerSRepModuleWidget.h"
#include "ui_qSlicerSRepModuleWidget.h"
#include "vtkSlicerSRepLogic.h"

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerSRepModuleWidgetPrivate: public Ui_qSlicerSRepModuleWidget
{
  Q_DECLARE_PUBLIC(qSlicerSRepModuleWidget)
public:
  qSlicerSRepModuleWidgetPrivate(qSlicerSRepModuleWidget& object);
  vtkSlicerSRepLogic* logic() const;
private:
    qSlicerSRepModuleWidget* const q_ptr;
};

//-----------------------------------------------------------------------------
// qSlicerSRepModuleWidgetPrivate methods

//-----------------------------------------------------------------------------
qSlicerSRepModuleWidgetPrivate::qSlicerSRepModuleWidgetPrivate(qSlicerSRepModuleWidget& object)
  : q_ptr(&object)
{
}

vtkSlicerSRepLogic* qSlicerSRepModuleWidgetPrivate::logic() const
{
    Q_Q(const qSlicerSRepModuleWidget);
    return vtkSlicerSRepLogic::SafeDownCast(q->logic());
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
  d->setupUi(this);
  this->Superclass::setup();
  QObject::connect(d->inputFileBrowseButton, SIGNAL(clicked()), this, SLOT(onInputFileBrowse()));
  QObject::connect(d->importButton, SIGNAL(clicked()), this, SLOT(onImport()));
}

//-----------------------------------------------------------------------------
void qSlicerSRepModuleWidget::onInputFileBrowse()
{
  Q_D(qSlicerSRepModuleWidget);
  QString fileName = QFileDialog::getOpenFileName(this, "Select input mesh");
  d->inputFileLineEdit->setText(fileName.toUtf8().constData());
}

void qSlicerSRepModuleWidget::onImport()
{
  Q_D(qSlicerSRepModuleWidget);
  QString inputFile = d->inputFileLineEdit->text();
  if (inputFile.isEmpty()) {
    QMessageBox::critical(this, QObject::tr("Error"), "Input file must not be empty");
    return;
  }

  d->logic()->ImportSRep(inputFile.toStdString());
}