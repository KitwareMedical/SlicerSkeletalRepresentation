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

// SlicerQt includes
#include "qSlicerSkeletalRepresentationRefinerModuleWidget.h"
#include "ui_qSlicerSkeletalRepresentationRefinerModuleWidget.h"

#include "vtkSlicerSkeletalRepresentationRefinerLogic.h"

#include <QFileDialog>
#include <QMessageBox>

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerSkeletalRepresentationRefinerModuleWidgetPrivate: public Ui_qSlicerSkeletalRepresentationRefinerModuleWidget
{
    Q_DECLARE_PUBLIC(qSlicerSkeletalRepresentationRefinerModuleWidget);
protected:
    qSlicerSkeletalRepresentationRefinerModuleWidget* const q_ptr;
public:
  qSlicerSkeletalRepresentationRefinerModuleWidgetPrivate(qSlicerSkeletalRepresentationRefinerModuleWidget &);
  ~qSlicerSkeletalRepresentationRefinerModuleWidgetPrivate();
  vtkSlicerSkeletalRepresentationRefinerLogic* logic() const;
};

//-----------------------------------------------------------------------------
// qSlicerSkeletalRepresentationRefinerModuleWidgetPrivate methods

//-----------------------------------------------------------------------------
qSlicerSkeletalRepresentationRefinerModuleWidgetPrivate::qSlicerSkeletalRepresentationRefinerModuleWidgetPrivate(qSlicerSkeletalRepresentationRefinerModuleWidget &object)
 : q_ptr(&object){
}

qSlicerSkeletalRepresentationRefinerModuleWidgetPrivate::~qSlicerSkeletalRepresentationRefinerModuleWidgetPrivate()
{

}

vtkSlicerSkeletalRepresentationRefinerLogic *qSlicerSkeletalRepresentationRefinerModuleWidgetPrivate::logic() const
{
    Q_Q(const qSlicerSkeletalRepresentationRefinerModuleWidget);
    return vtkSlicerSkeletalRepresentationRefinerLogic::SafeDownCast(q->logic());
}

//-----------------------------------------------------------------------------
// qSlicerSkeletalRepresentationRefinerModuleWidget methods

//-----------------------------------------------------------------------------
qSlicerSkeletalRepresentationRefinerModuleWidget::qSlicerSkeletalRepresentationRefinerModuleWidget(QWidget* _parent)
  : Superclass( _parent )
  , d_ptr( new qSlicerSkeletalRepresentationRefinerModuleWidgetPrivate(*this) )
{
}

//-----------------------------------------------------------------------------
qSlicerSkeletalRepresentationRefinerModuleWidget::~qSlicerSkeletalRepresentationRefinerModuleWidget()
{
}

void qSlicerSkeletalRepresentationRefinerModuleWidget::SelectImage()
{
    Q_D(qSlicerSkeletalRepresentationRefinerModuleWidget);
    QString fileName = QFileDialog::getOpenFileName(this, "Select image file");
    d->lb_imagePath->setText(fileName.toUtf8().constData());
    d->logic()->SetImageFileName(fileName.toUtf8().constData());
}

void qSlicerSkeletalRepresentationRefinerModuleWidget::SelectSrep()
{
    Q_D(qSlicerSkeletalRepresentationRefinerModuleWidget);
    QString fileName = QFileDialog::getOpenFileName(this, "Select s-rep file");
    d->lb_srepPath->setText(fileName.toUtf8().constData());
    d->logic()->SetSrepFileName(fileName.toUtf8().constData());
}

void qSlicerSkeletalRepresentationRefinerModuleWidget::Submit()
{
    Q_D(qSlicerSkeletalRepresentationRefinerModuleWidget);
    d->logic()->Refine();
}

void qSlicerSkeletalRepresentationRefinerModuleWidget::StartInterpolate()
{
    Q_D(qSlicerSkeletalRepresentationRefinerModuleWidget);
    int interpolationLevel = int(d->sl_interp->value());
    std::string srepFileName = d->lb_srepPath->text().toUtf8().constData();
    d->logic()->InterpolateSrep(interpolationLevel, srepFileName);
}

//-----------------------------------------------------------------------------
void qSlicerSkeletalRepresentationRefinerModuleWidget::setup()
{
  Q_D(qSlicerSkeletalRepresentationRefinerModuleWidget);
  d->setupUi(this);
  this->Superclass::setup();
  QObject::connect(d->btn_browseImage, SIGNAL(clicked()), this, SLOT(SelectImage()));
  QObject::connect(d->btn_browseSrep, SIGNAL(clicked()), this, SLOT(SelectSrep()));
  QObject::connect(d->btn_submit, SIGNAL(clicked()), this, SLOT(Submit()));
  QObject::connect(d->btn_interp, SIGNAL(clicked()), this, SLOT(StartInterpolate()));
}

