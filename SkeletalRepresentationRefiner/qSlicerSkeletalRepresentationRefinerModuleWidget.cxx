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

void qSlicerSkeletalRepresentationRefinerModuleWidget::SelectOutputPath()
{
    Q_D(qSlicerSkeletalRepresentationRefinerModuleWidget);
    QString fileName = QFileDialog::getExistingDirectory(this, "Select refinement output folder");
    d->lb_outputpath->setText(fileName.toUtf8().constData());
    d->logic()->SetOutputPath(fileName.toUtf8().constData());
}

void qSlicerSkeletalRepresentationRefinerModuleWidget::StartRefinement()
{
    Q_D(qSlicerSkeletalRepresentationRefinerModuleWidget);
    double stepSize = d->sl_stepSize->value();
    double tol = d->sl_tol->value();
    int maxIter = static_cast<int>(d->sl_maxIter->value());
    
    double wtImageMatch = d->sl_wtImageMatch->value();
    double wtNormalMatch = d->sl_wtNormal->value();
    double wtSrad = d->sl_wtSrad->value();
    int interpLevel = static_cast<int>(d->sl_interp->value());
    
    d->logic()->SetWeights(wtImageMatch, wtNormalMatch, wtSrad);
    d->logic()->Refine(stepSize, tol, maxIter, interpLevel);
}

void qSlicerSkeletalRepresentationRefinerModuleWidget::StartInterpolate()
{
    Q_D(qSlicerSkeletalRepresentationRefinerModuleWidget);
    int interpolationLevel = int(d->sl_interp->value());
    std::string srepFileName = d->lb_srepPath->text().toUtf8().constData();
    d->logic()->InterpolateSrep(interpolationLevel, srepFileName);
}

void qSlicerSkeletalRepresentationRefinerModuleWidget::GenerateImage()
{
    Q_D(qSlicerSkeletalRepresentationRefinerModuleWidget);
    d->logic()->AntiAliasSignedDistanceMap(d->lb_imagePath->text().toUtf8().constData());
}

void qSlicerSkeletalRepresentationRefinerModuleWidget::TransformSrep()
{
    Q_D(qSlicerSkeletalRepresentationRefinerModuleWidget);
    std::string headerFile = d->lb_srepPath->text().toUtf8().constData();
    d->logic()->TransformSrep(headerFile);
}

//-----------------------------------------------------------------------------
void qSlicerSkeletalRepresentationRefinerModuleWidget::setup()
{
  Q_D(qSlicerSkeletalRepresentationRefinerModuleWidget);
  d->setupUi(this);
  this->Superclass::setup();
  QObject::connect(d->btn_browseImage, SIGNAL(clicked()), this, SLOT(SelectImage()));
  QObject::connect(d->btn_browseSrep, SIGNAL(clicked()), this, SLOT(SelectSrep()));
  QObject::connect(d->btn_output, SIGNAL(clicked()), this, SLOT(SelectOutputPath()));
  QObject::connect(d->btn_submit, SIGNAL(clicked()), this, SLOT(StartRefinement()));
  QObject::connect(d->btn_interp, SIGNAL(clicked()), this, SLOT(StartInterpolate()));
  //QObject::connect(d->btn_image, SIGNAL(clicked()), this, SLOT(GenerateImage()));
  //QObject::connect(d->btn_transform, SIGNAL(clicked()), this, SLOT(TransformSrep()));
}

