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
#include "qSlicerSRepRefinerModuleWidget.h"
#include "ui_qSlicerSRepRefinerModuleWidget.h"

#include "vtkSlicerSRepRefinerLogic.h"

#include <QFileDialog>
#include <QMessageBox>

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerSRepRefinerModuleWidgetPrivate: public Ui_qSlicerSRepRefinerModuleWidget
{
    Q_DECLARE_PUBLIC(qSlicerSRepRefinerModuleWidget)
protected:
    qSlicerSRepRefinerModuleWidget* const q_ptr;
public:
  qSlicerSRepRefinerModuleWidgetPrivate(qSlicerSRepRefinerModuleWidget &);
  ~qSlicerSRepRefinerModuleWidgetPrivate();
  vtkSlicerSRepRefinerLogic* logic() const;
};

//-----------------------------------------------------------------------------
// qSlicerSRepRefinerModuleWidgetPrivate methods

//-----------------------------------------------------------------------------
qSlicerSRepRefinerModuleWidgetPrivate::qSlicerSRepRefinerModuleWidgetPrivate(qSlicerSRepRefinerModuleWidget &object)
 : q_ptr(&object){
}

qSlicerSRepRefinerModuleWidgetPrivate::~qSlicerSRepRefinerModuleWidgetPrivate()
{

}

vtkSlicerSRepRefinerLogic *qSlicerSRepRefinerModuleWidgetPrivate::logic() const
{
    Q_Q(const qSlicerSRepRefinerModuleWidget);
    return vtkSlicerSRepRefinerLogic::SafeDownCast(q->logic());
}

//-----------------------------------------------------------------------------
// qSlicerSRepRefinerModuleWidget methods

//-----------------------------------------------------------------------------
qSlicerSRepRefinerModuleWidget::qSlicerSRepRefinerModuleWidget(QWidget* _parent)
  : Superclass( _parent )
  , d_ptr( new qSlicerSRepRefinerModuleWidgetPrivate(*this) )
{
}

//-----------------------------------------------------------------------------
qSlicerSRepRefinerModuleWidget::~qSlicerSRepRefinerModuleWidget()
{
}

void qSlicerSRepRefinerModuleWidget::SelectImage()
{
    Q_D(qSlicerSRepRefinerModuleWidget);
    QString fileName = QFileDialog::getOpenFileName(this, "Select image file");
    d->lb_imagePath->setText(fileName.toUtf8().constData());
    d->logic()->SetImageFileName(fileName.toUtf8().constData());
}

void qSlicerSRepRefinerModuleWidget::SelectSrep()
{
    Q_D(qSlicerSRepRefinerModuleWidget);
    QString fileName = QFileDialog::getOpenFileName(this, "Select s-rep file");
    d->lb_srepPath->setText(fileName.toUtf8().constData());
    d->logic()->SetSrepFileName(fileName.toUtf8().constData());
}

void qSlicerSRepRefinerModuleWidget::SelectOutputPath()
{
    Q_D(qSlicerSRepRefinerModuleWidget);
    QString fileName = QFileDialog::getExistingDirectory(this, "Select refinement output folder");
    d->lb_outputpath->setText(fileName.toUtf8().constData());
    d->logic()->SetOutputPath(fileName.toUtf8().constData());
}

void qSlicerSRepRefinerModuleWidget::StartRefinement()
{
    Q_D(qSlicerSRepRefinerModuleWidget);
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

void qSlicerSRepRefinerModuleWidget::StartInterpolate()
{
    Q_D(qSlicerSRepRefinerModuleWidget);
    int interpolationLevel = int(d->sl_interp->value());
    std::string srepFileName = d->lb_srepPath->text().toUtf8().constData();
    d->logic()->InterpolateSrep(interpolationLevel, srepFileName);
}

void qSlicerSRepRefinerModuleWidget::GenerateImage()
{
    Q_D(qSlicerSRepRefinerModuleWidget);
    d->logic()->AntiAliasSignedDistanceMap(d->lb_imagePath->text().toUtf8().constData());
}

void qSlicerSRepRefinerModuleWidget::TransformSrep()
{
    Q_D(qSlicerSRepRefinerModuleWidget);
    std::string headerFile = d->lb_srepPath->text().toUtf8().constData();
    d->logic()->TransformSrep(headerFile);
}

void qSlicerSRepRefinerModuleWidget::showImpliedBoundary()
{
    Q_D(qSlicerSRepRefinerModuleWidget);
    int interpolationLevel = int(d->sl_interp->value());
    std::string srepFileName = d->lb_srepPath->text().toUtf8().constData();
    d->logic()->ShowImpliedBoundary(interpolationLevel, srepFileName);
}

//-----------------------------------------------------------------------------
void qSlicerSRepRefinerModuleWidget::setup()
{
  Q_D(qSlicerSRepRefinerModuleWidget);
  d->setupUi(this);
  this->Superclass::setup();
  QObject::connect(d->btn_browseImage, SIGNAL(clicked()), this, SLOT(SelectImage()));
  QObject::connect(d->btn_browseSrep, SIGNAL(clicked()), this, SLOT(SelectSrep()));
  QObject::connect(d->btn_output, SIGNAL(clicked()), this, SLOT(SelectOutputPath()));
  QObject::connect(d->btn_submit, SIGNAL(clicked()), this, SLOT(StartRefinement()));
  QObject::connect(d->btn_interp, SIGNAL(clicked()), this, SLOT(StartInterpolate()));
  QObject::connect(d->btn_initial_bdry, SIGNAL(clicked()), this, SLOT(showImpliedBoundary()));

  //QObject::connect(d->btn_image, SIGNAL(clicked()), this, SLOT(GenerateImage()));
  //QObject::connect(d->btn_transform, SIGNAL(clicked()), this, SLOT(TransformSrep()));
}

