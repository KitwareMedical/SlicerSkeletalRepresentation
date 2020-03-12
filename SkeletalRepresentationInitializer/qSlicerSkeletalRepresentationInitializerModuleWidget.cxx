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
#include "qSlicerSkeletalRepresentationInitializerModuleWidget.h"
#include "ui_qSlicerSkeletalRepresentationInitializerModuleWidget.h"

// module logic file
#include "vtkSlicerSkeletalRepresentationInitializerLogic.h"

#include <QFileDialog>
#include <QMessageBox>

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerSkeletalRepresentationInitializerModuleWidgetPrivate: public Ui_qSlicerSkeletalRepresentationInitializerModuleWidget
{
    Q_DECLARE_PUBLIC(qSlicerSkeletalRepresentationInitializerModuleWidget)
protected:
    qSlicerSkeletalRepresentationInitializerModuleWidget* const q_ptr;
public:
  qSlicerSkeletalRepresentationInitializerModuleWidgetPrivate(qSlicerSkeletalRepresentationInitializerModuleWidget &);
    ~qSlicerSkeletalRepresentationInitializerModuleWidgetPrivate();
    vtkSlicerSkeletalRepresentationInitializerLogic* logic() const;
};

//-----------------------------------------------------------------------------
// qSlicerSkeletalRepresentationInitializerModuleWidgetPrivate methods

//-----------------------------------------------------------------------------
qSlicerSkeletalRepresentationInitializerModuleWidgetPrivate::qSlicerSkeletalRepresentationInitializerModuleWidgetPrivate(qSlicerSkeletalRepresentationInitializerModuleWidget& object) : q_ptr(&object)
{
}

qSlicerSkeletalRepresentationInitializerModuleWidgetPrivate::~qSlicerSkeletalRepresentationInitializerModuleWidgetPrivate(){}

vtkSlicerSkeletalRepresentationInitializerLogic* qSlicerSkeletalRepresentationInitializerModuleWidgetPrivate::logic() const
{
    Q_Q(const qSlicerSkeletalRepresentationInitializerModuleWidget);
    return vtkSlicerSkeletalRepresentationInitializerLogic::SafeDownCast(q->logic());
}


//-----------------------------------------------------------------------------
// qSlicerSkeletalRepresentationInitializerModuleWidget methods

//-----------------------------------------------------------------------------
qSlicerSkeletalRepresentationInitializerModuleWidget::qSlicerSkeletalRepresentationInitializerModuleWidget(QWidget* _parent)
  : Superclass( _parent )
  , d_ptr( new qSlicerSkeletalRepresentationInitializerModuleWidgetPrivate(*this) )
{
}

//-----------------------------------------------------------------------------
qSlicerSkeletalRepresentationInitializerModuleWidget::~qSlicerSkeletalRepresentationInitializerModuleWidget()
{
}

//-----------------------------------------------------------------------------
void qSlicerSkeletalRepresentationInitializerModuleWidget::setup()
{
  Q_D(qSlicerSkeletalRepresentationInitializerModuleWidget);
  d->setupUi(this);
  this->Superclass::setup();
  QObject::connect(d->SelectInputButton, SIGNAL(clicked()), this, SLOT(selectInput()));
  QObject::connect(d->btn_flow, SIGNAL(clicked()), this, SLOT(flow()));
  //QObject::connect(d->btn_one_step_flow, SIGNAL(clicked()), this, SLOT(flowOneStep()));
  //QObject::connect(d->btn_match_ell, SIGNAL(clicked()), this, SLOT(pullUpFittingEllipsoid()));
  //QObject::connect(d->btn_inkling_flow, SIGNAL(clicked()), this, SLOT(inklingFlow()));
  QObject::connect(d->btn_back_flow, SIGNAL(clicked()), this, SLOT(backwardFlow()));
  QObject::connect(d->btn_output, SIGNAL(clicked()), this, SLOT(outputPath()));
  QObject::connect(d->btn_reorder_skeleton, SIGNAL(clicked()), this, SLOT(rotateSkeleton()));
}

void qSlicerSkeletalRepresentationInitializerModuleWidget::pullUpFittingEllipsoid()
{
//    Q_D(qSlicerSkeletalRepresentationInitializerModuleWidget);
    //d->logic()->DummyShowFittingEllipsoid();

}

void qSlicerSkeletalRepresentationInitializerModuleWidget::selectInput()
{
    Q_D(qSlicerSkeletalRepresentationInitializerModuleWidget);
    QString fileName = QFileDialog::getOpenFileName(this, "Select input mesh");
    d->lb_input_file_path->setText(fileName.toUtf8().constData());
    d->logic()->SetInputFileName(fileName.toUtf8().constData());
}

void qSlicerSkeletalRepresentationInitializerModuleWidget::outputPath()
{
    Q_D(qSlicerSkeletalRepresentationInitializerModuleWidget);
    QString fileName = QFileDialog::getExistingDirectory(this, "Select refinement output folder");
    d->lb_output_path->setText(fileName.toUtf8().constData());
    d->logic()->SetOutputPath(fileName.toUtf8().constData());
}

void qSlicerSkeletalRepresentationInitializerModuleWidget::rotateSkeleton()
{
    Q_D(qSlicerSkeletalRepresentationInitializerModuleWidget);
    int numCols = static_cast<int>(d->sl_num_cols->value());
    int numRows = static_cast<int>(d->sl_num_rows->value());
    d->logic()->SetRows(numRows);
    d->logic()->SetCols(numCols);
    d->logic()->RotateSkeleton(d->cb_flip_green->isChecked(),d->cb_flip_red->isChecked(), d->cb_flip_blue->isChecked());
}

void qSlicerSkeletalRepresentationInitializerModuleWidget::flow()
{
    Q_D(qSlicerSkeletalRepresentationInitializerModuleWidget);
    std::string fileName= d->lb_input_file_path->text().toUtf8().constData();

    double dt = d->sl_dt->value();
    double smoothAmount = d->sl_smooth_amount->value();
    int maxIter = int(d->sl_max_iter->value());
    int freq_output = int(d->sl_freq_output->value());
    int numCols = static_cast<int>(d->sl_num_cols->value());
    int numRows = static_cast<int>(d->sl_num_rows->value());

    // odd rows is required
    if(numRows % 2 == 0)
    {
        numRows += 1;
    }
    // odd cols is required
    if(numCols % 2 == 0)
    {
        numCols += 1;
    }
    d->sl_num_cols->setValue(numCols);
    d->sl_num_rows->setValue(numRows);
    d->logic()->SetRows(numRows);
    d->logic()->SetCols(numCols);

    d->logic()->FlowSurfaceMesh(fileName, dt, smoothAmount, maxIter, freq_output);
}

void qSlicerSkeletalRepresentationInitializerModuleWidget::flowOneStep()
{
    Q_D(qSlicerSkeletalRepresentationInitializerModuleWidget);
    std::string fileName= d->lb_input_file_path->text().toUtf8().constData();
    double dt = d->sl_dt->value();
    double smoothAmount = d->sl_smooth_amount->value();

    d->logic()->FlowSurfaceOneStep(fileName, dt, smoothAmount);
}
void qSlicerSkeletalRepresentationInitializerModuleWidget::inklingFlow()
{
    Q_D(qSlicerSkeletalRepresentationInitializerModuleWidget);
    std::string fileName= d->lb_input_file_path->text().toUtf8().constData();

    double dt = d->sl_dt->value();
    double smoothAmount = d->sl_smooth_amount->value();
    int maxIter = int(d->sl_max_iter->value());
    int freq_output = int(d->sl_freq_output->value());
//    double threshold = d->sl_threshold->value();
    double threshold = 13.0;// for test
    d->logic()->InklingFlow(fileName, dt, smoothAmount, maxIter, freq_output, threshold);
}

void qSlicerSkeletalRepresentationInitializerModuleWidget::backwardFlow()
{
    Q_D(qSlicerSkeletalRepresentationInitializerModuleWidget);
    int maxIter = static_cast<int>(d->sl_max_iter->value());
    d->logic()->BackwardFlow(maxIter);
}
