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
#include "qSlicerSRepInitializerModuleWidget.h"
#include "ui_qSlicerSRepInitializerModuleWidget.h"

// module logic file
#include "vtkSlicerSRepInitializerLogic.h"

#include <QFileDialog>
#include <QMessageBox>

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerSRepInitializerModuleWidgetPrivate: public Ui_qSlicerSRepInitializerModuleWidget
{
    Q_DECLARE_PUBLIC(qSlicerSRepInitializerModuleWidget)
protected:
    qSlicerSRepInitializerModuleWidget* const q_ptr;
public:
  qSlicerSRepInitializerModuleWidgetPrivate(qSlicerSRepInitializerModuleWidget &);
    ~qSlicerSRepInitializerModuleWidgetPrivate();
    vtkSlicerSRepInitializerLogic* logic() const;
};

//-----------------------------------------------------------------------------
// qSlicerSRepInitializerModuleWidgetPrivate methods

//-----------------------------------------------------------------------------
qSlicerSRepInitializerModuleWidgetPrivate::qSlicerSRepInitializerModuleWidgetPrivate(qSlicerSRepInitializerModuleWidget& object) : q_ptr(&object)
{
}

qSlicerSRepInitializerModuleWidgetPrivate::~qSlicerSRepInitializerModuleWidgetPrivate(){}

vtkSlicerSRepInitializerLogic* qSlicerSRepInitializerModuleWidgetPrivate::logic() const
{
    Q_Q(const qSlicerSRepInitializerModuleWidget);
    return vtkSlicerSRepInitializerLogic::SafeDownCast(q->logic());
}


//-----------------------------------------------------------------------------
// qSlicerSRepInitializerModuleWidget methods

//-----------------------------------------------------------------------------
qSlicerSRepInitializerModuleWidget::qSlicerSRepInitializerModuleWidget(QWidget* _parent)
  : Superclass( _parent )
  , d_ptr( new qSlicerSRepInitializerModuleWidgetPrivate(*this) )
{
}

//-----------------------------------------------------------------------------
qSlicerSRepInitializerModuleWidget::~qSlicerSRepInitializerModuleWidget()
{
}

//-----------------------------------------------------------------------------
void qSlicerSRepInitializerModuleWidget::setup()
{
  Q_D(qSlicerSRepInitializerModuleWidget);
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

void qSlicerSRepInitializerModuleWidget::pullUpFittingEllipsoid()
{
//    Q_D(qSlicerSRepInitializerModuleWidget);
    //d->logic()->DummyShowFittingEllipsoid();

}

void qSlicerSRepInitializerModuleWidget::selectInput()
{
    Q_D(qSlicerSRepInitializerModuleWidget);
    QString fileName = QFileDialog::getOpenFileName(this, "Select input mesh");
    d->lb_input_file_path->setText(fileName.toUtf8().constData());
    d->logic()->SetInputFileName(fileName.toUtf8().constData());
}

void qSlicerSRepInitializerModuleWidget::outputPath()
{
    Q_D(qSlicerSRepInitializerModuleWidget);
    QString fileName = QFileDialog::getExistingDirectory(this, "Select refinement output folder");
    d->lb_output_path->setText(fileName.toUtf8().constData());
    d->logic()->SetOutputPath(fileName.toUtf8().constData());
}

void qSlicerSRepInitializerModuleWidget::rotateSkeleton()
{
    Q_D(qSlicerSRepInitializerModuleWidget);
    int numBdryPoints = static_cast<int>(d->sl_num_cols->value());
    int numRadPoints = static_cast<int>(d->sl_num_rows->value());
    int numRows = numRadPoints * 2 + 1;
    int numCols = (numBdryPoints - 2 * numRows) / 2 + 2;
    d->logic()->SetRows(numRows);
    d->logic()->SetCols(numCols);
    d->logic()->RotateSkeleton(d->cb_flip_green->isChecked(),d->cb_flip_red->isChecked(), d->cb_flip_blue->isChecked());
}

void qSlicerSRepInitializerModuleWidget::flow()
{
    Q_D(qSlicerSRepInitializerModuleWidget);
    std::string fileName= d->lb_input_file_path->text().toUtf8().constData();

    double dt = d->sl_dt->value();
    double smoothAmount = d->sl_smooth_amount->value();
    int maxIter = int(d->sl_max_iter->value());
    int freq_output = int(d->sl_freq_output->value());
    int numBdryPoints = static_cast<int>(d->sl_num_cols->value());
    int numRadPoints = static_cast<int>(d->sl_num_rows->value());
    int numRows = numRadPoints * 2 + 1;
    int numCols = (numBdryPoints - 2 * numRows) / 2 + 2;
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
    d->logic()->SetRows(numRows);
    d->logic()->SetCols(numCols);

    d->logic()->FlowSurfaceMesh(fileName, dt, smoothAmount, maxIter, freq_output);
}

void qSlicerSRepInitializerModuleWidget::flowOneStep()
{
    Q_D(qSlicerSRepInitializerModuleWidget);
    std::string fileName= d->lb_input_file_path->text().toUtf8().constData();
    double dt = d->sl_dt->value();
    double smoothAmount = d->sl_smooth_amount->value();

    d->logic()->FlowSurfaceOneStep(fileName, dt, smoothAmount);
}
void qSlicerSRepInitializerModuleWidget::inklingFlow()
{
    Q_D(qSlicerSRepInitializerModuleWidget);
    std::string fileName= d->lb_input_file_path->text().toUtf8().constData();

    double dt = d->sl_dt->value();
    double smoothAmount = d->sl_smooth_amount->value();
    int maxIter = int(d->sl_max_iter->value());
    int freq_output = int(d->sl_freq_output->value());
//    double threshold = d->sl_threshold->value();
    double threshold = 13.0;// for test
    d->logic()->InklingFlow(fileName, dt, smoothAmount, maxIter, freq_output, threshold);
}

void qSlicerSRepInitializerModuleWidget::backwardFlow()
{
    Q_D(qSlicerSRepInitializerModuleWidget);
    int maxIter = static_cast<int>(d->sl_max_iter->value());
    d->logic()->BackwardFlow(maxIter);
}
