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
    Q_DECLARE_PUBLIC(qSlicerSkeletalRepresentationInitializerModuleWidget);
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
}

void qSlicerSkeletalRepresentationInitializerModuleWidget::selectInput()
{
    Q_D(qSlicerSkeletalRepresentationInitializerModuleWidget);
    QString filename = QFileDialog::getOpenFileName(this, "Select input mesh");
    d->logic()->FlowSurfaceMesh(filename.toUtf8().constData());
}
