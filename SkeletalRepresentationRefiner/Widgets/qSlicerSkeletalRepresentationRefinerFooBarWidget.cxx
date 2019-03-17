/*==============================================================================

  Program: 3D Slicer

  Copyright (c) Kitware Inc.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

  This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc.
  and was partially funded by NIH grant 3P41RR013218-12S1

==============================================================================*/

// FooBar Widgets includes
#include "qSlicerSkeletalRepresentationRefinerFooBarWidget.h"
#include "ui_qSlicerSkeletalRepresentationRefinerFooBarWidget.h"

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_SkeletalRepresentationRefiner
class qSlicerSkeletalRepresentationRefinerFooBarWidgetPrivate
  : public Ui_qSlicerSkeletalRepresentationRefinerFooBarWidget
{
  Q_DECLARE_PUBLIC(qSlicerSkeletalRepresentationRefinerFooBarWidget);
protected:
  qSlicerSkeletalRepresentationRefinerFooBarWidget* const q_ptr;

public:
  qSlicerSkeletalRepresentationRefinerFooBarWidgetPrivate(
    qSlicerSkeletalRepresentationRefinerFooBarWidget& object);
  virtual void setupUi(qSlicerSkeletalRepresentationRefinerFooBarWidget*);
};

// --------------------------------------------------------------------------
qSlicerSkeletalRepresentationRefinerFooBarWidgetPrivate
::qSlicerSkeletalRepresentationRefinerFooBarWidgetPrivate(
  qSlicerSkeletalRepresentationRefinerFooBarWidget& object)
  : q_ptr(&object)
{
}

// --------------------------------------------------------------------------
void qSlicerSkeletalRepresentationRefinerFooBarWidgetPrivate
::setupUi(qSlicerSkeletalRepresentationRefinerFooBarWidget* widget)
{
  this->Ui_qSlicerSkeletalRepresentationRefinerFooBarWidget::setupUi(widget);
}

//-----------------------------------------------------------------------------
// qSlicerSkeletalRepresentationRefinerFooBarWidget methods

//-----------------------------------------------------------------------------
qSlicerSkeletalRepresentationRefinerFooBarWidget
::qSlicerSkeletalRepresentationRefinerFooBarWidget(QWidget* parentWidget)
  : Superclass( parentWidget )
  , d_ptr( new qSlicerSkeletalRepresentationRefinerFooBarWidgetPrivate(*this) )
{
  Q_D(qSlicerSkeletalRepresentationRefinerFooBarWidget);
  d->setupUi(this);
}

//-----------------------------------------------------------------------------
qSlicerSkeletalRepresentationRefinerFooBarWidget
::~qSlicerSkeletalRepresentationRefinerFooBarWidget()
{
}
