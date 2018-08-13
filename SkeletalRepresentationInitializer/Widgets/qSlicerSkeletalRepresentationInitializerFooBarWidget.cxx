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
#include "qSlicerSkeletalRepresentationInitializerFooBarWidget.h"
#include "ui_qSlicerSkeletalRepresentationInitializerFooBarWidget.h"

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_SkeletalRepresentationInitializer
class qSlicerSkeletalRepresentationInitializerFooBarWidgetPrivate
  : public Ui_qSlicerSkeletalRepresentationInitializerFooBarWidget
{
  Q_DECLARE_PUBLIC(qSlicerSkeletalRepresentationInitializerFooBarWidget);
protected:
  qSlicerSkeletalRepresentationInitializerFooBarWidget* const q_ptr;

public:
  qSlicerSkeletalRepresentationInitializerFooBarWidgetPrivate(
    qSlicerSkeletalRepresentationInitializerFooBarWidget& object);
  virtual void setupUi(qSlicerSkeletalRepresentationInitializerFooBarWidget*);
};

// --------------------------------------------------------------------------
qSlicerSkeletalRepresentationInitializerFooBarWidgetPrivate
::qSlicerSkeletalRepresentationInitializerFooBarWidgetPrivate(
  qSlicerSkeletalRepresentationInitializerFooBarWidget& object)
  : q_ptr(&object)
{
}

// --------------------------------------------------------------------------
void qSlicerSkeletalRepresentationInitializerFooBarWidgetPrivate
::setupUi(qSlicerSkeletalRepresentationInitializerFooBarWidget* widget)
{
  this->Ui_qSlicerSkeletalRepresentationInitializerFooBarWidget::setupUi(widget);
}

//-----------------------------------------------------------------------------
// qSlicerSkeletalRepresentationInitializerFooBarWidget methods

//-----------------------------------------------------------------------------
qSlicerSkeletalRepresentationInitializerFooBarWidget
::qSlicerSkeletalRepresentationInitializerFooBarWidget(QWidget* parentWidget)
  : Superclass( parentWidget )
  , d_ptr( new qSlicerSkeletalRepresentationInitializerFooBarWidgetPrivate(*this) )
{
  Q_D(qSlicerSkeletalRepresentationInitializerFooBarWidget);
  d->setupUi(this);
}

//-----------------------------------------------------------------------------
qSlicerSkeletalRepresentationInitializerFooBarWidget
::~qSlicerSkeletalRepresentationInitializerFooBarWidget()
{
}
