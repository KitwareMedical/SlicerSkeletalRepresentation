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
#include "qSlicerSRepCreatorFooBarWidget.h"
#include "ui_qSlicerSRepCreatorFooBarWidget.h"

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_SRepCreator
class qSlicerSRepCreatorFooBarWidgetPrivate
  : public Ui_qSlicerSRepCreatorFooBarWidget
{
  Q_DECLARE_PUBLIC(qSlicerSRepCreatorFooBarWidget);
protected:
  qSlicerSRepCreatorFooBarWidget* const q_ptr;

public:
  qSlicerSRepCreatorFooBarWidgetPrivate(
    qSlicerSRepCreatorFooBarWidget& object);
  virtual void setupUi(qSlicerSRepCreatorFooBarWidget*);
};

// --------------------------------------------------------------------------
qSlicerSRepCreatorFooBarWidgetPrivate
::qSlicerSRepCreatorFooBarWidgetPrivate(
  qSlicerSRepCreatorFooBarWidget& object)
  : q_ptr(&object)
{
}

// --------------------------------------------------------------------------
void qSlicerSRepCreatorFooBarWidgetPrivate
::setupUi(qSlicerSRepCreatorFooBarWidget* widget)
{
  this->Ui_qSlicerSRepCreatorFooBarWidget::setupUi(widget);
}

//-----------------------------------------------------------------------------
// qSlicerSRepCreatorFooBarWidget methods

//-----------------------------------------------------------------------------
qSlicerSRepCreatorFooBarWidget
::qSlicerSRepCreatorFooBarWidget(QWidget* parentWidget)
  : Superclass( parentWidget )
  , d_ptr( new qSlicerSRepCreatorFooBarWidgetPrivate(*this) )
{
  Q_D(qSlicerSRepCreatorFooBarWidget);
  d->setupUi(this);
}

//-----------------------------------------------------------------------------
qSlicerSRepCreatorFooBarWidget
::~qSlicerSRepCreatorFooBarWidget()
{
}
