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
#include "qSlicerSRepRefinerFooBarWidget.h"
#include "ui_qSlicerSRepRefinerFooBarWidget.h"

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_SRepRefiner
class qSlicerSRepRefinerFooBarWidgetPrivate
  : public Ui_qSlicerSRepRefinerFooBarWidget
{
  Q_DECLARE_PUBLIC(qSlicerSRepRefinerFooBarWidget);
protected:
  qSlicerSRepRefinerFooBarWidget* const q_ptr;

public:
  qSlicerSRepRefinerFooBarWidgetPrivate(
    qSlicerSRepRefinerFooBarWidget& object);
  virtual void setupUi(qSlicerSRepRefinerFooBarWidget*);
};

// --------------------------------------------------------------------------
qSlicerSRepRefinerFooBarWidgetPrivate
::qSlicerSRepRefinerFooBarWidgetPrivate(
  qSlicerSRepRefinerFooBarWidget& object)
  : q_ptr(&object)
{
}

// --------------------------------------------------------------------------
void qSlicerSRepRefinerFooBarWidgetPrivate
::setupUi(qSlicerSRepRefinerFooBarWidget* widget)
{
  this->Ui_qSlicerSRepRefinerFooBarWidget::setupUi(widget);
}

//-----------------------------------------------------------------------------
// qSlicerSRepRefinerFooBarWidget methods

//-----------------------------------------------------------------------------
qSlicerSRepRefinerFooBarWidget
::qSlicerSRepRefinerFooBarWidget(QWidget* parentWidget)
  : Superclass( parentWidget )
  , d_ptr( new qSlicerSRepRefinerFooBarWidgetPrivate(*this) )
{
  Q_D(qSlicerSRepRefinerFooBarWidget);
  d->setupUi(this);
}

//-----------------------------------------------------------------------------
qSlicerSRepRefinerFooBarWidget
::~qSlicerSRepRefinerFooBarWidget()
{
}
