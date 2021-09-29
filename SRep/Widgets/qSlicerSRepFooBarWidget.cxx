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
#include "qSlicerSRepFooBarWidget.h"
#include "ui_qSlicerSRepFooBarWidget.h"

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_SRep
class qSlicerSRepFooBarWidgetPrivate
  : public Ui_qSlicerSRepFooBarWidget
{
  Q_DECLARE_PUBLIC(qSlicerSRepFooBarWidget);
protected:
  qSlicerSRepFooBarWidget* const q_ptr;

public:
  qSlicerSRepFooBarWidgetPrivate(
    qSlicerSRepFooBarWidget& object);
  virtual void setupUi(qSlicerSRepFooBarWidget*);
};

// --------------------------------------------------------------------------
qSlicerSRepFooBarWidgetPrivate
::qSlicerSRepFooBarWidgetPrivate(
  qSlicerSRepFooBarWidget& object)
  : q_ptr(&object)
{
}

// --------------------------------------------------------------------------
void qSlicerSRepFooBarWidgetPrivate
::setupUi(qSlicerSRepFooBarWidget* widget)
{
  this->Ui_qSlicerSRepFooBarWidget::setupUi(widget);
}

//-----------------------------------------------------------------------------
// qSlicerSRepFooBarWidget methods

//-----------------------------------------------------------------------------
qSlicerSRepFooBarWidget
::qSlicerSRepFooBarWidget(QWidget* parentWidget)
  : Superclass( parentWidget )
  , d_ptr( new qSlicerSRepFooBarWidgetPrivate(*this) )
{
  Q_D(qSlicerSRepFooBarWidget);
  d->setupUi(this);
}

//-----------------------------------------------------------------------------
qSlicerSRepFooBarWidget
::~qSlicerSRepFooBarWidget()
{
}
