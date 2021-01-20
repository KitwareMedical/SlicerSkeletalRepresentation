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
#include "qSlicerSRepInitializerFooBarWidget.h"
#include "ui_qSlicerSRepInitializerFooBarWidget.h"

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_SRepInitializer
class qSlicerSRepInitializerFooBarWidgetPrivate
  : public Ui_qSlicerSRepInitializerFooBarWidget
{
  Q_DECLARE_PUBLIC(qSlicerSRepInitializerFooBarWidget);
protected:
  qSlicerSRepInitializerFooBarWidget* const q_ptr;

public:
  qSlicerSRepInitializerFooBarWidgetPrivate(
    qSlicerSRepInitializerFooBarWidget& object);
  virtual void setupUi(qSlicerSRepInitializerFooBarWidget*);
};

// --------------------------------------------------------------------------
qSlicerSRepInitializerFooBarWidgetPrivate
::qSlicerSRepInitializerFooBarWidgetPrivate(
  qSlicerSRepInitializerFooBarWidget& object)
  : q_ptr(&object)
{
}

// --------------------------------------------------------------------------
void qSlicerSRepInitializerFooBarWidgetPrivate
::setupUi(qSlicerSRepInitializerFooBarWidget* widget)
{
  this->Ui_qSlicerSRepInitializerFooBarWidget::setupUi(widget);
}

//-----------------------------------------------------------------------------
// qSlicerSRepInitializerFooBarWidget methods

//-----------------------------------------------------------------------------
qSlicerSRepInitializerFooBarWidget
::qSlicerSRepInitializerFooBarWidget(QWidget* parentWidget)
  : Superclass( parentWidget )
  , d_ptr( new qSlicerSRepInitializerFooBarWidgetPrivate(*this) )
{
  Q_D(qSlicerSRepInitializerFooBarWidget);
  d->setupUi(this);
}

//-----------------------------------------------------------------------------
qSlicerSRepInitializerFooBarWidget
::~qSlicerSRepInitializerFooBarWidget()
{
}
