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

#ifndef _qSlicerSRepInitializerModuleWidget_h
#define _qSlicerSRepInitializerModuleWidget_h

// Slicer includes
#include <qSlicerAbstractModuleWidget.h>

// SRepInitializer includes
#include "qSlicerSRepInitializerModuleExport.h"

class qSlicerSRepInitializerModuleWidgetPrivate;
class vtkMRMLNode;

/// \ingroup Slicer_QtModules_ExtensionTemplate
class Q_SLICER_QTMODULES_SREPINITIALIZER_EXPORT qSlicerSRepInitializerModuleWidget :
  public qSlicerAbstractModuleWidget
{
  Q_OBJECT

public:

  typedef qSlicerAbstractModuleWidget Superclass;
  qSlicerSRepInitializerModuleWidget(QWidget *parent=nullptr);
  virtual ~qSlicerSRepInitializerModuleWidget();

public slots:
  // for select input mesh
  void selectInput();
  // connect the button Flow to the end
  void flow();
  // connect the button Flow step by step
  void flowOneStep();
  // connect the button Match ellipsoid
  void pullUpFittingEllipsoid();
  // connect the button flow with laplacian curvature
  void inklingFlow();

  // connect the button backward flow
  void backwardFlow();

  // set output path for s-rep
  void outputPath();

  // reorder / rotate the skeleton
  void rotateSkeleton();

protected:
  QScopedPointer<qSlicerSRepInitializerModuleWidgetPrivate> d_ptr;

  virtual void setup();

private:
  Q_DECLARE_PRIVATE(qSlicerSRepInitializerModuleWidget)
  Q_DISABLE_COPY(qSlicerSRepInitializerModuleWidget)
};

#endif
