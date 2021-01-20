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

#ifndef __qSlicerSRepRefinerModuleWidget_h
#define __qSlicerSRepRefinerModuleWidget_h

// SlicerQt includes
#include "qSlicerAbstractModuleWidget.h"

#include "qSlicerSRepRefinerModuleExport.h"

class qSlicerSRepRefinerModuleWidgetPrivate;
class vtkMRMLNode;

/// \ingroup Slicer_QtModules_ExtensionTemplate
class Q_SLICER_QTMODULES_SREPREFINER_EXPORT qSlicerSRepRefinerModuleWidget :
  public qSlicerAbstractModuleWidget
{
  Q_OBJECT

public:

  typedef qSlicerAbstractModuleWidget Superclass;
  qSlicerSRepRefinerModuleWidget(QWidget *parent=nullptr);
  virtual ~qSlicerSRepRefinerModuleWidget();

public slots:
  // select image
  void SelectImage();
  // select srep model
  void SelectSrep();
  // select output path
  void SelectOutputPath();
  // start refinement
  void StartRefinement();
  // interpolate
  void StartInterpolate();
  // generate anti-aliased image from surfacemesh
  void GenerateImage();
  // transform srep into unit cube
  void TransformSrep();

  // show initial implied boundary
  void showImpliedBoundary();

protected:
  QScopedPointer<qSlicerSRepRefinerModuleWidgetPrivate> d_ptr;

  virtual void setup();

private:
  Q_DECLARE_PRIVATE(qSlicerSRepRefinerModuleWidget);
  Q_DISABLE_COPY(qSlicerSRepRefinerModuleWidget);
};

#endif
