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

#ifndef __qSlicerSRepWarperModuleWidget_h
#define __qSlicerSRepWarperModuleWidget_h

// Slicer includes
#include "qSlicerAbstractModuleWidget.h"

#include "qSlicerSRepWarperModuleExport.h"

class qSlicerSRepWarperModuleWidgetPrivate;
class vtkMRMLNode;

/// \ingroup Slicer_QtModules_ExtensionTemplate
class Q_SLICER_QTMODULES_SREPWARPER_EXPORT qSlicerSRepWarperModuleWidget :
  public qSlicerAbstractModuleWidget
{
  Q_OBJECT

public:

  typedef qSlicerAbstractModuleWidget Superclass;
  qSlicerSRepWarperModuleWidget(QWidget *parent=0);
  virtual ~qSlicerSRepWarperModuleWidget();

public slots:
  void run();
  void setMRMLScene(vtkMRMLScene* scene) override;

protected:
  QScopedPointer<qSlicerSRepWarperModuleWidgetPrivate> d_ptr;

  virtual void setup();

private:
  Q_DECLARE_PRIVATE(qSlicerSRepWarperModuleWidget);
  Q_DISABLE_COPY(qSlicerSRepWarperModuleWidget);
};

#endif
