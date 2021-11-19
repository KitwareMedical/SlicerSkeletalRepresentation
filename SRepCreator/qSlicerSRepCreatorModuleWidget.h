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

#ifndef __qSlicerSRepCreatorModuleWidget_h
#define __qSlicerSRepCreatorModuleWidget_h

// Slicer includes
#include "qSlicerAbstractModuleWidget.h"

#include "qSlicerSRepCreatorModuleExport.h"

class qSlicerSRepCreatorModuleWidgetPrivate;
class vtkMRMLNode;

/// \ingroup Slicer_QtModules_ExtensionTemplate
class Q_SLICER_QTMODULES_SREPCREATOR_EXPORT qSlicerSRepCreatorModuleWidget :
  public qSlicerAbstractModuleWidget
{
  Q_OBJECT

public:
  typedef qSlicerAbstractModuleWidget Superclass;
  qSlicerSRepCreatorModuleWidget(QWidget *parent=0);
  virtual ~qSlicerSRepCreatorModuleWidget();

public slots:
  void onRunForward();
  void setMRMLScene(vtkMRMLScene* scene) override;

protected:
  QScopedPointer<qSlicerSRepCreatorModuleWidgetPrivate> d_ptr;

  virtual void setup();

private:
  Q_DECLARE_PRIVATE(qSlicerSRepCreatorModuleWidget);
  Q_DISABLE_COPY(qSlicerSRepCreatorModuleWidget);

  void onNumFoldPointsValueChanged();
};

#endif
