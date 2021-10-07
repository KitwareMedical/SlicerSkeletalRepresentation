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

#ifndef __qSlicerSRepModuleWidget_h
#define __qSlicerSRepModuleWidget_h

// Slicer includes
#include "qSlicerAbstractModuleWidget.h"

#include "qSlicerSRepModuleExport.h"

class qSlicerSRepModuleWidgetPrivate;
class vtkMRMLNode;
class vtkMRMLSRepNode;

/// \ingroup Slicer_QtModules_ExtensionTemplate
class Q_SLICER_QTMODULES_SREP_EXPORT qSlicerSRepModuleWidget :
  public qSlicerAbstractModuleWidget
{
  Q_OBJECT
public:
  typedef qSlicerAbstractModuleWidget Superclass;
  qSlicerSRepModuleWidget(QWidget *parent=0);
  virtual ~qSlicerSRepModuleWidget();

  /// Set up the GUI from mrml when entering
  void enter() override;
  /// Disconnect from scene when exiting
  void exit() override;

  void setMRMLSRepNode(vtkMRMLSRepNode* srepnode, bool forceReconnect = false);
  void updateWidgetFromMRML();

  bool setEditedNode(vtkMRMLNode* node, QString role = QString(), QString context = QString()) override;
  double nodeEditable(vtkMRMLNode* node) override;

public slots:
  // this UI related slots
  void onImport();
  void onInputFileBrowse();
  void onOpacitySliderChanged();
  void onOpacitySpinboxChanged();
  void onVisibilityChanged();

  // MRML change related slots
  void onActiveSRepItemChanged(vtkIdType);
  void onActiveSRepMRMLNodeChanged(vtkMRMLNode* node);
  void onActiveSRepNodeDisplayModifiedEvent();
  void onActiveSRepNodeModifiedEvent();
  void onMRMLSceneEndBatchProcessEvent();
  void onMRMLSceneEndCloseEvent();
  void onMRMLSceneEndImportEvent();
  void onMRMLSceneEndRestoreEvent();
  void onNodeAddedEvent(vtkObject*, vtkObject* node);

protected:
  QScopedPointer<qSlicerSRepModuleWidgetPrivate> d_ptr;

  virtual void setup();

private:
  Q_DECLARE_PRIVATE(qSlicerSRepModuleWidget);
  Q_DISABLE_COPY(qSlicerSRepModuleWidget);
};

#endif
