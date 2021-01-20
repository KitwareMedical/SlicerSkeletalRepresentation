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

#ifndef _qSlicerSRepInitializerModule_h
#define _qSlicerSRepInitializerModule_h

// Slicer includes
#include "qSlicerLoadableModule.h"

// SRepInitializer includes
#include "qSlicerSRepInitializerModuleExport.h"

class qSlicerSRepInitializerModulePrivate;

/// \ingroup Slicer_QtModules_ExtensionTemplate
class Q_SLICER_QTMODULES_SREPINITIALIZER_EXPORT
qSlicerSRepInitializerModule
  : public qSlicerLoadableModule
{
  Q_OBJECT
#ifdef Slicer_HAVE_QT5
  Q_PLUGIN_METADATA(IID "org.slicer.modules.loadable.qSlicerLoadableModule/1.0");
#endif
  Q_INTERFACES(qSlicerLoadableModule)

public:

  typedef qSlicerLoadableModule Superclass;
  explicit qSlicerSRepInitializerModule(QObject *parent=nullptr);
  virtual ~qSlicerSRepInitializerModule() override;

  qSlicerGetTitleMacro(QTMODULE_TITLE)

  virtual QString helpText()const override;
  virtual QString acknowledgementText()const override;
  virtual QStringList contributors()const override;

  virtual QIcon icon()const override;

  virtual QStringList categories()const override;
  virtual QStringList dependencies() const override;

protected:

  /// Initialize the module. Register the volumes reader/writer
  virtual void setup() override;

  /// Create and return the widget representation associated to this module
  virtual qSlicerAbstractModuleRepresentation * createWidgetRepresentation() override;

  /// Create and return the logic associated to this module
  virtual vtkMRMLAbstractLogic* createLogic() override;

protected:
  QScopedPointer<qSlicerSRepInitializerModulePrivate> d_ptr;

private:
  Q_DECLARE_PRIVATE(qSlicerSRepInitializerModule)
  Q_DISABLE_COPY(qSlicerSRepInitializerModule)

};

#endif
