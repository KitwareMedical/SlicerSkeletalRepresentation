#ifndef __vtkMRMLSRepDisplayableManager_h
#define __vtkMRMLSRepDisplayableManager_h

// SRepModule includes
#include "vtkSlicerSRepModuleMRMLDisplayableManagerExport.h"
#include <vtkMRMLSRepNode.h>
#include <vtkMRMLSRepDisplayNode.h>
#include <vtkSlicerSRepWidget.h>

// MRMLDisplayableManager includes
#include <vtkMRMLAbstractDisplayableManager.h>

class VTK_SLICER_SREP_MODULE_MRMLDISPLAYABLEMANAGER_EXPORT vtkMRMLSRepDisplayableManager
    : public vtkMRMLAbstractDisplayableManager
{
public:
  static vtkMRMLSRepDisplayableManager *New();
  vtkTypeMacro(vtkMRMLSRepDisplayableManager, vtkMRMLAbstractDisplayableManager);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  void OnMRMLSceneNodeAdded(vtkMRMLNode* node) override;

protected:
  vtkMRMLSRepDisplayableManager();
  ~vtkMRMLSRepDisplayableManager() override;

  void ProcessMRMLNodesEvents(vtkObject *caller, unsigned long event, void *callData) override;

private:
  void AddSRepNode(vtkMRMLSRepNode* node);
  void AddDisplayNode(vtkMRMLSRepDisplayNode* displayNode);

  void AddObservations(vtkMRMLSRepNode* node);
  void RemoveObservations(vtkMRMLSRepNode* node);

  /// Create a widget.
  vtkSmartPointer<vtkSlicerSRepWidget> CreateWidget(vtkMRMLSRepDisplayNode* node);

  using SRepNodesSet = std::set<vtkSmartPointer<vtkMRMLSRepNode>>;
  SRepNodesSet SRepNodes;

  using DisplayNodesToWidgetsMap =
    std::map<vtkSmartPointer<vtkMRMLSRepDisplayNode>, vtkSmartPointer<vtkSlicerSRepWidget>>;
  DisplayNodesToWidgetsMap DisplayNodesToWidgets;

  std::vector<unsigned long> ObservedSRepNodeEvents;
};

#endif