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

protected:
  vtkMRMLSRepDisplayableManager();
  ~vtkMRMLSRepDisplayableManager() override;

  void ProcessMRMLNodesEvents(vtkObject *caller, unsigned long event, void *callData) override;
  void OnMRMLSceneNodeAdded(vtkMRMLNode* node) override;
  void OnMRMLSceneNodeRemoved(vtkMRMLNode* node) override;

  void BatchSafeRequestRender();
  void UpdateFromMRML() override;

private:
  using SRepNodesSet = std::set<vtkSmartPointer<vtkMRMLSRepNode>>;
  using DisplayNodesToWidgetsMap =
    std::map<vtkSmartPointer<vtkMRMLSRepDisplayNode>, vtkSmartPointer<vtkSlicerSRepWidget>>;

  void AddSRepNode(vtkMRMLSRepNode* node);
  void RemoveSRepNode(vtkMRMLSRepNode* node);
  SRepNodesSet::iterator RemoveSRepNode(SRepNodesSet::iterator it);

  void AddDisplayNode(vtkMRMLSRepDisplayNode* displayNode);
  void RemoveDisplayNode(vtkMRMLSRepDisplayNode* displayNode);
  DisplayNodesToWidgetsMap::iterator RemoveDisplayNode(DisplayNodesToWidgetsMap::iterator wit);

  void AddObservations(vtkMRMLSRepNode* node);
  void RemoveObservations(vtkMRMLSRepNode* node);

  vtkSmartPointer<vtkSlicerSRepWidget> CreateWidget(vtkMRMLSRepDisplayNode* node);

  //Members
  SRepNodesSet SRepNodes;
  DisplayNodesToWidgetsMap DisplayNodesToWidgets;
  std::vector<unsigned long> ObservedSRepNodeEvents;
};

#endif