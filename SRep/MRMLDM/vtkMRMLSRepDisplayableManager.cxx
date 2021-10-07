#include "vtkMRMLSRepDisplayableManager.h"
#include <vtkEventBroker.h>
#include <vtkObjectFactory.h>
#include <vtkMRMLScene.h>
#include <vtkMRMLAbstractViewNode.h>

vtkStandardNewMacro(vtkMRMLSRepDisplayableManager);

//---------------------------------------------------------------------------
void vtkMRMLSRepDisplayableManager::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);

  os << indent << "SRepDisplayableManager" << std::endl;
}

vtkMRMLSRepDisplayableManager::vtkMRMLSRepDisplayableManager()
    : vtkMRMLAbstractDisplayableManager()
    , SRepNodes()
    , DisplayNodesToWidgets()
    , ObservedSRepNodeEvents({vtkCommand::ModifiedEvent
                            , vtkMRMLTransformableNode::TransformModifiedEvent
                            , vtkMRMLDisplayableNode::DisplayModifiedEvent
    })
{}
vtkMRMLSRepDisplayableManager::~vtkMRMLSRepDisplayableManager() = default;

void vtkMRMLSRepDisplayableManager::ProcessMRMLNodesEvents(vtkObject *caller, unsigned long event, void *callData) {
  vtkMRMLSRepNode * srepNode = vtkMRMLSRepNode::SafeDownCast(caller);
  if (srepNode) {
    bool renderRequested = false;
    for (int i = 0; i < srepNode->GetNumberOfDisplayNodes(); ++i) {
      vtkMRMLSRepDisplayNode* displayNode = vtkMRMLSRepDisplayNode::SafeDownCast(srepNode->GetNthDisplayNode(i));
      auto wit = this->DisplayNodesToWidgets.find(displayNode);
      if (this->DisplayNodesToWidgets.end() == wit) {
        this->AddDisplayNode(displayNode);
        wit = this->DisplayNodesToWidgets.find(displayNode);
        if (this->DisplayNodesToWidgets.end() == wit) {
          vtkErrorMacro("vtkMRMLSRepDisplayableManager: can't find widget right after adding it");
          continue;
        }
      }
      auto& widget = wit->second;
      widget->UpdateFromMRML(srepNode, event, callData);
      if (widget->GetNeedToRender()) {
        renderRequested = true;
        widget->NeedToRenderOff();
      }
    }

    if (renderRequested) {
      this->BatchSafeRequestRender();
    }
  } else {
    this->Superclass::ProcessMRMLNodesEvents(caller, event, callData);
  }
}

void vtkMRMLSRepDisplayableManager::BatchSafeRequestRender() {
  if (this->GetMRMLScene() && !this->GetMRMLScene()->IsBatchProcessing()) {
    this->Superclass::RequestRender();
  }
}
void vtkMRMLSRepDisplayableManager::UpdateFromMRML() {
  // this gets called from RequestRender, so make sure to jump out quickly if possible
  if (!this->GetMRMLScene()) {
    return;
  }

  // turn off update from mrml requested, as we're doing it now, and create
  // widget requests a render which checks this flag before calling update
  // from mrml again
  this->SetUpdateFromMRMLRequested(false);

  // add any new srep nodes
  std::vector<vtkMRMLNode*> srepNodesInScene;
  this->GetMRMLScene()->GetNodesByClass("vtkMRMLSRepNode", srepNodesInScene);
  for (const auto node : srepNodesInScene) {
    auto srepNode = vtkMRMLSRepNode::SafeDownCast(node);
    if (srepNode && !this->SRepNodes.count(srepNode)) {
      this->AddSRepNode(srepNode);
    }
  }

  // add any new srep display nodes
  std::vector<vtkMRMLNode*> srepDisplayNodesInScene;
  this->GetMRMLScene()->GetNodesByClass("vtkMRMLSRepDisplayNode", srepDisplayNodesInScene);
  for (const auto node : srepDisplayNodesInScene) {
    auto srepDisplayNode = vtkMRMLSRepDisplayNode::SafeDownCast(node);
    if (srepDisplayNode && !this->DisplayNodesToWidgets.count(srepDisplayNode)) {
      this->AddDisplayNode(srepDisplayNode);
    }
  }

  // remove any srep nodes that have been removed from the mrml scene
  for (auto it = this->SRepNodes.begin(); it != this->SRepNodes.end();) {
    if (srepNodesInScene.end() == std::find(srepNodesInScene.begin(), srepNodesInScene.end(), *it)) {
      it = this->RemoveSRepNode(it);
    } else {
      ++it;
    }
  }
}

void vtkMRMLSRepDisplayableManager::OnMRMLSceneNodeAdded(vtkMRMLNode* node) {
  if (!node || !this->GetMRMLScene()) {
      return;
  }

    // if the scene is still updating, jump out
  if (this->GetMRMLScene()->IsBatchProcessing()) {
    this->SetUpdateFromMRMLRequested(true);
    return;
  }

  if (node->IsA("vtkMRMLSRepNode")) {
    this->AddSRepNode(vtkMRMLSRepNode::SafeDownCast(node));
    this->BatchSafeRequestRender();
  }
}

void vtkMRMLSRepDisplayableManager::OnMRMLSceneNodeRemoved(vtkMRMLNode* node) {
  if (!node) {
    return;
  }

  auto srepNode = vtkMRMLSRepNode::SafeDownCast(node);
  if (srepNode) {
    this->RemoveSRepNode(srepNode);
    this->BatchSafeRequestRender();
    return;
  }

  auto srepDisplayNode = vtkMRMLSRepDisplayNode::SafeDownCast(node);
  if (srepDisplayNode) {
    this->RemoveDisplayNode(srepDisplayNode);
    this->BatchSafeRequestRender();
    return;
  }
}

void vtkMRMLSRepDisplayableManager::AddSRepNode(vtkMRMLSRepNode* node) {
  if (!node) {
    vtkErrorMacro("Failed to add null SRep node");
    return;
  }

  vtkMRMLAbstractViewNode* viewNode = vtkMRMLAbstractViewNode::SafeDownCast(this->GetMRMLDisplayableNode());
  if (!viewNode) {
    vtkErrorMacro("Failed to add SRep node due to null view node");
    return;
  }

  this->AddObservations(node);
  this->SRepNodes.insert(node);

  // Add Display Nodes
  const int numDisplayNodes = node->GetNumberOfDisplayNodes();
  for (int i = 0; i < numDisplayNodes; ++i) {
    vtkMRMLSRepDisplayNode *display = vtkMRMLSRepDisplayNode::SafeDownCast(node->GetNthDisplayNode(i));
    if (display
      && display->IsDisplayableInView(viewNode->GetID()))
    {
      this->AddDisplayNode(display);
    }
    else {
      vtkErrorMacro("Not adding display node because not displayable in view");
    }
  }
}

void vtkMRMLSRepDisplayableManager::RemoveSRepNode(vtkMRMLSRepNode* node) {
  if (!node) {
    return;
  }

  this->RemoveSRepNode(this->SRepNodes.find(node));
}

vtkMRMLSRepDisplayableManager::SRepNodesSet::iterator
vtkMRMLSRepDisplayableManager::RemoveSRepNode(SRepNodesSet::iterator it) {
  if (this->SRepNodes.end() != it) {
    auto node = *it;
    // Remove associated display nodes
    for (auto wit = this->DisplayNodesToWidgets.begin(); wit != this->DisplayNodesToWidgets.end();) {
      auto& displayNode = wit->first;
      if (displayNode->GetDisplayableNode() == node) {
        wit = this->RemoveDisplayNode(wit);
      } else {
        ++wit;
      }
    }

    this->RemoveObservations(node);
    return this->SRepNodes.erase(it);
  }
  return it;
}

vtkMRMLSRepDisplayableManager::DisplayNodesToWidgetsMap::iterator
vtkMRMLSRepDisplayableManager::RemoveDisplayNode(DisplayNodesToWidgetsMap::iterator wit) {
  if (wit != this->DisplayNodesToWidgets.end()) {
    auto& widget = wit->second;
    widget->SetRenderer(nullptr);
    widget->SetRepresentation(nullptr);
    return this->DisplayNodesToWidgets.erase(wit);
  }
  return wit;
}

void vtkMRMLSRepDisplayableManager::RemoveDisplayNode(vtkMRMLSRepDisplayNode* displayNode) {
  if (!displayNode) {
    return;
  }
  this->RemoveDisplayNode(this->DisplayNodesToWidgets.find(displayNode));
}

void vtkMRMLSRepDisplayableManager::AddDisplayNode(vtkMRMLSRepDisplayNode* displayNode) {
  if (!displayNode) {
    return;
  }

  //see if we already have a widget for this display node. If so, return
  if (this->DisplayNodesToWidgets.count(displayNode)) {
    return;
  }

  auto newWidget = this->CreateWidget(displayNode);
  if (!newWidget) {
    vtkErrorMacro("vtkMRMLSRepDisplayableManager: Failed to create widget");
    return;
  }

  const auto ret = this->DisplayNodesToWidgets.insert(std::make_pair(displayNode, newWidget));
  if (!ret.second) {
    vtkErrorMacro("vtkMRMLSRepDisplayableManager: Error adding widget to map");
    return;
  }

  newWidget->UpdateFromMRML(displayNode, 0);

  this->BatchSafeRequestRender();
}

vtkSmartPointer<vtkSlicerSRepWidget> vtkMRMLSRepDisplayableManager::CreateWidget(vtkMRMLSRepDisplayNode* srepDisplayNode) {
  vtkMRMLSRepNode* srepNode = srepDisplayNode->GetSRepNode();
  if (!srepNode) {
    vtkErrorMacro("vtkMRMLSRepDisplayableManager: Error cannot create widget without srep node");
    return nullptr;
  }

  vtkMRMLAbstractViewNode* viewNode = vtkMRMLAbstractViewNode::SafeDownCast(this->GetMRMLDisplayableNode());
  vtkRenderer* renderer = this->GetRenderer();

  auto widget = vtkSmartPointer<vtkSlicerSRepWidget>::New();
  widget->SetMRMLApplicationLogic(this->GetMRMLApplicationLogic());
  widget->CreateDefaultRepresentation(srepDisplayNode, viewNode, renderer);
  return widget;
}

//---------------------------------------------------------------------------
void vtkMRMLSRepDisplayableManager::AddObservations(vtkMRMLSRepNode* node)
{
  vtkCallbackCommand* callbackCommand = this->GetMRMLNodesCallbackCommand();
  vtkEventBroker* broker = vtkEventBroker::GetInstance();
  for (auto observedSRepNodeEvent : this->ObservedSRepNodeEvents) {
    if (!broker->GetObservationExist(node, observedSRepNodeEvent, this, callbackCommand)) {
      broker->AddObservation(node, observedSRepNodeEvent, this, callbackCommand);
    }
  }
}

//---------------------------------------------------------------------------
void vtkMRMLSRepDisplayableManager::RemoveObservations(vtkMRMLSRepNode* node)
{
  vtkCallbackCommand* callbackCommand = this->GetMRMLNodesCallbackCommand();
  vtkEventBroker* broker = vtkEventBroker::GetInstance();
  for (auto observedSRepNodeEvent : this->ObservedSRepNodeEvents) {
    vtkEventBroker::ObservationVector observations;
    observations = broker->GetObservations(node, observedSRepNodeEvent, this, callbackCommand);
    broker->RemoveObservations(observations);
  }
}
