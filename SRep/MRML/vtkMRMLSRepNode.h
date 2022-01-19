#ifndef __vtkMRMLSRepNode_h
#define __vtkMRMLSRepNode_h

// MRML includes
#include "vtkMRMLDisplayableNode.h"
#include "vtkMRMLSRepDisplayNode.h"

#include "vtkSlicerSRepModuleMRMLExport.h"

#include <vtkPoints.h>

#include <memory>

#include "vtkMeshSRepInterface.h"

class vtkAbstractTransform;
class vtkGeneralTransform;
class vtkMRMLTransformNode;

class VTK_SLICER_SREP_MODULE_MRML_EXPORT vtkMRMLSRepNode : public vtkMRMLDisplayableNode
{
public:
  vtkMRMLSRepNode();
  vtkMRMLSRepNode(const vtkMRMLSRepNode&) = delete;
  vtkMRMLSRepNode& operator=(const vtkMRMLSRepNode&) = delete;
  vtkMRMLSRepNode(vtkMRMLSRepNode&&) = delete;
  vtkMRMLSRepNode& operator=(vtkMRMLSRepNode&&) = delete;
  virtual ~vtkMRMLSRepNode();

  vtkTypeMacro(vtkMRMLSRepNode,vtkMRMLDisplayableNode);

  void PrintSelf(ostream& os, vtkIndent indent) override;

  /// Gets the world coordinate bounds of this SRep.
  ///
  /// This is the bounds of the SRep after all transforms have been applied.
  void GetRASBounds(double bounds[6]) override;

  /// Gets the bounds of this SRep.
  ///
  /// This is the bounds of the SRep before any non-hardened transforms have been applied.
  void GetBounds(double bounds[6]) override;

  /// @{
  /// Gets the bounds of an SRep.
  static void GetSRepBounds(const vtkMeshSRepInterface* srep, double bounds[6]);
  static void GetSRepBounds(const vtkMeshSRepInterface& srep, double bounds[6]);
  /// @}

  //--------------------------------------------------------------------------
  // MRMLNode methods
  //--------------------------------------------------------------------------
  /// Get node XML tag name (like Volume, Model)
  const char* GetNodeTagName() override {return "SRep";};

  /// Create and observe default display node(s)
  void CreateDefaultDisplayNodes() override;

  /// Creates the default storage node.
  ///
  /// \returns An owning pointer to a storage node capable of reading/writing SReps.
  vtkMRMLStorageNode* CreateDefaultStorageNode() override;

  /// Return a cast display node, returns nullptr if none
  vtkMRMLSRepDisplayNode* GetSRepDisplayNode();

  /// Returns true since can apply non linear transforms
  /// \sa ApplyTransform
  bool CanApplyNonLinearTransforms() const override;

  void OnTransformNodeReferenceChanged(vtkMRMLTransformNode* transformNode) override;

  void ProcessMRMLEvents(vtkObject * caller,
                         unsigned long event,
                         void * /*callData*/ ) override;

  //--------------------------------------------------------------------------
  // SRep specific methods
  //--------------------------------------------------------------------------

  /// Returns true if the MRML node has an SRep, false otherwise
  bool HasSRep() const;

  /// Gets the SRep, if any, before any transforms are applied.
  /// \sa GetSRepWorld, HasSRep
  virtual const vtkMeshSRepInterface* GetSRep() const = 0;

  /// Gets the SRep, if any, after all transforms are applied.
  /// If there are no transforms, this will be the same as GetSRep
  /// \sa GetSRep, HasSRep
  virtual const vtkMeshSRepInterface* GetSRepWorld() const = 0;

  /// Copy node content (excludes basic data, such as name and node references).
  /// \sa vtkMRMLNode::CopyContent
  vtkMRMLCopyContentMacro(vtkMRMLSRepNode);

protected:
  /// Call this to update the srep world with the stored transform.
  void UpdateSRepWorld();

  /// Transform the SRep.
  ///
  /// \param transform The transform to use. If nullptr, same as identity (i.e. no transform).
  virtual void DoUpdateSRepWorld(vtkAbstractTransform* transform) = 0;

private:
  vtkNew<vtkGeneralTransform> SRepTransform; // nullptr means no transform.
};

#endif
