#ifndef __vtkMRMLEllipticalSRepNode_h
#define __vtkMRMLEllipticalSRepNode_h

#include "vtkMRMLSRepNode.h"
#include "vtkEllipticalSRep.h"

vtkEllipticalSRep* TransformSRep(vtkEllipticalSRep* srep, vtkAbstractTransform* transform);
vtkSmartPointer<vtkEllipticalSRep> SmartTransformSRep(vtkEllipticalSRep* srep, vtkAbstractTransform* transform);

class VTK_SLICER_SREP_MODULE_MRML_EXPORT vtkMRMLEllipticalSRepNode
  : public vtkMRMLSRepNode
{
public:
  vtkMRMLEllipticalSRepNode();
  vtkMRMLEllipticalSRepNode(const vtkMRMLEllipticalSRepNode&) = delete;
  vtkMRMLEllipticalSRepNode& operator=(const vtkMRMLEllipticalSRepNode&) = delete;
  vtkMRMLEllipticalSRepNode(vtkMRMLEllipticalSRepNode&&) = delete;
  vtkMRMLEllipticalSRepNode& operator=(vtkMRMLEllipticalSRepNode&&) = delete;
  virtual ~vtkMRMLEllipticalSRepNode();

  //--------------------------------------------------------------------------
  // MRMLNode methods
  //--------------------------------------------------------------------------
  vtkMRMLNode* CreateNodeInstance() override;
  /// Get node XML tag name (like Volume, Model)
  const char* GetNodeTagName() override {return "EllipticalSRep";};

  static vtkMRMLEllipticalSRepNode *New();
  vtkTypeMacro(vtkMRMLEllipticalSRepNode,vtkMRMLSRepNode);

  /// Apply the passed transformation to the SRep
  /// \sa CanApplyNonLinearTransforms
  void ApplyTransform(vtkAbstractTransform* transform) override;

  /// Copy node content (excludes basic data, such as name and node references).
  /// \sa vtkMRMLNode::CopyContent
  vtkMRMLCopyContentMacro(vtkMRMLSRepNode);

  //--------------------------------------------------------------------------
  // SRep methods
  //--------------------------------------------------------------------------
  /// @{
  /// Gets the SRep, if any, before any transforms are applied.
  /// \sa GetSRepWorld, HasSRep
  const vtkEllipticalSRep* GetEllipticalSRep() const;
  vtkEllipticalSRep* GetEllipticalSRep();
  /// @}

  /// Gets the SRep, if any, after all transforms are applied.
  /// If there are no transforms, this will be the same as GetEllipticalSRep
  /// \sa GetSRep, HasSRep
  const vtkEllipticalSRep* GetEllipticalSRepWorld() const;

  /// Sets the SRep. Takes sole ownership.
  void SetEllipticalSRep(vtkEllipticalSRep* srep);

  /// Gets the SRep before any non-hardened transforms are applied.
  const vtkMeshSRepInterface* GetSRep() const override;

  /// Gets the SRep after all non-hardened transforms are applied.
  const vtkMeshSRepInterface* GetSRepWorld() const override;

protected:
  void DoUpdateSRepWorld(vtkAbstractTransform* transform) override;

private:
  // using shared_ptr to allow easy shallow copy in CopyContent
  vtkSmartPointer<vtkEllipticalSRep> SRep;
  unsigned long SRepObservationTag;
  vtkSmartPointer<vtkEllipticalSRep> SRepWorld;

  void onSRepModified(vtkObject *caller, unsigned long event, void* callData);
};

#endif
