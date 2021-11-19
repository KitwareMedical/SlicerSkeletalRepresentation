#ifndef __vtkMRMLEllipticalSRepNode_h
#define __vtkMRMLEllipticalSRepNode_h

#include "vtkMRMLSRepNode.h"
#include "srep/EllipticalSRep.h"

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
  /// Gets the SRep, if any, before any transforms are applied.
  /// \sa GetSRepWorld, HasSRep
  const srep::EllipticalSRep* GetEllipticalSRep() const;

  /// Gets the SRep, if any, after all transforms are applied.
  /// If there are no transforms, this will be the same as GetSRep
  /// \sa GetSRep, HasSRep
  const srep::EllipticalSRep* GetEllipticalSRepWorld() const;

  /// Sets the SRep. Takes sole ownership.
  void SetEllipticalSRep(std::unique_ptr<srep::EllipticalSRep> srep);

  const srep::MeshSRepInterface* GetSRep() const override;
  const srep::MeshSRepInterface* GetSRepWorld() const override;

protected:
  void DoUpdateSRepWorld(vtkAbstractTransform* transform) override;

private:
  // using shared_ptr to allow easy shallow copy in CopyContent
  std::shared_ptr<srep::EllipticalSRep> SRep;
  std::shared_ptr<srep::EllipticalSRep> SRepWorld;
};

#endif
