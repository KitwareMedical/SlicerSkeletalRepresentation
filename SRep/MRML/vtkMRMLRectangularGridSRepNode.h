#ifndef __vtkMRMLRectangularGridSRepNode_h
#define __vtkMRMLRectangularGridSRepNode_h

#include "vtkMRMLSRepNode.h"

class VTK_SLICER_SREP_MODULE_MRML_EXPORT vtkMRMLRectangularGridSRepNode
  : public vtkMRMLSRepNode
{
public:
  vtkMRMLRectangularGridSRepNode();
  vtkMRMLRectangularGridSRepNode(const vtkMRMLRectangularGridSRepNode&) = delete;
  vtkMRMLRectangularGridSRepNode& operator=(const vtkMRMLRectangularGridSRepNode&) = delete;
  vtkMRMLRectangularGridSRepNode(vtkMRMLRectangularGridSRepNode&&) = delete;
  vtkMRMLRectangularGridSRepNode& operator=(vtkMRMLRectangularGridSRepNode&&) = delete;
  virtual ~vtkMRMLRectangularGridSRepNode();

  //--------------------------------------------------------------------------
  // MRMLNode methods
  //--------------------------------------------------------------------------
  vtkMRMLNode* CreateNodeInstance() override;
  /// Get node XML tag name (like Volume, Model)
  const char* GetNodeTagName() override {return "RectangularGridSRep";};

  static vtkMRMLRectangularGridSRepNode *New();
  vtkTypeMacro(vtkMRMLRectangularGridSRepNode,vtkMRMLSRepNode);

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
  const srep::RectangularGridSRep* GetRectangularGridSRep() const;

  /// Gets the SRep, if any, after all transforms are applied.
  /// If there are no transforms, this will be the same as GetSRep
  /// \sa GetSRep, HasSRep
  const srep::RectangularGridSRep* GetRectangularGridSRepWorld() const;

  const srep::MeshSRepInterface* GetSRep() const override;
  const srep::MeshSRepInterface* GetSRepWorld() const override;

  /// Loads SRep from a file. 
  /// \sa WriteSRepToFiles
  void LoadRectangularGridSRepFromFile(const std::string& filename);

  /// Writes SRep to a file.
  ///
  /// Will non-transformed SRep. (Any hardened transformation will still apply).
  /// \sa LoadSRepFromFile, ApplyTransform
  bool WriteRectangularGridSRepToFiles(const std::string& headerFilename,
                        const std::string& upFilename,
                        const std::string& downFilename,
                        const std::string& crestFilename);

protected:
  void DoUpdateSRepWorld(vtkAbstractTransform* transform) override;

private:
  // using shared_ptr to allow easy shallow copy in CopyContent
  std::shared_ptr<srep::RectangularGridSRep> SRep;
  std::shared_ptr<srep::RectangularGridSRep> SRepWorld;
};

#endif
