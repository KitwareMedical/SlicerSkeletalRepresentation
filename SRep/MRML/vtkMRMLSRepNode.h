#ifndef __vtkMRMLSRepNode_h
#define __vtkMRMLSRepNode_h

// MRML includes
#include "vtkMRMLDisplayableNode.h"
#include "vtkMRMLSRepDisplayNode.h"

#include "vtkSlicerSRepModuleMRMLExport.h"

#include <vtkPoints.h>

#include <memory>

#include <srep/RectangularGridSRep.h>

class vtkAbstractTransform;
class vtkGeneralTransform;
class vtkMRMLTransformNode;

class VTK_SLICER_SREP_MODULE_MRML_EXPORT vtkMRMLSRepNode : public vtkMRMLDisplayableNode
{
public:
  vtkMRMLSRepNode();
  virtual ~vtkMRMLSRepNode();

  static vtkMRMLSRepNode *New();
  vtkTypeMacro(vtkMRMLSRepNode,vtkMRMLDisplayableNode);

  void PrintSelf(ostream& os, vtkIndent indent) override;

  void GetRASBounds(double bounds[6]) override;
  void GetBounds(double bounds[6]) override;
  static void GetSRepBounds(const srep::RectangularGridSRep* srep, double bounds[6]);

  //--------------------------------------------------------------------------
  // MRMLNode methods
  //--------------------------------------------------------------------------
  vtkMRMLNode* CreateNodeInstance() override;
  /// Get node XML tag name (like Volume, Model)
  const char* GetNodeTagName() override {return "SRep";};

  /// Create and observe default display node(s)
  void CreateDefaultDisplayNodes() override;

  /// Return a cast display node, returns nullptr if none
  vtkMRMLSRepDisplayNode* GetSRepDisplayNode();

  /// Returns true since can apply non linear transforms
  /// \sa ApplyTransform
  bool CanApplyNonLinearTransforms() const override;

  /// Apply the passed transformation to the SRep
  /// \sa CanApplyNonLinearTransforms
  void ApplyTransform(vtkAbstractTransform* transform) override;

  void OnTransformNodeReferenceChanged(vtkMRMLTransformNode* transformNode) override;

  void ProcessMRMLEvents(vtkObject * caller,
                         unsigned long event,
                         void * /*callData*/ ) override;

  //--------------------------------------------------------------------------
  // SRep specific methods
  //--------------------------------------------------------------------------
  /// Loads SRep from a file. 
  /// \sa WriteSRepToFiles
  void LoadSRepFromFile(const std::string& filename);

  /// Writes SRep to a file.
  ///
  /// Will non-transformed SRep. (Any hardened transformation will still apply).
  /// \sa LoadSRepFromFile, ApplyTransform
  bool WriteSRepToFiles(const std::string& headerFilename,
                        const std::string& upFilename,
                        const std::string& downFilename,
                        const std::string& crestFilename);

  /// Returns true if the MRML node has an SRep, false otherwise
  bool HasSRep() const;

  /// Gets the SRep, if any, before any transforms are applied.
  /// \sa GetSRepWorld, HasSRep
  const srep::RectangularGridSRep* GetSRep() const;

  /// Gets the SRep, if any, after all transforms are applied.
  /// If there are no transforms, this will be the same as GetSRep
  /// \sa GetSRep, HasSRep
  const srep::RectangularGridSRep* GetSRepWorld() const;

  // TODO: void SetSRep(const srep::RectangularGridSRep&);
  // TODO: void SetSRep(srep::RectangularGridSRep&&);

  /// Copy node content (excludes basic data, such as name and node references).
  /// \sa vtkMRMLNode::CopyContent
  vtkMRMLCopyContentMacro(vtkMRMLSRepNode);

private:
  void UpdateSRepWorld(vtkAbstractTransform* transform);

  // using shared_ptr to allow easy shallow copy in CopyContent
  std::shared_ptr<srep::RectangularGridSRep> SRep;
  std::shared_ptr<srep::RectangularGridSRep> SRepWorld;
  vtkNew<vtkGeneralTransform> SRepTransform; // nullptr means no transform.
};

#endif
