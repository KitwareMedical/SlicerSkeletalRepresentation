#ifndef __vtkMRMLSRepNode_h
#define __vtkMRMLSRepNode_h

// MRML includes
#include "vtkMRMLDisplayableNode.h"
#include "vtkMRMLSRepDisplayNode.h"

#include "vtkSlicerSRepModuleMRMLExport.h"

#include <vtkPoints.h>

#include <memory>

namespace srep {
  class SRep;
}

class VTK_SLICER_SREP_MODULE_MRML_EXPORT vtkMRMLSRepNode : public vtkMRMLDisplayableNode
{
public:
  vtkMRMLSRepNode();
  virtual ~vtkMRMLSRepNode();

  static vtkMRMLSRepNode *New();
  vtkTypeMacro(vtkMRMLSRepNode,vtkMRMLDisplayableNode);

  void PrintSelf(ostream& os, vtkIndent indent) override;

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


  //--------------------------------------------------------------------------
  // SRep specific methods
  //--------------------------------------------------------------------------
  void LoadSRepFromFile(const std::string& filename);

  bool HasSRep() const;

  /// @param append If true, will append to points. If false, will clear points before writing to them
  void PopulateUpBoundaryPoints(vtkPoints* points, bool append = false) const;

  /// @param append If true, will append to points. If false, will clear points before writing to them
  void PopulateDownBoundaryPoints(vtkPoints* points, bool append = false) const;



private:
    std::unique_ptr<srep::SRep> SRep;
};

#endif