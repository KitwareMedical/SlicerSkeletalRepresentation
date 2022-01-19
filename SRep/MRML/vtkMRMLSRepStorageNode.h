#ifndef __vtkMRMLSRepJsonStorageNode_h
#define __vtkMRMLSRepJsonStorageNode_h

#include "vtkSlicerSRepModuleMRMLExport.h"
#include "vtkMRMLStorageNode.h"
#include "vtkMRMLSRepNode.h"

class VTK_SLICER_SREP_MODULE_MRML_EXPORT vtkMRMLSRepStorageNode : public vtkMRMLStorageNode
{
public:
  static vtkMRMLSRepStorageNode *New();
  vtkTypeMacro(vtkMRMLSRepStorageNode, vtkMRMLStorageNode);

  vtkMRMLNode* CreateNodeInstance() override;

  /// Get node XML tag name (like Storage, Model)
  const char* GetNodeTagName() override {return "SRepStorage";};

  bool CanReadInReferenceNode(vtkMRMLNode *refNode) override;

  /// Creates a new SRep node with the given name and the file set by SetFileName.
  ///
  /// SRep is created in the same scene as this storage node.
  vtkMRMLSRepNode* CreateSRepNode(const char* nodeName);

  /// Gets the MRML node type of the SRep with the given file name
  ///
  /// The return value, if not empty, is suitable to be passed into
  /// vtkMRMLScene::AddNewNodeByClass as the class name for the same
  /// scene this node is in.
  /// \returns MRML node type of srep, empty string if no file name is set.
  std::string GetSRepType();

  using SRepCoordinateSystemType = int;

  /// @{
  /// Get set the coordinate system to write in.
  ///
  /// Choose either vtkMRMLStorageNode::CoordinateSystemRAS or vtkMRMLStorageNode::CoordinateSystemLPS
  /// Default is LPS.
  void SetCoordinateSystemWrite(SRepCoordinateSystemType system);
  SRepCoordinateSystemType GetCoordinateSystemWrite() const;
  void CoordinateSystemWriteRASOn();
  void CoordinateSystemWriteLPSOn();
  /// @}

protected:
  vtkMRMLSRepStorageNode();
  ~vtkMRMLSRepStorageNode() override;
  vtkMRMLSRepStorageNode(const vtkMRMLSRepStorageNode&) = delete;
  vtkMRMLSRepStorageNode(vtkMRMLSRepStorageNode&&) = delete;
  void operator=(const vtkMRMLSRepStorageNode&) = delete;
  void operator=(vtkMRMLSRepStorageNode&&) = delete;

  /// Initialize all the supported write file types
  void InitializeSupportedReadFileTypes() override;

  /// Initialize all the supported write file types
  void InitializeSupportedWriteFileTypes() override;

  /// Read data and set it in the referenced node
  int ReadDataInternal(vtkMRMLNode *refNode) override;

  /// Write data from a  referenced node.
  int WriteDataInternal(vtkMRMLNode *refNode) override;
private:
  int CoordinateSystemWrite;
};

#endif
