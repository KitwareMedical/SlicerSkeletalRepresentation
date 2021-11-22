#ifndef __vtkMRMLSRepDisplayNode_h
#define __vtkMRMLSRepDisplayNode_h

#include "vtkSlicerSRepModuleMRMLExport.h"

#include "vtkMRMLDisplayNode.h"
#include "vtkNamedColors.h"

class vtkMRMLSRepNode;

class VTK_SLICER_SREP_MODULE_MRML_EXPORT vtkMRMLSRepDisplayNode : public vtkMRMLDisplayNode {
public:
  static vtkMRMLSRepDisplayNode *New();
  vtkTypeMacro(vtkMRMLSRepDisplayNode,vtkMRMLDisplayNode);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  vtkMRMLSRepDisplayNode();
  virtual ~vtkMRMLSRepDisplayNode();

  // wrapper around GetDisplayableNode that does the cast
  vtkMRMLSRepNode* GetSRepNode();

  //--------------------------------------------------------------------------
  // MRMLNode methods
  //--------------------------------------------------------------------------

  vtkMRMLNode* CreateNodeInstance() override;

  /// Get node XML tag name (like Volume, Markups)
  const char* GetNodeTagName() override { return "SRepDisplay"; };

  /// Get name of the default interaction context (typically the mouse)
  static const std::string GetDefaultContextName() { return ""; };

  // vtkMRMLDisplayNode has a bunch of protected members we don't want to use.
  // It also offers a bunch of functions that use these members.
  // We are are overriding all the functions that use the visibility members we don't like
  // This is why we are using ints instead of bools
  /// Currently 2D visibility is not supported.
  int GetVisibility2D() override;
  // Currently 2D visibility is not supported.
  void SetVisibility2D(int visible) override;

  int GetVisibility3D() override;
  void SetVisibility3D(int visible) override;

  void SetUpSpokeColor(const vtkColor3ub& color);
  const vtkColor3ub& GetUpSpokeColor() const;

  void SetDownSpokeColor(const vtkColor3ub& color);
  const vtkColor3ub& GetDownSpokeColor() const;

  void SetCrestSpokeColor(const vtkColor3ub& color);
  const vtkColor3ub& GetCrestSpokeColor() const;

  void SetCrestCurveColor(const vtkColor3ub& color);
  const vtkColor3ub& GetCrestCurveColor() const;

  void SetSkeletalSheetColor(const vtkColor3ub& color);
  const vtkColor3ub& GetSkeletalSheetColor() const;

  void SetSkeletonToCrestConnectionColor(const vtkColor3ub& color);
  const vtkColor3ub& GetSkeletonToCrestConnectionColor() const;

private:
  struct DisplayHelper {
    bool visible;
    vtkColor3ub color;

    DisplayHelper(bool visible_, const vtkColor3ub& color_);
    DisplayHelper();
    ~DisplayHelper() = default;
    DisplayHelper(const DisplayHelper&) = default;
    DisplayHelper& operator=(const DisplayHelper&) = default;
    DisplayHelper(DisplayHelper&&) = default;
    DisplayHelper& operator=(DisplayHelper&&) = default;
  };

  bool OverallVisibility;
  DisplayHelper UpSpokes;
  DisplayHelper DownSpokes;
  DisplayHelper CrestSpokes;
  DisplayHelper CrestCurve;
  DisplayHelper SkeletalSheet;
  DisplayHelper SkeletonToCrestConnection;
};

#endif
