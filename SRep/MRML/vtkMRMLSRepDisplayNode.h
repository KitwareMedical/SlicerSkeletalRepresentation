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

  //--------------------------------------------------------------------------
  // display methods
  //--------------------------------------------------------------------------

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

  void SetUpSpokeVisibility(bool visible);
  bool GetUpSpokeVisibility() const;
  void SetUpSpokeColor(const vtkColor3ub& color);
  const vtkColor3ub& GetUpSpokeColor() const;

  void SetDownSpokeVisibility(bool visible);
  bool GetDownSpokeVisibility() const;
  void SetDownSpokeColor(const vtkColor3ub& color);
  const vtkColor3ub& GetDownSpokeColor() const;

  void SetCrestSpokeVisibility(bool visible);
  bool GetCrestSpokeVisibility() const;
  void SetCrestSpokeColor(const vtkColor3ub& color);
  const vtkColor3ub& GetCrestSpokeColor() const;

  void SetCrestCurveVisibility(bool visible);
  bool GetCrestCurveVisibility() const;
  void SetCrestCurveColor(const vtkColor3ub& color);
  const vtkColor3ub& GetCrestCurveColor() const;

  void SetSkeletalSheetVisibility(bool visible);
  bool GetSkeletalSheetVisibility() const;
  void SetSkeletalSheetColor(const vtkColor3ub& color);
  const vtkColor3ub& GetSkeletalSheetColor() const;

  void SetSkeletonToCrestConnectionVisibility(bool visible);
  bool GetSkeletonToCrestConnectionVisibility() const;
  void SetSkeletonToCrestConnectionColor(const vtkColor3ub& color);
  const vtkColor3ub& GetSkeletonToCrestConnectionColor() const;

  void SetSpineVisibility(bool visible);
  bool GetSpineVisibility() const;
  void SetSpineColor(const vtkColor3ub& color);
  const vtkColor3ub& GetSpineColor() const;

  void SetAbsoluteThickness(double absoluteThickness);
  double GetAbsoluteThickness() const;

  void SetRelativeThickness(double relativeThickness);
  double GetRelativeThickness() const;

  bool GetUseAbsoluteThickness() const;
  void SetUseAbsoluteThickness(bool use);
  void UseAbsoluteThicknessOn();
  void UseAbsoluteThicknessOff();

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
  DisplayHelper UpSpoke;
  DisplayHelper DownSpoke;
  DisplayHelper CrestSpoke;
  DisplayHelper CrestCurve;
  DisplayHelper SkeletalSheet;
  DisplayHelper SkeletonToCrestConnection;
  DisplayHelper Spine;
  double RelativeThickness;
  double AbsoluteThickness;
  bool UseAbsoluteThickness;
};

#endif
