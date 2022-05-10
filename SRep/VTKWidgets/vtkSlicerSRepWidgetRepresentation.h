#ifndef vtkSlicerSRepRepresentation_h
#define vtkSlicerSRepRepresentation_h

#include "vtkSlicerSRepModuleVTKWidgetsExport.h"
#include "vtkMRMLAbstractWidgetRepresentation.h"

#include "vtkMRMLSRepDisplayNode.h"
#include "vtkMRMLSRepNode.h"

#include <vtkSmartPointer.h>

class vtkActor;
class vtkCellArray;
class vtkGlyph3D;
class vtkPoints;
class vtkPolyData;
class vtkPolyDataMapper;
class vtkProperty;
class vtkSphereSource;
class vtkTubeFilter;
class vtkUnsignedCharArray;

class VTK_SLICER_SREP_MODULE_VTKWIDGETS_EXPORT vtkSlicerSRepWidgetRepresentation
  : public vtkMRMLAbstractWidgetRepresentation
{
public:
  static vtkSlicerSRepWidgetRepresentation *New();

  /// Standard methods for instances of this class.
  vtkTypeMacro(vtkSlicerSRepWidgetRepresentation, vtkMRMLAbstractWidgetRepresentation);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  void SetSRepDisplayNode(vtkMRMLSRepDisplayNode* srepDisplayNode);
  vtkMRMLSRepDisplayNode* GetSRepDisplayNode() const;
  vtkMRMLSRepNode* GetSRepNode();

  bool IsDisplayable();

  /// Update the representation from srep node
  void UpdateFromMRML(vtkMRMLNode* caller, unsigned long event, void *callData = nullptr) override;

  /// Methods to make this class behave as a vtkProp.
  void GetActors(vtkPropCollection*) override;
  void ReleaseGraphicsResources(vtkWindow*) override;
  int RenderOverlay(vtkViewport* viewport) override;
  int RenderOpaqueGeometry(vtkViewport* viewport) override;
  int RenderTranslucentPolygonalGeometry(vtkViewport* viewport) override;
  vtkTypeBool HasTranslucentPolygonalGeometry() override;

protected:
  vtkSlicerSRepWidgetRepresentation();
  ~vtkSlicerSRepWidgetRepresentation();

private:
  struct PointsRep {
      vtkSmartPointer<vtkSphereSource>      GlyphSourceSphere;
      vtkSmartPointer<vtkGlyph3D>           Glypher;
      vtkSmartPointer<vtkPolyData>          PointsPolyData;
      vtkSmartPointer<vtkProperty>          Property;
      vtkSmartPointer<vtkPolyDataMapper>    Mapper;
      vtkSmartPointer<vtkActor>             Actor;

      vtkSmartPointer<vtkTubeFilter>     TubeFilter;
      vtkSmartPointer<vtkPolyDataMapper> TubeMapper;
      vtkSmartPointer<vtkActor>          TubeActor;

      PointsRep();
      ~PointsRep();
      void SetPolyData(vtkSmartPointer<vtkPolyData> polyData);
  };

  PointsRep Skeleton;
  vtkMRMLSRepDisplayNode* SRepDisplayNode;
};

#endif
