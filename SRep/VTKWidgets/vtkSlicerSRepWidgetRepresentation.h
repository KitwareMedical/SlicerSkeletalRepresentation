#ifndef vtkSlicerSRepRepresentation_h
#define vtkSlicerSRepRepresentation_h

#include "vtkSlicerSRepModuleVTKWidgetsExport.h"
#include "vtkMRMLAbstractWidgetRepresentation.h"

#include "vtkMRMLSRepDisplayNode.h"
#include "vtkMRMLSRepNode.h"

#include <vtkSmartPointer.h>

class vtkSphereSource;
class vtkGlyph3D;
class vtkPoints;
class vtkPolyData;
class vtkProperty;
class vtkPolyDataMapper;
class vtkActor;

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
  struct SpokesRep {
      vtkSmartPointer<vtkSphereSource>   GlyphSourceSphere;
      vtkSmartPointer<vtkGlyph3D>        Glypher;
      vtkSmartPointer<vtkPoints>         BoundaryPoints;
      vtkSmartPointer<vtkPolyData>       BoundaryPointsPolyData;
      vtkSmartPointer<vtkProperty>       Property;
      vtkSmartPointer<vtkPolyDataMapper> Mapper;
      vtkSmartPointer<vtkActor>          Actor;

      SpokesRep();
      ~SpokesRep();
  };

  vtkMRMLSRepDisplayNode* SRepDisplayNode;
  SpokesRep UpSpokes;
  SpokesRep DownSpokes;
};

#endif