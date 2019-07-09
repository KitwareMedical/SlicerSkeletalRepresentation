/*==============================================================================

  Program: 3D Slicer

  Portions (c) Copyright Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

==============================================================================*/

// .NAME vtkSlicerSkeletalRepresentationInitializerLogic - slicer logic class for volumes manipulation
// .SECTION Description
// This class manages the logic associated with reading, saving,
// and changing propertied of the volumes


#ifndef __vtkSlicerSkeletalRepresentationInitializerLogic_h
#define __vtkSlicerSkeletalRepresentationInitializerLogic_h

// Slicer includes
#include <vtkSlicerModuleLogic.h>

// SkeletalRepresentationInitializer Logic includes
#include "vtkSlicerSkeletalRepresentationInitializerModuleLogicExport.h"
#include <itkThinPlateSplineExtended.h>

// ITK includes
#include <itkPointSet.h>

// VTK includes
class vtkPolyData;
class vtkPoints;
class vtkCellArray;

// STD includes
#include <cstdlib>

class VTK_SLICER_SKELETALREPRESENTATIONINITIALIZER_MODULE_LOGIC_EXPORT vtkSlicerSkeletalRepresentationInitializerLogic :
  public vtkSlicerModuleLogic
{
public:
  typedef double CoordinateRepType;
  typedef itk::Point< CoordinateRepType, 3 > PointType;
  static vtkSlicerSkeletalRepresentationInitializerLogic *New();
  vtkTypeMacro(vtkSlicerSkeletalRepresentationInitializerLogic, vtkSlicerModuleLogic);
  void PrintSelf(ostream& os, vtkIndent indent);

  // flow surface to the end: either it's ellipsoidal enough or reach max_itre
  // input[dt]: delta t in each move
  // input[smooth_amount]: 0-2 double value for smooth filter
  // input[max_iter]: maximum of iteration number
  // input[freq_output]: the frequence of output model(intermediate surface mesh of flow) node to scene
  int FlowSurfaceMesh(const std::string &filename, double dt, double smooth_amount, int max_iter, int freq_output);

  // flow one step only
  // input[dt]: delta t in each move
  // input[smooth_amount]: 0-2 double value for smooth filter
  int FlowSurfaceOneStep(const std::string &filename, double dt, double smooth_amount);

  // Select input mesh and render it in scene
  // input[filename]: whole path of vtk file
  void SetInputFileName(const std::string &filename);

  // Show fitting ellipsoid in 3D window
  // Can be called after one step flow or overall flow
  // if called after one step flow,  render the ellipsoid generated just now
  // otherwise render the ellipsoid at the end
  // input the forward flow deformed mesh, output radii ( double &rx, double &ry, double &rz)
  void ShowFittingEllipsoid(vtkPolyData* mesh, double &rx, double &ry, double &rz);

  // generate srep given an ellipsoid and expected rows and columns of medial sheet.
  void GenerateSrepForEllipsoid(vtkPolyData* mesh, int rows, int cols, int forwardCount);

  int InklingFlow(const std::string &filename, double dt, double smooth_amount, int max_iter, int freq_output, double threshold);

  // Real backward flow
  // input: totalNum of surface files from forward flow
  //
  // output: files and srep of initial object
  void BackwardFlow(int totalNum);
  
protected:
  vtkSlicerSkeletalRepresentationInitializerLogic();
  virtual ~vtkSlicerSkeletalRepresentationInitializerLogic();

  virtual void SetMRMLSceneInternal(vtkMRMLScene* newScene);
  /// Register MRML Node classes to Scene. Gets called automatically when the MRMLScene is attached to this logic class.
  virtual void RegisterNodes();
  virtual void UpdateFromMRMLScene();
  virtual void OnMRMLSceneNodeAdded(vtkMRMLNode* node);
  virtual void OnMRMLSceneNodeRemoved(vtkMRMLNode* node);

private:
  void AddModelNodeToScene(vtkPolyData* mesh, const char* modelName, bool isModelVisible, double r = 0.25, double g = 0.25, double b = 0.25);
  void HideNodesByNameByClass(const std::string & nodeName, const std::string &className);
  void HideNodesByClass(const std::string &className);
  void AddPointToScene(double x, double y, double z, int glyphType, double r = 1, double g = 0, double b = 0);

  void ComputePairwiseTps(int totalNum);
  int ApplyTps(int totalNum);
  void DisplayResultSrep();
  void TransformNOutput(itkThinPlateSplineExtended::Pointer tps,
                       vtkPolyData* spokes, const std::string& outputFileName);
  void TransformPoints(itkThinPlateSplineExtended::Pointer tps,
                       vtkPolyData* poly, const std::string& outputFileName);
  double CalculateSpokeLength(PointType tail, PointType tip);
  void CalculateSpokeDirection(PointType tail, PointType tip, double *x, double *y, double *z);
  void GetNeighborCells(vtkPolyData* mesh, int ptId, int newId, vtkCellArray* output, vtkPoints* morePts);
private:

  vtkSlicerSkeletalRepresentationInitializerLogic(const vtkSlicerSkeletalRepresentationInitializerLogic&); // Not implemented
  void operator=(const vtkSlicerSkeletalRepresentationInitializerLogic&); // Not implemented

private:
  int forwardCount = 0;
};

#endif
