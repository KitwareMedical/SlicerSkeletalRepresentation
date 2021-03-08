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

// .NAME vtkSlicerSRepInitializerLogic - slicer logic class for volumes manipulation
// .SECTION Description
// This class manages the logic associated with reading, saving,
// and changing propertied of the volumes


#ifndef _vtkSlicerSRepInitializerLogic_h
#define _vtkSlicerSRepInitializerLogic_h

// Slicer includes
#include <vtkSlicerModuleLogic.h>

// SRepInitializer Logic includes
#include "vtkSlicerSRepInitializerModuleLogicExport.h"
#include <itkThinPlateSplineExtended.h>

// ITK includes
#include <itkPointSet.h>

// VTK includes
class vtkPolyData;
class vtkPoints;
class vtkCellArray;

// STD includes
#include <cstdlib>

class VTK_SLICER_SREPINITIALIZER_MODULE_LOGIC_EXPORT vtkSlicerSRepInitializerLogic :
  public vtkSlicerModuleLogic
{
public:
  typedef double CoordinateRepType;
  typedef itk::Point< CoordinateRepType, 3 > PointType;
  static vtkSlicerSRepInitializerLogic *New();
  vtkTypeMacro(vtkSlicerSRepInitializerLogic, vtkSlicerModuleLogic)
  void PrintSelf(ostream& os, vtkIndent indent) override;

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
  void GenerateSrepForEllipsoid(vtkPolyData* mesh, int rows, int cols, int forwardCount,
                                bool rotateX = false, bool rotateY = false,
                                bool rotateZ = false);

  int InklingFlow(const std::string &filename, double dt, double smooth_amount, int max_iter, int freq_output, double threshold);

  // Real backward flow
  // input: totalNum of surface files from forward flow
  //
  // output: files and srep of initial object
  void BackwardFlow(int totalNum);

  // set output path for initialized s-rep
  void SetOutputPath(const std::string &outputPath);

  // set rows and cols (resolution) of s-rep
  void SetRows(int r);

  void SetCols(int c);

  // display and save current initial srep
  // input: flip is true if users want to flip the orientation of the srep
  void DisplayResultSrep(bool flip = false);

  // reorder skeletal points
  void RotateSkeleton(bool rotateX, bool rotateY, bool rotateZ);
  void ReorderSpokes(vtkPolyData* input, vtkPoints* outputPts, vtkCellArray* outputPolys);
protected:
  vtkSlicerSRepInitializerLogic();
  virtual ~vtkSlicerSRepInitializerLogic() override;

  virtual void SetMRMLSceneInternal(vtkMRMLScene* newScene) override;
  /// Register MRML Node classes to Scene. Gets called automatically when the MRMLScene is attached to this logic class.
  virtual void RegisterNodes() override;
  virtual void UpdateFromMRMLScene() override;
  virtual void OnMRMLSceneNodeAdded(vtkMRMLNode* node) override;
  virtual void OnMRMLSceneNodeRemoved(vtkMRMLNode* node) override;

private:
  void AddModelNodeToScene(vtkPolyData* mesh, const char* modelName, bool isModelVisible, double r = 0.25, double g = 0.25, double b = 0.25);
  void HideNodesByNameByClass(const std::string & nodeName, const std::string &className);
  void HideNodesByClass(const std::string &className);
  void AddPointToScene(double x, double y, double z, int glyphType, double r = 1, double g = 0, double b = 0);

  void ComputePairwiseTps(int totalNum);

  void TransformNOutput(itkThinPlateSplineExtended::Pointer tps,
                       vtkPolyData* spokes, const std::string& outputFileName);
  void TransformPoints(itkThinPlateSplineExtended::Pointer tps,
                       vtkPolyData* poly, const std::string& outputFileName);
  double CalculateSpokeLength(PointType tail, PointType tip);
  void GetNeighborCells(vtkPolyData* mesh, int ptId, int newId, vtkCellArray* output, vtkPoints* morePts);
  void CompletePolyData(vtkPolyData *poly, vtkPolyData *output, vtkPolyData *medialMesh, bool isCrest = false);

private:

  vtkSlicerSRepInitializerLogic(const vtkSlicerSRepInitializerLogic&); // Not implemented
  void operator=(const vtkSlicerSRepInitializerLogic&); // Not implemented

private:
  int forwardCount = 0;
  int mRows = 5;
  int mCols = 9;
  std::string mOutputPath;
};

#endif
