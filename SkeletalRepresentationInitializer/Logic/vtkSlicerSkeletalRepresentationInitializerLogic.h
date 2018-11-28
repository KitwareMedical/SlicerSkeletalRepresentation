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
#include "vtkSlicerModuleLogic.h"

// MRML includes

// STD includes
#include <cstdlib>

#include "vtkSlicerSkeletalRepresentationInitializerModuleLogicExport.h"

class vtkPolyData;
class vtkPoints;

/// \ingroup Slicer_QtModules_ExtensionTemplate
class VTK_SLICER_SKELETALREPRESENTATIONINITIALIZER_MODULE_LOGIC_EXPORT vtkSlicerSkeletalRepresentationInitializerLogic :
  public vtkSlicerModuleLogic
{
public:

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
  int FlowSurfaceOneStep(double dt, double smooth_amount);

  // Select input mesh and render it in scene
  // input[filename]: whole path of vtk file
  int SetInputFileName(const std::string &filename);

  // Show fitting ellipsoid in 3D window
  // Can be called after one step flow or overall flow
  // if called after one step flow,  render the ellipsoid generated just now
  // otherwise render the ellipsoid at the end
  // input the forward flow deformed mesh, output radii ( double &rx, double &ry, double &rz)
  int ShowFittingEllipsoid(vtkPolyData* mesh, double &rx, double &ry, double &rz);

  // generate srep given an ellipsoid and expected rows and columns of medial sheet.
  int GenerateSrepForEllipsoid(vtkPolyData* mesh, int rows, int cols);

  int InklingFlow(const std::string &filename, double dt, double smooth_amount, int max_iter, int freq_output, double threshold);

  int BackwardFlow();

  // For the sake of completion of backward flow,
  // add this function to show what the process like.
  // This will be replaced by BackwardFlow later.
  int DummyBackwardFlow(std::string& output);
  int GenerateSrep(std::string& output);
  
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
  void AddPointToScene(double x, double y, double z, int glyphType, double r = 1, double g = 0, double b = 0);

private:

  vtkSlicerSkeletalRepresentationInitializerLogic(const vtkSlicerSkeletalRepresentationInitializerLogic&); // Not implemented
  void operator=(const vtkSlicerSkeletalRepresentationInitializerLogic&); // Not implemented

private:
  int forwardCount = 0;
};

#endif
