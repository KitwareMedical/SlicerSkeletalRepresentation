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

// .NAME vtkSlicerSkeletalRepresentationRefinerLogic - slicer logic class for volumes manipulation
// .SECTION Description
// This class manages the logic associated with reading, saving,
// and changing propertied of the volumes


#ifndef __vtkSlicerSkeletalRepresentationRefinerLogic_h
#define __vtkSlicerSkeletalRepresentationRefinerLogic_h

// Slicer includes
#include "vtkSlicerModuleLogic.h"

// MRML includes

// STD includes
#include <cstdlib>

#include "vtkSlicerSkeletalRepresentationRefinerModuleLogicExport.h"
#include "vtkSlicerSkeletalRepresentationInterpolater.h"

class vtkPolyData;
class vtkSpoke;
/// \ingroup Slicer_QtModules_ExtensionTemplate
class VTK_SLICER_SKELETALREPRESENTATIONREFINER_MODULE_LOGIC_EXPORT vtkSlicerSkeletalRepresentationRefinerLogic :
  public vtkSlicerModuleLogic
{
public:

  static vtkSlicerSkeletalRepresentationRefinerLogic *New();
  vtkTypeMacro(vtkSlicerSkeletalRepresentationRefinerLogic, vtkSlicerModuleLogic);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Select image file
  void SetImageFileName(const std::string &imageFilePath);

  // Select srep file
  void SetSrepFileName(const std::string &srepFilePath);

  // Start refinement
  void Refine();
  
  // Interpolate srep
  void InterpolateSrep(int interpolationLevel, const std::string& srepFileName);

protected:
  vtkSlicerSkeletalRepresentationRefinerLogic();
  virtual ~vtkSlicerSkeletalRepresentationRefinerLogic();

  virtual void SetMRMLSceneInternal(vtkMRMLScene* newScene);
  /// Register MRML Node classes to Scene. Gets called automatically when the MRMLScene is attached to this logic class.
  virtual void RegisterNodes();
  virtual void UpdateFromMRMLScene();
  virtual void OnMRMLSceneNodeAdded(vtkMRMLNode* node);
  virtual void OnMRMLSceneNodeRemoved(vtkMRMLNode* node);

private:
  // interpolate s-rep
  void Interpolate();

  // parse the s-rep
  // put the spoke length and direction into coeffArray
  void Parse(const std::string &modelFileName, std::vector<double> &coeffArray, 
             std::vector<double> &radii, std::vector<double> &dirs, std::vector<double> &skeletalPoints);
  
  // parse the header of s-rep including the rows and cols in the s-rep
  void ParseHeader(const std::string &headerFileName, int *nRows, int *nCols);
  
  void computeDiff(double *head, double *tail, double factor, double *output);
  
  // derivative of skeletal point
  void computeDerivative(std::vector<double> skeletalPoints, int r, int c, int nRows, int nCols, double *dXdu, double *dXdv);
  
  void convertSpokes2PolyData(std::vector<vtkSpoke*> input, vtkPolyData* output);

  // visualize model in MRMLScene
  void Visualize(vtkPolyData* model, const std::string &modelName, double r, double g, double b);

private:
  std::string mImageFilePath;
  std::string mSrepFilePath;
private:

  vtkSlicerSkeletalRepresentationRefinerLogic(const vtkSlicerSkeletalRepresentationRefinerLogic&); // Not implemented
  void operator=(const vtkSlicerSkeletalRepresentationRefinerLogic&); // Not implemented
};

#endif
