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
#include <set>
#include <utility>

#include "vtkSlicerSkeletalRepresentationRefinerModuleLogicExport.h"
#include "vtkSlicerSkeletalRepresentationInterpolater.h"
#include "vtkImageData.h"
#include "vtkSmartPointer.h"

class vtkPolyData;
class vtkSpoke;
class vtkSrep;
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
  void Refine(double stepSize, double endCriterion, int maxIter);
  
  // Interpolate srep
  void InterpolateSrep(int interpolationLevel, std::string& srepFileName);
  
  // set weights for three items in the objective function
  void SetWeights(double wtImageMatch, double wtNormal, double wtSrad);
  
  // Description: Override operator (). Required by min_newuoa.
  // Parameter: @coeff: the pointer to coefficients
  double operator () (double *coeff);
  
  // This function returns the cost value result from defined objective function
  // given the current coeff array
  double EvaluateObjectiveFunction(double *coeff);
  
  // Generate anti-aliased signed distance map from surface mesh
  // Input: vtk file that contains target surface mesh
  // Output: image file that can be used in refinement
  void AntiAliasSignedDistanceMap(const std::string &meshFileName);

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
  
  // compute the difference between two vectors, factor can be used to compute center-difference
  void computeDiff(double *head, double *tail, double factor, double *output);
  
  // compute distance from implied boundary in signed distance map
  void computeDistance(vtkSpoke *theSpoke);
  
  // derivative of skeletal point
  void computeDerivative(std::vector<double> skeletalPoints, int r, int c, int nRows, int nCols, double *dXdu, double *dXdv);
  
  void convertSpokes2PolyData(std::vector<vtkSpoke*> input, vtkPolyData* output);

  // visualize model in MRMLScene
  void Visualize(vtkPolyData* model, const std::string &modelName, double r, double g, double b);
  
  void HideNodesByClass(const std::string &className);
  

private:
  std::string mImageFilePath;
  std::string mSrepFilePath;
  // weights in optimization algorithm
  double mWtImageMatch;
  double mWtNormalMatch;
  double mWtSrad;
  
  // output the first terms in object func can help to set weights
  bool mFirstCost = true;
  
  // store the data at the beginning of refinement
  int mNumRows;
  int mNumCols;
  vtkSrep* mSrep;
  std::vector<double> mCoeffArray;
  std::set<std::pair<double, double> > mInterpolatePositions;
  vtkSmartPointer<vtkImageData> mAntiAliasedImage = vtkSmartPointer<vtkImageData>::New();
private:

  vtkSlicerSkeletalRepresentationRefinerLogic(const vtkSlicerSkeletalRepresentationRefinerLogic&); // Not implemented
  void operator=(const vtkSlicerSkeletalRepresentationRefinerLogic&); // Not implemented
};

#endif
