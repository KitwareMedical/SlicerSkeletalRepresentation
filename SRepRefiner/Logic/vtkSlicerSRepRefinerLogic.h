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

// .NAME vtkSlicerSRepRefinerLogic - slicer logic class for volumes manipulation
// .SECTION Description
// This class manages the logic associated with reading, saving,
// and changing propertied of the volumes


#ifndef __vtkSlicerSRepRefinerLogic_h
#define __vtkSlicerSRepRefinerLogic_h

// Slicer includes
#include "vtkSlicerModuleLogic.h"

// MRML includes

// STD includes
#include <cstdlib>
#include <set>
#include <utility>

#include "vtkSlicerSRepRefinerModuleLogicExport.h"
#include "SRepInterpolater.h"
#include "vtkImageData.h"
#include "vtkSmartPointer.h"
#include "itkImage.h"
#include "itkCovariantVector.h"

class vtkPolyData;
class vtkPoints;
class vtkCellArray;
class Spoke;
class SRep;
class vtkImplicitPolyDataDistance;
class vtkAppendPolyData;
/// \ingroup Slicer_QtModules_ExtensionTemplate
class VTK_SLICER_SREPREFINER_MODULE_LOGIC_EXPORT vtkSlicerSRepRefinerLogic :
  public vtkSlicerModuleLogic
{
public:
  typedef itk::Image<float, 3> RealImage;
  typedef itk::Image<itk::CovariantVector<float, 3>, 3> VectorImage;
  typedef std::pair<double, double> pairs;

  static vtkSlicerSRepRefinerLogic *New();
  vtkTypeMacro(vtkSlicerSRepRefinerLogic, vtkSlicerModuleLogic)
  void PrintSelf(ostream& os, vtkIndent indent) override;

  // Select image file
  void SetImageFileName(const std::string &imageFilePath);

  // Select srep file
  void SetSrepFileName(const std::string &srepFilePath);

  // Select output path
  void SetOutputPath(const std::string &outputPath);

  // Start refinement
  // Input: stepSize is step in NEWUOA
  // Input: endCriterion is tol in NEWUOA
  // Input: maxIter is the max number of iteration of NEWUOA
  // Input: interpolationLevel is the density when computing image match term
  void Refine(double stepSize, double endCriterion, int maxIter, int interpolationLevel);

  // Interpolate srep
  void InterpolateSrep(int interpolationLevel, std::string& srepFileName);
  void InterpolateSrep(int interpolationLevel, int nRows, int nCols,
                       std::string& srepFileName, std::string& crestFileName, std::vector<Spoke*> &interpolatedSpokes);
  void InterpolateSrep(int interpolationLevel, int nRows, int nCols,
                       std::string& srepFileName, std::string& crestFileName,
                       std::vector<Spoke*> &interpolatedSpokes, std::vector<Spoke*> &interpolatedCrestSpokes);

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

  // Compute transformation matrix from srep to image coordinate system, namely, unit cube cs.
  void TransformSrep(const std::string &headerFile);

  // show wired frame of implied boundary
  void ShowImpliedBoundary(int interpolationLevel, const std::string& srepFileName);

  // command line interface for refinement
  void CLIRefine(const std::string &srepFileName, const std::string &imgFileName, const std::string &outputPath,
                 double stepSize = 0.01, double endCriterion = 0.001, int maxIter = 2000,
                 double wtImg = 0.004, double wtNormal = 20, double wtSrad = 50, int interpolationLevel = 3);


protected:
  vtkSlicerSRepRefinerLogic();
  ~vtkSlicerSRepRefinerLogic() override;

  void SetMRMLSceneInternal(vtkMRMLScene* newScene) override;
  /// Register MRML Node classes to Scene. Gets called automatically when the MRMLScene is attached to this logic class.
  void RegisterNodes() override;
  void UpdateFromMRMLScene() override;
  void OnMRMLSceneNodeAdded(vtkMRMLNode* node) override;
  void OnMRMLSceneNodeRemoved(vtkMRMLNode* node) override;

protected:
  // interpolate s-rep
  void Interpolate();
  void VisualizeHeatMap(vtkPolyData* inputMesh);
  // parse the s-rep
  // put the spoke length and direction into coeffArray
  void Parse(const std::string &modelFileName, std::vector<double> &coeffArray,
             std::vector<double> &radii, std::vector<double> &dirs, std::vector<double> &skeletalPoints);

  // parse the header of s-rep including the rows and cols in the s-rep
  void ParseHeader(const std::string &headerFileName, int *nRows, int *nCols, double *shift,
                   std::string* upFileName, std::string* downFileName, std::string* crestFileName);

  // update header file after refinement
  void UpdateHeader(const std::string &headerFileName, const std::string &outputFileName, std::string *newHeaderFileName);

  // compute the difference between two vectors, factor can be used to compute center-difference
  void ComputeDiff(double *head, double *tail, double factor, double *output);

  // compute distance from an implied boundary to the target boundary
  double ComputeDistance(Spoke *theSpoke, double *normalMatch);

  // derivative of skeletal point
  void ComputeDerivative(std::vector<double> skeletalPoints, int r, int c, int nRows, int nCols, double *dXdu, double *dXdv);

  void ConvertSpokes2PolyData(std::vector<Spoke*> input, vtkPolyData* output);

  void SaveSpokes2Vtp(std::vector<Spoke*> input, const std::string &path);

  void TransSpokes2PolyData(std::vector<Spoke *>input, vtkPolyData *output);

  // show points as fiducial markups
  void VisualizePoints(vtkPoints* input);

  // visualize model in MRMLScene
  void Visualize(vtkPolyData* model, const std::string &modelName,
                 double r, double g, double b, bool isVisible = true);

  void HideNodesByClass(const std::string &className);

  // transform model cs to unit cube cs then to image cs
  // Output should be formed as 4x4 matrix
  void TransformSrep2ImageCS(SRep *input, double output[][4]);

  // Get all interpolated as well as primary spokes including top, bottom and down
  void ConnectImpliedBoundaryPts(int interpolationLevel, int nRows, int nCols, const std::string &srepFileName,
                                 std::vector<Spoke *> &borderSpokes,
                                 std::vector<Spoke*>& interpolated,
                                 std::vector<Spoke *> &repeatedInterps,
                                 std::vector<Spoke*>& primary,
                                 vtkPolyData* impliedPolyData);

  // Interpolate crest region and connect them with quads
  void ConnectImpliedCrest(int interpolationLevel, int nRows, int nCols,
                           const std::string &crest, std::vector<Spoke*> &upSpokes,std::vector<Spoke*> &downSpokes,
                           vtkPolyData* crestPoly);
  void ConnectImpliedCrest(int interpolationLevel,
                           std::vector<Spoke *> upInterpSpokes,
                           std::vector<Spoke *> downInterpSpokes,
                           std::vector<Spoke *> crestInterpSpokes,
                           vtkPolyData* crestPoly);

  // connect fold curve macro
  void ConnectFoldCurve(const std::vector<Spoke *>& edgeSpokes, vtkPoints *foldCurvePts, vtkCellArray *foldCurveCell);

  // e.g. Refine up spokes saved in upFileName
  // return refined collection of spokes
  std::vector<Spoke*>& RefinePartOfSpokes(const std::string& srepFileName, double stepSize, double endCriterion, int maxIter);

  void RefineCrestSpokes(const std::string &crest,
                         double stepSize, double endCriterion, int maxIter);
  // compute total distance of left top spoke to the quad
  double TotalDistOfLeftTopSpoke(SRep* input, double u, double v, int r, int c, double *normalMatch);

  // compute total distance of right top spoke to the quad
  double TotalDistOfRightTopSpoke(SRep* input, double u, double v, int r, int c, double *normalMatch);

  // compute total distance of left bottom spoke to the quad
  double TotalDistOfLeftBotSpoke(SRep* input, double u, double v, int r, int c, double *normalMatch);

  // compute total distance of Right bottom spoke to the quad
  double TotalDistOfRightBotSpoke(SRep* input, double u, double v, int r, int c, double *normalMatch);

  // compute rSrad penalty
  double ComputeRSradPenalty(SRep* input);

  void FindTopLeftNeigbors(int r, int c,
                           SRep* input,
                           std::vector<Spoke *> &neighborU,
                           std::vector<Spoke *> &neighborV);

  void FindTopRightNeigbors(int r, int c,
                           SRep* input,
                           std::vector<Spoke *> &neighborU,
                           std::vector<Spoke *> &neighborV);

  void FindBotLeftNeigbors(int r, int c,
                           SRep* input,
                           std::vector<Spoke *> &neighborU,
                           std::vector<Spoke *> &neighborV);

  void FindBotRightNeigbors(int r, int c,
                           SRep* input,
                           std::vector<Spoke *> &neighborU,
                           std::vector<Spoke *> &neighborV);
  void ParseCrest(const std::string &crestFileName, std::vector<Spoke*> &crestSpokes);

  void InterpolateCrest(std::vector<Spoke*> &crestSpoke, std::vector<Spoke*> &interiorSpokes,
                        int interpolationLevel,
                        int nRows, int nCols,
                        std::vector<Spoke*> &crest, std::vector<Spoke*> &interior);

  void ReorderCrestSpokes(int nRows, int nCols, std::vector<Spoke*> &input, std::vector<Spoke*> &output);

  void OptimizeCrestSpokeLength(vtkImplicitPolyDataDistance *distanceFunction, Spoke *targetSpoke, double stepSize, int maxIter);

  void Transform2ImageCS(double *ptInput, int *ptOutput);

  void RetileMesh(vtkPolyData* targetMesh, vtkPolyData* impliedMesh, vtkPolyData* retiledMesh);
  void InterpolateCrest(int interpolationLevel,
                        std::vector<Spoke *> upInterpSpokes,
                        std::vector<Spoke *> downInterpSpokes,
                        std::vector<Spoke *> crestInterpSpokes,
                        vtkPolyData* crestPoly);

protected:
  std::vector<std::pair<double, double> > mInterpolatePositions;
  //vtkSmartPointer<vtkImageData> mAntiAliasedImage = vtkSmartPointer<vtkImageData>::New();


private:
  // when apply this transformation: [x, y, z, 1] * mTransformationMat
  double mTransformationMat[4][4]; // homogeneous matrix transfrom from srep to unit cube cs.

  RealImage::Pointer mAntiAliasedImage = RealImage::New();
  VectorImage::Pointer mGradDistImage = VectorImage::New();
  // weights in optimization algorithm
  double mWtImageMatch;
  double mWtNormalMatch;
  double mWtSrad;

  std::string mTargetMeshFilePath;
  std::string mSrepFilePath;
  std::string mOutputPath;
  std::vector<Spoke*> mUpSpokes;
  std::vector<Spoke*> mDownSpokes;
  int mInterpolationLevel;
  int mNumRows;
  int mNumCols;
  // output the first terms in object func can help to set weights
  bool mFirstCost = true;

  // store the data at the beginning of refinement
  SRep* mSrep;
  std::vector<double> mCoeffArray;
private:

  vtkSlicerSRepRefinerLogic(const vtkSlicerSRepRefinerLogic&); // Not implemented
  void operator=(const vtkSlicerSRepRefinerLogic&); // Not implemented
};

#endif
