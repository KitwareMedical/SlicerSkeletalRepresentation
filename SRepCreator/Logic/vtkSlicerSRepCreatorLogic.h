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

#ifndef __vtkSlicerSRepCreatorLogic_h
#define __vtkSlicerSRepCreatorLogic_h

// Slicer includes
#include "vtkSlicerModuleLogic.h"

// MRML includes
#include "vtkMRMLModelNode.h"
#include "vtkMRMLEllipticalSRepNode.h"

// STD includes
#include <cstdlib>

// Eigen includes
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "vtkSlicerSRepCreatorModuleLogicExport.h"
#include <srep/EllipticalSRep.h>


/// \ingroup Slicer_QtModules_ExtensionTemplate
class VTK_SLICER_SREPCREATOR_MODULE_LOGIC_EXPORT vtkSlicerSRepCreatorLogic :
  public vtkSlicerModuleLogic
{
public:
  static vtkSlicerSRepCreatorLogic *New();
  vtkTypeMacro(vtkSlicerSRepCreatorLogic, vtkSlicerModuleLogic);
  void PrintSelf(ostream& os, vtkIndent indent);

  /// Creates a best fit ellipsoid for the input model
  bool RunForward(vtkMRMLModelNode* model, size_t numFoldPoints, size_t numStepsToFold, double dt, double smoothAmount, size_t maxIterations);

protected:
  vtkSlicerSRepCreatorLogic();
  virtual ~vtkSlicerSRepCreatorLogic();

  void SetMRMLSceneInternal(vtkMRMLScene* newScene) override;
private:
  vtkSlicerSRepCreatorLogic(const vtkSlicerSRepCreatorLogic&) = delete;
  void operator=(const vtkSlicerSRepCreatorLogic&) = delete;
  vtkSlicerSRepCreatorLogic(vtkSlicerSRepCreatorLogic&&) = delete;
  void operator=(vtkSlicerSRepCreatorLogic&&) = delete;

  struct EllipsoidParameters {
    Eigen::Vector3d radii;
    Eigen::Matrix3d rotation; ///< 3 by 3 rotation relative to deformed object
    Eigen::RowVector3d center;

    // not quite sure exactly what this and mry_o are.
    inline double mrx_o() const {
      // rz = radii(0)
      // ry = radii(1)
      // rx = radii(2)
      return (this->radii(2)*this->radii(2) - this->radii(0)*this->radii(0)) / this->radii(2);
    }

    inline double mry_o() const {
      // rz = radii(0)
      // ry = radii(1)
      // rx = radii(2)
      return (this->radii(1)*this->radii(1) - this->radii(0)*this->radii(0)) / this->radii(1);
    }
  };

  // It is on the user of this class to ensure that skeletalPoints, upSpokeBoundaryPoints,
  // and downSpokeBoundaryPoints are all the same size and crestSpokeBoundaryPoints and
  // crestSkeletalPoints are the same size and that those sizes match up with numFoldPoints
  // and numStepsToCrest. resize and the constructor can help, but there is no real protection.
  struct EigenSRep {
    Eigen::MatrixXd skeletalPoints; // up and down spoke skeletal points (they share the point)
    Eigen::MatrixXd upSpokeBoundaryPoints;
    Eigen::MatrixXd downSpokeBoundaryPoints;
    Eigen::MatrixXd crestSpokeBoundaryPoints;
    Eigen::MatrixXd crestSkeletalPoints;
    size_t numFoldPoints;
    size_t numStepsToCrest;

    EigenSRep()
      : numFoldPoints(0)
      , numStepsToCrest(0)
    {}
    EigenSRep(const EigenSRep&) = default;
    EigenSRep& operator=(const EigenSRep&) = default;
    EigenSRep(EigenSRep&&) = default;
    EigenSRep& operator=(EigenSRep&&) = default;
    ~EigenSRep() = default;

    inline EigenSRep(size_t numFoldPoints_, size_t numStepsToCrest_)
    {
      this->resize(numFoldPoints_, numStepsToCrest_);
    }

    inline void resize(size_t numFoldPoints, size_t numStepsToCrest) {
      const auto numSkeletalPoints = numFoldPoints*(numStepsToCrest+1);
      this->skeletalPoints.resize(numSkeletalPoints, 3);
      this->upSpokeBoundaryPoints.resize(numSkeletalPoints, 3);
      this->downSpokeBoundaryPoints.resize(numSkeletalPoints, 3);
      this->crestSpokeBoundaryPoints.resize(numFoldPoints, 3);
      this->crestSkeletalPoints.resize(numFoldPoints, 3);
      this->numFoldPoints = numFoldPoints;
      this->numStepsToCrest = numStepsToCrest;
    }
  };
  static std::unique_ptr<srep::EllipticalSRep> ConvertEigenSRepToEllipticalSRep(const EigenSRep& eigenSRep);

  std::string TempFolder();

  // Take surface mesh and "flows" it toward a more elliptical shape
  vtkSmartPointer<vtkPolyData> FlowSurfaceMesh(
    vtkMRMLModelNode* model,
    double dt,
    double smoothAmount,
    size_t maxIterations);

  static EllipsoidParameters CalculateBestFitEllipsoid(vtkPolyData& alreadyFlowedMesh);

  std::unique_ptr<srep::EllipticalSRep> GenerateSRep(
    const EllipsoidParameters& ellipsoid,
    size_t numFoldPoints,
    size_t numStepsToCrest);

  //return is <x part of sheet, y part of sheet>
  static std::pair<Eigen::MatrixXd, Eigen::MatrixXd> GenerateMedialSkeletalSheet(
    const EllipsoidParameters& ellipsoid,
    const size_t numFoldPoints,
    const size_t numStepsToCrest);

  static EigenSRep GenerateEigenSRep(
    const EllipsoidParameters& ellipsoid,
    const size_t numFoldPoints,
    const size_t numStepsToCrest);

  vtkMRMLModelNode* MakeEllipsoidModelNode(
    const EllipsoidParameters& ellipsoid,
    const std::string& name,
    bool visible = true,
    const double* color = nullptr);

  vtkMRMLModelNode* MakeModelNode(
    vtkPolyData* mesh,
    const std::string& name,
    bool visible = true,
    const double* color = nullptr);

  vtkMRMLEllipticalSRepNode* MakeEllipticalSRepNode(
    std::unique_ptr<srep::EllipticalSRep> srep,
    const std::string& name,
    bool visible = true);

  static constexpr double ellipse_scale = 0.9;
  static constexpr double eps = 1e-6;
  static constexpr double crestShift = 0.1; //10%
};

#endif
