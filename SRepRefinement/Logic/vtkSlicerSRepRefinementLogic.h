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

// .NAME vtkSlicerSRepRefinementLogic - slicer logic class for volumes manipulation
// .SECTION Description
// This class manages the logic associated with reading, saving,
// and changing propertied of the volumes


#ifndef __vtkSlicerSRepRefinementLogic_h
#define __vtkSlicerSRepRefinementLogic_h

// Slicer includes
#include "vtkSlicerModuleLogic.h"

// MRML includes
#include <vtkMRMLModelNode.h>
#include <vtkMRMLEllipticalSRepNode.h>

#include "vtkSlicerSRepRefinementModuleLogicExport.h"

/// \ingroup Slicer_QtModules_ExtensionTemplate
class VTK_SLICER_SREPREFINEMENT_MODULE_LOGIC_EXPORT vtkSlicerSRepRefinementLogic :
  public vtkSlicerModuleLogic
{
public:

  static vtkSlicerSRepRefinementLogic *New();
  vtkTypeMacro(vtkSlicerSRepRefinementLogic, vtkSlicerModuleLogic);
  void PrintSelf(ostream& os, vtkIndent indent);

  /// @{
  /// Refines the given SRep to a Model.
  /// \param model The model to refine to.
  /// \param srep The srep to refine to the model.
  /// \param initialRegionSize Initial size of the trust region radius. Must be positive.
  /// \param finalRegionSize Final size of the trust region radius. Must be positive and less than or equal to initialRegionSize.
  ///        This serves as an ending criteria.
  /// \param maxIterations The maximum iterations to run. Less than this number of iterations will be run if the finalRegionSize is reached.
  /// \param interpolationLevel The level of interpolation to use when determining the "fit" to the model.
  /// \param L0Weight The weight to put on the L0 parameter. The L0 parameter is the squared distance
  ///        between the SRep implied boundary and the model boundary.
  /// \param L1Weight The weight to put on the L1 parameter. The L1 parameter is the deviation of spokes
  ///        from being perpendicular to the boundary.
  /// \param L2Weight The weight to put on the L2 parameter. The L2 parameter is the geometric illegality
  ///        of spokes. This parameter is intended to prevent spokes from crossing each other.
  /// \returns The refined SRep.
  vtkMRMLEllipticalSRepNode* Run(
    vtkMRMLModelNode* model,
    vtkMRMLEllipticalSRepNode* srep,
    double initialRegionSize,
    double finalRegionSize,
    int maxIterations,
    int interpolationLevel,
    double L0Weight,
    double L1Weight,
    double L2Weight);
  void Run(
    vtkMRMLModelNode* model,
    vtkMRMLEllipticalSRepNode* srep,
    double initialRegionSize,
    double finalRegionSize,
    int maxIterations,
    int interpolationLevel,
    double L0Weight,
    double L1Weight,
    double L2Weight,
    vtkMRMLEllipticalSRepNode* destination);
  /// @}

protected:
  vtkSlicerSRepRefinementLogic();
  virtual ~vtkSlicerSRepRefinementLogic();
private:
  void ProgressCallback(double progress);

  vtkSlicerSRepRefinementLogic(const vtkSlicerSRepRefinementLogic&); // Not implemented
  void operator=(const vtkSlicerSRepRefinementLogic&); // Not implemented
};

#endif
