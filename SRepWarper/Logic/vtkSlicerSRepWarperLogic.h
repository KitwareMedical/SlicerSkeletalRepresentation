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

// .NAME vtkSlicerSRepWarperLogic - slicer logic class for volumes manipulation
// .SECTION Description
// This class manages the logic associated with reading, saving,
// and changing propertied of the volumes


#ifndef __vtkSlicerSRepWarperLogic_h
#define __vtkSlicerSRepWarperLogic_h

// Slicer includes
#include "vtkSlicerModuleLogic.h"

// MRML includes
#include <vtkMRMLModelNode.h>
#include <vtkMRMLEllipticalSRepNode.h>

#include "vtkSlicerSRepWarperModuleLogicExport.h"

/// \ingroup Slicer_QtModules_ExtensionTemplate
class VTK_SLICER_SREPWARPER_MODULE_LOGIC_EXPORT vtkSlicerSRepWarperLogic :
  public vtkSlicerModuleLogic
{
public:

  static vtkSlicerSRepWarperLogic *New();
  vtkTypeMacro(vtkSlicerSRepWarperLogic, vtkSlicerModuleLogic);
  void PrintSelf(ostream& os, vtkIndent indent);

  /// @{
  /// Warps the given SRep to a target model based on correspondences with the source model.
  /// \param sourceModelNode The source model
  /// \param sourceSRepNode The SRep fit to the source model
  /// \param targetModelNode The model to warp the SRep to
  /// \returns The warped SRep.
  vtkMRMLEllipticalSRepNode* Run(
    vtkMRMLModelNode* sourceModelNode,
    vtkMRMLEllipticalSRepNode* sourceSRepNode,
    vtkMRMLModelNode* targetModelNode
    );
  void Run(
    vtkMRMLModelNode* sourceModelNode,
    vtkMRMLEllipticalSRepNode* sourceSRepNode,
    vtkMRMLModelNode* targetModelNode,
    vtkMRMLEllipticalSRepNode* outputSRepNode);
  /// @}

protected:
  vtkSlicerSRepWarperLogic();
  virtual ~vtkSlicerSRepWarperLogic();
private:
  void ProgressCallback(double progress);

  vtkSlicerSRepWarperLogic(const vtkSlicerSRepWarperLogic&); // Not implemented
  void operator=(const vtkSlicerSRepWarperLogic&); // Not implemented
};

#endif
