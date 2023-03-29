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

// SRepWarper Logic includes
#include "vtkSlicerSRepWarperLogic.h"
#include "vtkSlicerSRepLogic.h"

// MRML includes
#include <vtkMRMLScene.h>

// VTK includes
#include <vtkCurvatures.h>
#include <vtkDoubleArray.h>
#include <vtkImageData.h>
#include <vtkImageMagnitude.h>
#include <vtkImageStencil.h>
#include <vtkImplicitPolyDataDistance.h>
#include <vtkIntArray.h>
#include <vtkMatrix4x4.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPointLocator.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataToImageStencil.h>

// ITK includes
#include <itkThinPlateSplineExtended.h>

// STD includes
#include <array>
#include <cassert>
#include <cstdlib>
#include <tuple>
#include <vector>

namespace {
  //---------------------------------------------------------------------------
  srep::Point3d ApplyTPS(const srep::Point3d& point, itkThinPlateSplineExtended::Pointer tps) {
    const auto transformed = tps->TransformPoint(point.AsArray());
    return srep::Point3d(transformed[0], transformed[1], transformed[2]);
  }

  //---------------------------------------------------------------------------
  vtkSmartPointer<vtkSRepSpoke> ApplyTPS(const vtkSRepSpoke& spoke, itkThinPlateSplineExtended::Pointer tps) {
    return vtkSRepSpoke::SmartCreate(ApplyTPS(spoke.GetSkeletalPoint(), tps), ApplyTPS(spoke.GetBoundaryPoint(), tps));
  }
  
  //---------------------------------------------------------------------------
  void ApplyTPSInPlace(vtkSRepSkeletalPoint& skeletalPoint, itkThinPlateSplineExtended::Pointer tps) {
    skeletalPoint.SetUpSpoke(ApplyTPS(*skeletalPoint.GetUpSpoke(), tps));
    skeletalPoint.SetDownSpoke(ApplyTPS(*skeletalPoint.GetDownSpoke(), tps));
    if (skeletalPoint.IsCrest()) {
      skeletalPoint.SetCrestSpoke(ApplyTPS(*skeletalPoint.GetCrestSpoke(), tps));
    }
  }

  //---------------------------------------------------------------------------
  void ApplyTPSInPlace(vtkEllipticalSRep& srep, itkThinPlateSplineExtended::Pointer tps) {
    using IndexType = vtkEllipticalSRep::IndexType;
    for (IndexType l = 0; l < srep.GetNumberOfLines(); ++l) {
      for (IndexType s = 0; s < srep.GetNumberOfSteps(); ++s) {
        ApplyTPSInPlace(*srep.GetSkeletalPoint(l, s), tps);
      }
    }
  }

} // namespace {}

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkSlicerSRepWarperLogic);

//----------------------------------------------------------------------------
vtkSlicerSRepWarperLogic::vtkSlicerSRepWarperLogic() = default;

//----------------------------------------------------------------------------
vtkSlicerSRepWarperLogic::~vtkSlicerSRepWarperLogic() = default;

//----------------------------------------------------------------------------
void vtkSlicerSRepWarperLogic::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

//---------------------------------------------------------------------------
void vtkSlicerSRepWarperLogic::ProgressCallback(double progress) {
  this->InvokeEvent(vtkCommand::ProgressEvent, &progress);
}

//---------------------------------------------------------------------------
vtkMRMLEllipticalSRepNode* vtkSlicerSRepWarperLogic::Run(
  vtkMRMLModelNode* sourceModelNode,
  vtkMRMLEllipticalSRepNode* sourceSRepNode,
  vtkMRMLModelNode* targetModelNode)
{
  vtkSmartPointer<vtkMRMLScene> scene = this->GetMRMLScene();
  if (!scene) {
    throw std::runtime_error("Can't add new vtkMRMLEllipticalSRepNode with null scene");
  }
  vtkMRMLEllipticalSRepNode* destination = vtkMRMLEllipticalSRepNode::SafeDownCast(scene->AddNewNodeByClass("vtkMRMLEllipticalSRepNode"));
  try {
    this->Run(sourceModelNode, sourceSRepNode, targetModelNode, destination);
    return destination;
  } catch (...) {
    scene->RemoveNode(destination);
    throw;
  }
}

//---------------------------------------------------------------------------
void vtkSlicerSRepWarperLogic::Run(
  vtkMRMLModelNode* sourceModelNode,
  vtkMRMLEllipticalSRepNode* sourceSRepNode,
  vtkMRMLModelNode* targetModelNode,
  vtkMRMLEllipticalSRepNode* targetSRepNode)
{
  try {
    if (!sourceModelNode) {
      throw std::invalid_argument("Cannot waro an SRep from a null model");
    }
    if (!sourceSRepNode || !sourceSRepNode->GetSRep() || sourceSRepNode->GetSRep()->IsEmpty()) {
      throw std::invalid_argument("Cannot warp an SRep with a null srep");
    }
    if (!targetModelNode) {
      throw std::invalid_argument("Cannot warp an SRep to a null model");
    }

    using TransformType = itkThinPlateSplineExtended;
    using PointType = itk::Point<double, 3>;
    using PointSetType = TransformType::PointSetType;
    using PointIdType = PointSetType::PointIdentifier;

    const auto polyDataPointToPointType = [](vtkPolyData& poly, unsigned int index) {
      PointType pt;
      double p[3];
      poly.GetPoint(index, p);
      pt[0] = p[0];
      pt[1] = p[1];
      pt[2] = p[2];
      return pt;
    };

    auto sourceModel = sourceModelNode->GetPolyData();
    auto sourceSRep = sourceSRepNode->GetEllipticalSRep();
    auto targetModel = targetModelNode->GetPolyData();
    auto targetSRep = sourceSRep->SmartClone();

    PointSetType::Pointer sourceLandMarks = PointSetType::New();
    PointSetType::Pointer targetLandMarks = PointSetType::New();
    PointSetType::PointsContainer::Pointer sourceLandMarkContainer
                = sourceLandMarks->GetPoints();
    PointSetType::PointsContainer::Pointer targetLandMarkContainer
                = targetLandMarks->GetPoints();

    // Read in the source points set
    for(unsigned int i = 0; i < sourceModel->GetNumberOfPoints(); ++i) {
        sourceLandMarkContainer->InsertElement(i, polyDataPointToPointType(*sourceModel, i));
    }

    // Read in the target points set
    for(unsigned int i = 0; i < targetModel->GetNumberOfPoints(); ++i) {
        targetLandMarkContainer->InsertElement(i, polyDataPointToPointType(*targetModel, i));
    }

    TransformType::Pointer tps = TransformType::New();
    tps->SetSourceLandmarks(sourceLandMarks);
    tps->SetTargetLandmarks(targetLandMarks);
    tps->ComputeWMatrix();

    ApplyTPSInPlace(*targetSRep, tps);

    targetSRepNode->SetEllipticalSRep(targetSRep);

  } catch (const std::exception& e) {
    vtkErrorMacro("Error running SRep Warper: " << e.what());
    throw;
  }
  catch (...) {
    vtkErrorMacro("Unknown error running SRep Warper");
    throw;
  }
}
