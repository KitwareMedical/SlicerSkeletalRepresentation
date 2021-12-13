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

// Logic includes
#include "vtkSlicerSRepCreatorLogic.h"
#include "vtkSlicerSRepLogic.h"

// MRML includes
#include <vtkMRMLScene.h>
#include <vtkMRMLDisplayNode.h>

// VTK includes
#include <vtkCellLocator.h>
#include <vtkCurvatures.h>
#include <vtkDoubleArray.h>
#include <vtkGenericCell.h>
#include <vtkIntArray.h>
#include <vtkMassProperties.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkParametricEllipsoid.h>
#include <vtkParametricFunctionSource.h>
#include <vtkPointData.h>
#include <vtkPolyDataNormals.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkWindowedSincPolyDataFilter.h>

#include <vtksys/SystemTools.hxx>

#include <itkThinPlateSplineExtended.h>
#include <itkPointSet.h>

// STD includes
#include <cassert>
#include <sstream>

namespace {
  //---------------------------------------------------------------------------
  Eigen::MatrixXd ConvertVTKPointsToEigen(vtkPoints& points) {
    Eigen::MatrixXd point_matrix(points.GetNumberOfPoints(), 3);
    for(int i = 0; i < points.GetNumberOfPoints(); ++i) {
      double p[3];
      points.GetPoint(i, p);
      point_matrix.row(i) << p[0], p[1], p[2];
    }
    return point_matrix;
  }

  //---------------------------------------------------------------------------
  vtkSmartPointer<vtkPoints> ConvertEigenToVTKPoints(const Eigen::MatrixXd& matrix) {
    if (matrix.cols() != 3) {
      throw std::invalid_argument("Expected 3 columns to convert matrix to vtkPoints");
    }
    auto result = vtkSmartPointer<vtkPoints>::New();
    for (int row = 0; row < matrix.rows(); ++row) {
      const double p[] = {matrix(row, 0), matrix(row, 1), matrix(row, 2)};
      result->InsertNextPoint(p);
    }
    return result;
  }

  //---------------------------------------------------------------------------
  double VolumeOfEllipsoid(const double radii1, const double radii2, const double radii3) {
    return 4 / 3.0 * vtkMath::Pi() * radii1 * radii2 * radii3;
  }

  //---------------------------------------------------------------------------
  srep::Point3d PointFromEigen(const Eigen::Vector3d& v) {
    return srep::Point3d(v(0), v(1), v(2));
  }

  //---------------------------------------------------------------------------
  srep::Point3d ApplyTPS(const srep::Point3d& point, itkThinPlateSplineExtended::Pointer tps) {
    const auto transformed = tps->TransformPoint(point.AsArray());
    return srep::Point3d(transformed[0], transformed[1], transformed[2]);
  }

  //---------------------------------------------------------------------------
  srep::Spoke ApplyTPS(const srep::Spoke& spoke, itkThinPlateSplineExtended::Pointer tps) {
    return srep::Spoke(ApplyTPS(spoke.GetSkeletalPoint(), tps), ApplyTPS(spoke.GetBoundaryPoint(), tps));
  }

  //---------------------------------------------------------------------------
  void ApplyTPSInPlace(srep::SkeletalPoint& skeletalPoint, itkThinPlateSplineExtended::Pointer tps) {
    skeletalPoint.SetUpSpoke(ApplyTPS(skeletalPoint.GetUpSpoke(), tps));
    skeletalPoint.SetDownSpoke(ApplyTPS(skeletalPoint.GetDownSpoke(), tps));
    if (skeletalPoint.IsCrest()) {
      skeletalPoint.SetCrestSpoke(ApplyTPS(skeletalPoint.GetCrestSpoke(), tps));
    }
  }

  //---------------------------------------------------------------------------
  void ApplyTPSInPlace(srep::EllipticalSRep::UnrolledEllipticalGrid& grid, itkThinPlateSplineExtended::Pointer tps) {
    for (size_t i = 0; i < grid.size(); ++i) {
      for (size_t j = 0; j < grid[i].size(); ++j) {
        ApplyTPSInPlace(grid[i][j], tps);
      }
    }
  }

} // namespace {}


//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkSlicerSRepCreatorLogic);

//----------------------------------------------------------------------------
vtkSlicerSRepCreatorLogic::vtkSlicerSRepCreatorLogic()
  : ActualForwardIterations(0)
  , SRepNodeId()
  , ModelName()
{}

//----------------------------------------------------------------------------
vtkSlicerSRepCreatorLogic::~vtkSlicerSRepCreatorLogic() = default;

//----------------------------------------------------------------------------
void vtkSlicerSRepCreatorLogic::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

//---------------------------------------------------------------------------
void vtkSlicerSRepCreatorLogic::SetMRMLSceneInternal(vtkMRMLScene * newScene)
{
  vtkNew<vtkIntArray> events;
  events->InsertNextValue(vtkMRMLScene::NodeAddedEvent);
  events->InsertNextValue(vtkMRMLScene::NodeRemovedEvent);
  events->InsertNextValue(vtkMRMLScene::EndBatchProcessEvent);
  this->SetAndObserveMRMLSceneEventsInternal(newScene, events.GetPointer());
}

//---------------------------------------------------------------------------
std::string vtkSlicerSRepCreatorLogic::TempFolder() {
  std::stringstream ssTempFolder;
  // putting this pointer in the folder name so multiple logics can exist without interfering with each other
  ssTempFolder << this->GetApplicationLogic()->GetTemporaryPath() << "/SRepCreator-" << this;
  const auto tempFolder = ssTempFolder.str();
  if (!vtksys::SystemTools::FileExists(tempFolder, false)) {
    if (!vtksys::SystemTools::MakeDirectory(tempFolder)) {
      std::cerr << "Failed to create folder : " << tempFolder << std::endl;
      return "";
    }
  }
  return tempFolder;
}

//---------------------------------------------------------------------------
vtkSlicerSRepCreatorLogic::EllipsoidParameters vtkSlicerSRepCreatorLogic::FlowSurfaceMeshToEllipsoid(
  vtkMRMLModelNode* model,
  const double dt,
  const double smoothAmount,
  const size_t maxIterations,
  const size_t outputEveryNumIterations)
{
  auto flowedMesh = this->FlowSurfaceMesh(model, dt, smoothAmount, maxIterations, outputEveryNumIterations);
  if (!flowedMesh) {
    throw std::runtime_error("Error creating flowed mesh");
  }

  const auto ellipsoidParameters = CalculateBestFitEllipsoid(*flowedMesh);
  auto ellipsoidalMesh = this->SnapFlowedMeshToEllipsoid(*flowedMesh, ellipsoidParameters);

  vtkNew<vtkPolyDataWriter> writer;
  writer->SetFileName(this->ForwardIterationFilename(maxIterations+1).c_str());
  writer->SetInputData(ellipsoidalMesh);
  writer->Update();
  ++this->ActualForwardIterations;

  if (outputEveryNumIterations != 0) {
    this->MakeModelNode(ellipsoidalMesh,
      model->GetName() + std::string("-final-flowed-ellipsoidal-mesh-") + std::to_string(maxIterations+1),
      true, model->GetDisplayNode()->GetColor());
  }

  return ellipsoidParameters;
}

//---------------------------------------------------------------------------
vtkSmartPointer<vtkPolyData> vtkSlicerSRepCreatorLogic::FlowSurfaceMesh(
  vtkMRMLModelNode* model,
  const double dt,
  const double smoothAmount,
  const size_t maxIterations,
  const size_t outputEveryNumIterations)
{
  if (!model) {
    return nullptr;
  }

  auto mesh = vtkSmartPointer<vtkPolyData>::New();
  mesh->DeepCopy(model->GetMesh());
  if (!mesh) {
    return nullptr;
  }

  //create a temp folder to store stuff for backwards flow in
  const auto tempFolder = this->TempFolder();
  if (tempFolder.empty()) {
    return nullptr;
  }

  //TODO: delete if don't need volume
  // auto massFilter = vtkSmartPointer<vtkMassProperties>::New();
  // massFilter->SetInputData(mesh);
  // massFilter->Update();

  // const double originalVolume = massFilter->GetVolume();

  const bool smoothing = smoothAmount > 0;

  vtkNew<vtkPolyDataNormals> normalFilter;
  normalFilter->SplittingOff();
  normalFilter->ComputeCellNormalsOff();
  normalFilter->ComputePointNormalsOn();
  if (!smoothing) {
    normalFilter->SetInputData(mesh); //mesh only changes if we are smoothing
  }

  vtkNew<vtkCurvatures> curvatureFilter;
  curvatureFilter->SetInputConnection(normalFilter->GetOutputPort());
  curvatureFilter->SetCurvatureTypeToMean();

  vtkNew<vtkPolyDataWriter> writer;

  for (size_t i = 0; i < maxIterations; ++i) {
    vtkSmartPointer<vtkWindowedSincPolyDataFilter> smoothFilter;
    if (smoothing) {
      //there is something weird about this filter where it doesn't work
      //if you don't make it new every time
      smoothFilter = vtkSmartPointer<vtkWindowedSincPolyDataFilter>::New();
      smoothFilter->SetPassBand(smoothAmount);
      smoothFilter->NonManifoldSmoothingOn();
      smoothFilter->NormalizeCoordinatesOn();
      smoothFilter->SetNumberOfIterations(20);
      smoothFilter->FeatureEdgeSmoothingOff();
      smoothFilter->BoundarySmoothingOff();

      smoothFilter->SetInputData(mesh);
      normalFilter->SetInputConnection(smoothFilter->GetOutputPort());
    }
    curvatureFilter->Update();
    if (smoothing) {
      mesh = smoothFilter->GetOutput();
    }

    auto curvature = vtkDoubleArray::SafeDownCast(curvatureFilter->GetOutput()->GetPointData()->GetArray("Mean_Curvature"));
    auto normals = normalFilter->GetOutput()->GetPointData()->GetNormals();
    if (!curvature || !normals) {
      return nullptr;
    }

    // perform the flow
    auto points = mesh->GetPoints();
    for (int i = 0; i < points->GetNumberOfPoints(); ++i) {
      double p[3];
      points->GetPoint(i, p);
      const double* n = normals->GetTuple(i);
      const auto h = curvature->GetValue(i);
      for (int j = 0; j < 3; ++j) {
        p[j] -= dt * h * n[j];
      }
      points->SetPoint(i, p);
    }
    points->Modified();

    // TODO: do we need the mass stuff?
    writer->SetFileName(this->ForwardIterationFilename(i+1).c_str());
    writer->SetInputData(mesh);
    writer->Update();

    if (outputEveryNumIterations != 0 && i % outputEveryNumIterations == 0) {
      this->MakeModelNode(mesh,
        model->GetName() + std::string("-forwardflow-") + std::to_string(i),
        true, model->GetDisplayNode()->GetColor());
    }
  }
  this->ActualForwardIterations = maxIterations;

  if (outputEveryNumIterations != 0) {
    this->MakeModelNode(mesh,
      model->GetName() + std::string("-final-flowed-mesh-") + std::to_string(maxIterations),
      true, model->GetDisplayNode()->GetColor());
  }

  return mesh;
}

//---------------------------------------------------------------------------
std::string vtkSlicerSRepCreatorLogic::ForwardIterationFilename(const long iteration) {
  return this->TempFolder() + "/" + std::to_string(iteration) + ".vtk";
}

//---------------------------------------------------------------------------
vtkSmartPointer<vtkPolyData> vtkSlicerSRepCreatorLogic::SnapFlowedMeshToEllipsoid(vtkPolyData& alreadyFlowedMesh, const EllipsoidParameters& ellipsoid) {
  auto ellipsoidPolyData = vtkSlicerSRepCreatorLogic::MakeEllipsoidPolyData(ellipsoid);
  if (!ellipsoidPolyData) {
    throw std::runtime_error("Error creating ellipsoid PolyData");
  }

  vtkNew<vtkCellLocator> cellLocator;
  cellLocator->SetDataSet(ellipsoidPolyData);
  cellLocator->BuildLocator();

  double point[3];
  double closestPoint[3];
  vtkNew<vtkGenericCell> cell;
  vtkIdType cellId;
  int subId;
  double distanceSquared;

  vtkNew<vtkPoints> snappedMeshPoints;
  for (auto i = 0; i < alreadyFlowedMesh.GetNumberOfPoints(); ++i) {
    alreadyFlowedMesh.GetPoint(i, point);
    cellLocator->FindClosestPoint(point, closestPoint, cell, cellId, subId, distanceSquared);
    snappedMeshPoints->InsertNextPoint(closestPoint);
  }

  vtkNew<vtkPolyData> snappedMeshPolyData;
  snappedMeshPolyData->SetPoints(snappedMeshPoints);
  snappedMeshPolyData->SetPolys(alreadyFlowedMesh.GetPolys());
  snappedMeshPolyData->Modified();
  return snappedMeshPolyData;
}

//---------------------------------------------------------------------------
vtkSlicerSRepCreatorLogic::EllipsoidParameters
vtkSlicerSRepCreatorLogic::CalculateBestFitEllipsoid(vtkPolyData& alreadyFlowedMesh) {
  EllipsoidParameters result;

  const auto point_matrix = ConvertVTKPointsToEigen(*(alreadyFlowedMesh.GetPoints()));
  result.center = point_matrix.colwise().mean();
  Eigen::MatrixXd centered_point_mat = point_matrix - result.center.replicate(point_matrix.rows(), 1);
  Eigen::MatrixXd point_matrix_transposed = centered_point_mat.transpose();
  Eigen::Matrix3d second_moment = point_matrix_transposed * centered_point_mat;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(second_moment);
  result.radii = es.eigenvalues();
  result.radii(0) = sqrt(result.radii(0));
  result.radii(1) = sqrt(result.radii(1));
  result.radii(2) = sqrt(result.radii(2));

  const auto ellipsoid_volume = VolumeOfEllipsoid(result.radii(0), result.radii(1), result.radii(2));

  vtkNew<vtkMassProperties> mass;
  mass->SetInputData(&alreadyFlowedMesh);
  mass->Update();

  double volume_factor = pow(mass->GetVolume() / ellipsoid_volume, 1.0 / 3.0);
  result.radii(0) *= volume_factor;
  result.radii(1) *= volume_factor;
  result.radii(2) *= volume_factor;

  result.rotation = es.eigenvectors();
  return result;
}

//---------------------------------------------------------------------------
vtkSmartPointer<vtkPolyData> vtkSlicerSRepCreatorLogic::MakeEllipsoidPolyData(const EllipsoidParameters& ellipsoid) {
  vtkNew<vtkParametricEllipsoid> parametricEllipsoid;
  parametricEllipsoid->SetXRadius(ellipsoid.radii(0));
  parametricEllipsoid->SetYRadius(ellipsoid.radii(1));
  parametricEllipsoid->SetZRadius(ellipsoid.radii(2));

  vtkNew<vtkParametricFunctionSource> parametric_function;
  parametric_function->SetParametricFunction(parametricEllipsoid);
  parametric_function->SetUResolution(30);
  parametric_function->SetVResolution(30);
  parametric_function->Update();

  vtkSmartPointer<vtkPolyData> ellipsoid_polydata = parametric_function->GetOutput();

  auto ellipsoid_points_matrix = ConvertVTKPointsToEigen(*(ellipsoid_polydata->GetPoints()));
  // rotate the points
  Eigen::MatrixXd rotated_ellipsoid_points = ellipsoid.rotation * (ellipsoid_points_matrix.transpose());
  rotated_ellipsoid_points.transposeInPlace(); // n x 3
  // translate the points
  const Eigen::MatrixXd translated_points = rotated_ellipsoid_points + ellipsoid.center.replicate(rotated_ellipsoid_points.rows(),1);

  auto bestFitEllipsoidPoints = ConvertEigenToVTKPoints(translated_points);
  vtkNew<vtkPolyData> bestFitEllipsoidPolyData;
  bestFitEllipsoidPolyData->SetPoints(bestFitEllipsoidPoints);
  bestFitEllipsoidPolyData->SetPolys(ellipsoid_polydata->GetPolys());
  bestFitEllipsoidPolyData->Modified();

  return bestFitEllipsoidPolyData;
}

//---------------------------------------------------------------------------
vtkMRMLModelNode* vtkSlicerSRepCreatorLogic::MakeEllipsoidModelNode(
  const EllipsoidParameters& ellipsoid,
  const std::string& name,
  const bool visible,
  const double* color)
{
  auto ellipsoidPolyData = this->MakeEllipsoidPolyData(ellipsoid);
  return MakeModelNode(ellipsoidPolyData, name, visible, color);
}

//---------------------------------------------------------------------------
vtkMRMLModelNode* vtkSlicerSRepCreatorLogic::MakeModelNode(
  vtkPolyData* mesh,
  const std::string& name,
  const bool visible,
  const double* color)
{
  auto scene = this->GetMRMLScene();
  if (!scene) {
    return nullptr;
  }

  auto model = vtkMRMLModelNode::SafeDownCast(scene->AddNewNodeByClass("vtkMRMLModelNode"));
  if (!model) {
    return nullptr;
  }
  model->SetScene(scene);
  model->SetName(name.c_str());
  model->SetAndObserveMesh(mesh);
  model->CreateDefaultDisplayNodes();

  auto displayNode = model->GetDisplayNode();
  if (color) {
    // SetColor isn't going to change the value of the color param
    // and we want out interface to be const correct where possible,
    // hence the const_cast
    displayNode->SetColor(const_cast<double*>(color));
  }
  displayNode->SetBackfaceCulling(0);
  displayNode->SetRepresentation(vtkMRMLDisplayNode::WireframeRepresentation);
  displayNode->SetVisibility(visible);

  return model;
}

//---------------------------------------------------------------------------
vtkMRMLEllipticalSRepNode* vtkSlicerSRepCreatorLogic::MakeEllipticalSRepNode(
    std::unique_ptr<srep::EllipticalSRep> srep,
    const std::string& name,
    const bool visible)
{
  auto scene = this->GetMRMLScene();
  if (!scene) {
    return nullptr;
  }

  const auto srepNodeId = vtkSmartPointer<vtkSlicerSRepLogic>::New()->AddNewEllipticalSRepNode(name, scene);
  if (srepNodeId.empty()) {
    return nullptr;
  }
  auto srepNode = vtkMRMLEllipticalSRepNode::SafeDownCast(scene->GetNodeByID(srepNodeId));
  srepNode->SetEllipticalSRep(std::move(srep));
  srepNode->GetDisplayNode()->SetVisibility(visible);
  return srepNode;
}

//---------------------------------------------------------------------------
std::pair<Eigen::MatrixXd, Eigen::MatrixXd> vtkSlicerSRepCreatorLogic::GenerateMedialSkeletalSheet(
  const EllipsoidParameters& ellipsoid,
  const size_t numFoldPoints,
  const size_t numStepsToCrest)
{
  const double mra = ellipsoid.mrx_o() * vtkSlicerSRepCreatorLogic::ellipse_scale; // radius-ish A, but a little shorter
  const double mrb = ellipsoid.mry_o() * vtkSlicerSRepCreatorLogic::ellipse_scale; // radius-ish B, but a little shorter

  const double deltaTheta = 2*vtkMath::Pi()/(numFoldPoints);
  const double stepSize = 1.0 / numStepsToCrest;

  // the +1 is for the point on the spine
  Eigen::MatrixXd reformed_points_x(numFoldPoints, numStepsToCrest+1);
  Eigen::MatrixXd reformed_points_y(numFoldPoints, numStepsToCrest+1);

  for (size_t i = 0; i < numFoldPoints; ++i) {
    // go around the entire ellipse radially, starting at pi radians
    const double theta = vtkMath::Pi() - deltaTheta*i;
    const double x = mra * cos(theta); // x of the final step on the boundary for this line
    const double y = mrb * sin(theta); // y of the final step on the boundary for this line
    //interesting properties of this computation of mx_
    // 1) For a perfect circle, the length of the spine is 0
    // 2) For a degenerate ellipse where mrb = 0 (ellipse flattens to a line), the length of the spine is mra
    const double mx_ = (mra * mra - mrb * mrb) * cos(theta) / mra; // this is the middle line (aka the spine)
    const double my_ = .0; //y is always zero because this is the middle line (aka the spine)

    // get distances in x,y directions between the first step on the spine and the last step
    // on the boundary
    const double dx_ = x - mx_;
    const double dy_ = y - my_;

    for (size_t j = 0; j < numStepsToCrest + 1; ++j) {
      const double tempX_ = mx_ + stepSize * j * dx_; // x of this step on the line
      const double tempY_ = my_ + stepSize * j * dy_; // y of this step on the line
      reformed_points_x(i, j) = tempX_;
      reformed_points_y(i, j) = tempY_;
    }
  }

  return std::make_pair(reformed_points_x, reformed_points_y);
}

//---------------------------------------------------------------------------
vtkSlicerSRepCreatorLogic::EigenSRep
vtkSlicerSRepCreatorLogic::GenerateEigenSRep(
  const EllipsoidParameters& ellipsoid,
  const size_t numFoldPoints,
  const size_t numStepsToCrest)
{
  const auto medial_sheet = vtkSlicerSRepCreatorLogic::GenerateMedialSkeletalSheet(ellipsoid, numFoldPoints, numStepsToCrest);
  const auto& reformed_points_x = medial_sheet.first;
  const auto& reformed_points_y = medial_sheet.second;

  // compute head points of spokes
  EigenSRep preTransformSRep(numFoldPoints, numStepsToCrest);

  const double mrx_o = ellipsoid.mrx_o();
  const double mry_o = ellipsoid.mry_o();
  const double rz = ellipsoid.radii(0);
  const double ry = ellipsoid.radii(1);
  const double rx = ellipsoid.radii(2);

  // the last row of the reformed points x/y is the fold
  for (int i = 0; i < numFoldPoints; ++i) {
    for (int j = 0; j < numStepsToCrest + 1; ++j) { //+1 for the spine point
      const double mx = reformed_points_x(i,j);
      const double my = reformed_points_y(i,j);

      const double sB = my * mrx_o;
      const double cB = mx * mry_o;
      const double l = sqrt(sB*sB + cB*cB);
      // sin(theta)
      const double sB_n = l < vtkSlicerSRepCreatorLogic::eps ? sB : sB / l;
      // cos(theta)
      const double cB_n = l < vtkSlicerSRepCreatorLogic::eps ? cB : cB / l;

      const double cA = l / (mrx_o * mry_o); // cos(phi)
      const double sA = sqrt(1 - cA*cA); // sin(phi)
      const double sx = rx * cA * cB_n - mx;
      const double sy = ry * cA * sB_n - my;
      const double sz = rz * sA;

      const double bx = sx + mx; // up/down spoke boundary point x
      const double by = sy + my; // up/down spoke boundary point y
      const double bz = sz; // up spoke boundary point z, negative of down spoke boundary point z

      const int id_resampled = i * (numStepsToCrest + 1) + j;
      preTransformSRep.skeletalPoints.row(id_resampled) << mx, my, 0.0;
      preTransformSRep.upSpokeBoundaryPoints.row(id_resampled) << bx, by, bz;
      preTransformSRep.downSpokeBoundaryPoints.row(id_resampled) << bx, by, -bz;

      if (j == numStepsToCrest) {
        // we are on the crest (aka fold)
        const double cx = rx * cB_n - mx;
        const double cy = ry * sB_n - my;
        const double cz = 0;
        const Eigen::Vector3d v {cx, cy, cz};
        const double v_n = v.norm();
        Eigen::Vector3d v2 {sx, sy, 0.0};
        v2.normalize(); // v2 is the unit vector pointing out to norm dir
        Eigen::Vector3d v3 = v_n * v2;

        const double bx = (v3(0) + mx);
        const double by = (v3(1) + my);

        // shift the skeletal side of the crest spoke off the interior skeleton toward
        // the boundary
        const double cmx = mx + (bx - mx) * vtkSlicerSRepCreatorLogic::crestShift;
        const double cmy = my + (by - my) * vtkSlicerSRepCreatorLogic::crestShift;
        
        preTransformSRep.crestSpokeBoundaryPoints.row(i) << bx, by, 0.0;
        preTransformSRep.crestSkeletalPoints.row(i) << cmx, cmy, 0.0;
      }
    }
  }

  //rotation and translation
  Eigen::MatrixXd transpose_srep = preTransformSRep.skeletalPoints.transpose(); // 3xn
  Eigen::Matrix3d srep_secondMoment = transpose_srep * preTransformSRep.skeletalPoints; // 3x3
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es_srep(srep_secondMoment);
  const auto rot_srep = es_srep.eigenvectors().transpose();

  auto rotation = ellipsoid.rotation;

  //TODO: Add rotate X/Y/Z stuff here?

  rotation = rotation * rot_srep;

  //transform points
  const auto transformMat = [&ellipsoid, &rotation](const Eigen::MatrixXd& mat) {
    auto transRotatedMat = mat * rotation.transpose();
    return transRotatedMat + ellipsoid.center.replicate(transRotatedMat.rows(), 1);
  };
  const auto transformSRep = [transformMat](const EigenSRep& srep) {
    EigenSRep transformed;
    transformed.skeletalPoints = transformMat(srep.skeletalPoints);
    transformed.upSpokeBoundaryPoints = transformMat(srep.upSpokeBoundaryPoints);
    transformed.downSpokeBoundaryPoints = transformMat(srep.downSpokeBoundaryPoints);
    transformed.crestSkeletalPoints = transformMat(srep.crestSkeletalPoints);
    transformed.crestSpokeBoundaryPoints = transformMat(srep.crestSpokeBoundaryPoints);
    transformed.numFoldPoints = srep.numFoldPoints;
    transformed.numStepsToCrest = srep.numStepsToCrest;
    return transformed;
  };

  return transformSRep(preTransformSRep);
}

//---------------------------------------------------------------------------
std::unique_ptr<srep::EllipticalSRep> vtkSlicerSRepCreatorLogic::ConvertEigenSRepToEllipticalSRep(
  const EigenSRep& eigenSRep)
{
  srep::EllipticalSRep::UnrolledEllipticalGrid grid;

  grid.reserve(eigenSRep.numFoldPoints);
  for (size_t i = 0; i < eigenSRep.numFoldPoints; ++i) {
    grid.resize(grid.size() + 1);
    grid.back().reserve(eigenSRep.numStepsToCrest + 1);
    for (size_t j = 0; j < eigenSRep.numStepsToCrest + 1; ++j) {
      const auto upDownSkeletonIndex = i * (eigenSRep.numStepsToCrest + 1) + j;
      const auto skeletalPoint = PointFromEigen(eigenSRep.skeletalPoints.row(upDownSkeletonIndex));
      const auto upBoundaryPoint = PointFromEigen(eigenSRep.upSpokeBoundaryPoints.row(upDownSkeletonIndex));
      const auto downBoundaryPoint = PointFromEigen(eigenSRep.downSpokeBoundaryPoints.row(upDownSkeletonIndex));
      srep::Spoke upSpoke = srep::Spoke(skeletalPoint, srep::Vector3d(skeletalPoint, upBoundaryPoint));
      srep::Spoke downSpoke = srep::Spoke(skeletalPoint, srep::Vector3d(skeletalPoint, downBoundaryPoint));

      if (j == eigenSRep.numStepsToCrest) {
        //crest
        const auto crestSkeletalPoint = PointFromEigen(eigenSRep.crestSkeletalPoints.row(i));
        const auto crestBoundaryPoint = PointFromEigen(eigenSRep.crestSpokeBoundaryPoints.row(i));
        srep::Spoke crestSpoke = srep::Spoke(crestSkeletalPoint, srep::Vector3d(crestSkeletalPoint, crestBoundaryPoint));
        grid.back().emplace_back(std::move(upSpoke), std::move(downSpoke), std::move(crestSpoke));
      } else {
        //interior, not crest
        grid.back().emplace_back(std::move(upSpoke), std::move(downSpoke));
      }
    }
  }
  return std::unique_ptr<srep::EllipticalSRep>(new srep::EllipticalSRep(std::move(grid)));
}

//---------------------------------------------------------------------------
std::unique_ptr<srep::EllipticalSRep> vtkSlicerSRepCreatorLogic::GenerateSRep(
  const EllipsoidParameters& ellipsoid,
  const size_t numFoldPoints,
  const size_t numStepsToCrest)
{
  auto eigenSRep = vtkSlicerSRepCreatorLogic::GenerateEigenSRep(ellipsoid, numFoldPoints, numStepsToCrest);
  return vtkSlicerSRepCreatorLogic::ConvertEigenSRepToEllipticalSRep(eigenSRep);
}

//---------------------------------------------------------------------------
void vtkSlicerSRepCreatorLogic::Reset() {
  this->ActualForwardIterations = 0;
  this->SRepNodeId.clear();
  this->ModelName.clear();
}

//---------------------------------------------------------------------------
vtkMRMLEllipticalSRepNode* vtkSlicerSRepCreatorLogic::RunForward(
  vtkMRMLModelNode* model,
  const size_t numFoldPoints,
  const size_t numStepsToCrest,
  const double dt,
  const double smoothAmount,
  const size_t maxIterations,
  const bool outputEllipsoidModel,
  const size_t outputEveryNumIterations)
{
  this->Reset();
  try {
    this->ModelName = model->GetName();
    const auto ellipsoidParameters = this->FlowSurfaceMeshToEllipsoid(model, dt, smoothAmount, maxIterations, outputEveryNumIterations);

    if (outputEllipsoidModel) {
      this->MakeEllipsoidModelNode(ellipsoidParameters, model->GetName() + std::string("-best-fit-ellipsoid"));
    }

    auto srepNode = this->MakeEllipticalSRepNode(
      this->GenerateSRep(ellipsoidParameters, numFoldPoints, numStepsToCrest),
      model->GetName() + std::string("-forward-only-srep"));
    if (srepNode) {
      this->SRepNodeId = srepNode->GetID();
      auto originalModelDisplayNode = model->GetDisplayNode();
      if (originalModelDisplayNode) {
        originalModelDisplayNode->SetVisibility(false);
      }
      return srepNode;
    }
    this->Reset();
    return nullptr;
  } catch (const std::exception& e) {
    vtkErrorMacro("Exception caught creating SRep: " << e.what());
    this->Reset();
    return nullptr;
  } catch (...) {
    vtkErrorMacro("Unknown exception caught creating SRep");
    this->Reset();
    return nullptr;
  }
}

//---------------------------------------------------------------------------
vtkMRMLEllipticalSRepNode* vtkSlicerSRepCreatorLogic::RunBackward(const size_t outputEveryNumIterations) {
  try {
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

    auto mrmlScene = this->GetMRMLScene();
    if (!mrmlScene) {
      vtkErrorMacro("vtkSlicerSRepCreatorLogic::RunBackward() cannot find mrmlScene");
      return nullptr;
    }
    auto srepNode = vtkMRMLEllipticalSRepNode::SafeDownCast(mrmlScene->GetNodeByID(this->SRepNodeId));
    if (!srepNode) {
      vtkErrorMacro("vtkSlicerSRepCreatorLogic::RunBackward() cannot find srepNode: " + this->SRepNodeId);
      return nullptr;
    }

    auto srep = srepNode->GetEllipticalSRep();
    if (!srep) {
      vtkErrorMacro("vtkSlicerSRepCreatorLogic::RunBackward() cannot find srep: " + this->SRepNodeId);
      return nullptr;
    }
    //copy the grid
    auto grid = srep->GetSkeleton();

    vtkNew<vtkPolyDataReader> reader1;
    vtkNew<vtkPolyDataReader> reader2;

    vtkPolyDataReader* sourceSurfaceReader = reader1;
    vtkPolyDataReader* targetSurfaceReader = reader2;

    sourceSurfaceReader->SetFileName(this->ForwardIterationFilename(this->ActualForwardIterations).c_str());
    sourceSurfaceReader->Update();

    std::cout << "Computing backflow transformation for " << this->ActualForwardIterations << " iterations" << std::endl;
    for (long iteration = this->ActualForwardIterations; iteration > 1; --iteration) {
      //swap source and target at bottom because target becomes source
      targetSurfaceReader->SetFileName(this->ForwardIterationFilename(iteration - 1).c_str());
      targetSurfaceReader->Update();

      vtkSmartPointer<vtkPolyData> polyData_source = sourceSurfaceReader->GetOutput();
      vtkSmartPointer<vtkPolyData> polyData_target = targetSurfaceReader->GetOutput();

      PointSetType::Pointer sourceLandMarks = PointSetType::New();
      PointSetType::Pointer targetLandMarks = PointSetType::New();
      PointSetType::PointsContainer::Pointer sourceLandMarkContainer
                  = sourceLandMarks->GetPoints();
      PointSetType::PointsContainer::Pointer targetLandMarkContainer
                  = targetLandMarks->GetPoints();

      // Read in the source points set
      for(unsigned int i = 0, s = 0; i < polyData_source->GetNumberOfPoints(); i += 10, ++s) {
          sourceLandMarkContainer->InsertElement(s, polyDataPointToPointType(*polyData_source, i));
      }

      // Read in the target points set
      for(unsigned int i = 0, s = 0; i < polyData_target->GetNumberOfPoints(); i += 10, ++s) {
          targetLandMarkContainer->InsertElement(s, polyDataPointToPointType(*polyData_target, i));
      }

      TransformType::Pointer tps = TransformType::New();
      tps->SetSourceLandmarks(sourceLandMarks);
      tps->SetTargetLandmarks(targetLandMarks);
      tps->ComputeWMatrix();

      ApplyTPSInPlace(grid, tps);

      if (outputEveryNumIterations != 0 && iteration % outputEveryNumIterations == 0) {
        // deep copy the grid
        std::unique_ptr<srep::EllipticalSRep> partiallyTransformedSRep(new srep::EllipticalSRep(grid));
        this->MakeEllipticalSRepNode(std::move(partiallyTransformedSRep), this->ModelName + "-backflow-srep-" + std::to_string(iteration));
      }

      std::swap(sourceSurfaceReader, targetSurfaceReader);
    }

    // can move the grid now since we aren't transforming it anymore
    std::unique_ptr<srep::EllipticalSRep> transformedSRep(new srep::EllipticalSRep(std::move(grid)));
    auto transformedSRepNode = this->MakeEllipticalSRepNode(std::move(transformedSRep), this->ModelName + "-srep");
    return transformedSRepNode;
  } catch (const std::exception& e) {
    vtkErrorMacro("Exception caught backflowing SRep: " << e.what());
    return nullptr;
  } catch (...) {
    vtkErrorMacro("Unknown exception caught backflowing SRep");
    return nullptr;
  }
}

//---------------------------------------------------------------------------
vtkMRMLEllipticalSRepNode* vtkSlicerSRepCreatorLogic::Run(
  vtkMRMLModelNode* model,
  const size_t numFoldPoints,
  const size_t numStepsToCrest,
  const double dt,
  const double smoothAmount,
  const size_t maxIterations,
  bool outputEllipsoidModel,
  size_t forwardOutputEveryNumIterations,
  size_t backwardOutputEveryNumIterations)
{
  auto ellipsoidSRep = this->RunForward(model, numFoldPoints, numStepsToCrest, dt,
      smoothAmount, maxIterations, outputEllipsoidModel, forwardOutputEveryNumIterations);

  if (ellipsoidSRep) {
    auto initialFitSRep = this->RunBackward(backwardOutputEveryNumIterations);
    if (!outputEllipsoidModel) {
      this->GetMRMLScene()->RemoveNode(ellipsoidSRep);
    }
    return initialFitSRep;
  }
  return nullptr;
}
