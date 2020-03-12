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

// SkeletalRepresentationInitializer Logic includes
#include "vtkSlicerSkeletalRepresentationInitializerLogic.h"

// MRML includes
#include <vtkMRMLScene.h>
#include <vtkMRMLModelNode.h>
#include <vtkMRMLModelDisplayNode.h>
#include <vtkMRMLDisplayNode.h>
#include <vtkMRMLMarkupsDisplayNode.h>
#include <vtkMRMLMarkupsFiducialNode.h>
#include <vtkMRMLMarkupsNode.h>
#include "vtkSlicerMarkupsLogic.h"

// VTK includes
#include <vtkCenterOfMass.h>
#include <vtkCurvatures.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkImplicitPolyDataDistance.h>
#include <vtkIntArray.h>
#include <vtkLine.h>
#include <vtkMassProperties.h>
#include <vtkMath.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkParametricEllipsoid.h>
#include <vtkParametricFunctionSource.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkQuad.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkVector.h>
#include <vtkWindowedSincPolyDataFilter.h>
#include <vtkAppendPolyData.h>
#include <vtkConeSource.h>
// Eigen includes
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

// STD includes
#include <cassert>
#include <iostream>

// vtk system tools
#include <vtksys/SystemTools.hxx>

#define MAX_FILE_NAME  256
//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkSlicerSkeletalRepresentationInitializerLogic);

//----------------------------------------------------------------------------
vtkSlicerSkeletalRepresentationInitializerLogic::vtkSlicerSkeletalRepresentationInitializerLogic()
{
}

//----------------------------------------------------------------------------
vtkSlicerSkeletalRepresentationInitializerLogic::~vtkSlicerSkeletalRepresentationInitializerLogic()
{
}

//----------------------------------------------------------------------------
void vtkSlicerSkeletalRepresentationInitializerLogic::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

//---------------------------------------------------------------------------
void vtkSlicerSkeletalRepresentationInitializerLogic::SetMRMLSceneInternal(vtkMRMLScene * newScene)
{
  vtkNew<vtkIntArray> events;
  events->InsertNextValue(vtkMRMLScene::NodeAddedEvent);
  events->InsertNextValue(vtkMRMLScene::NodeRemovedEvent);
  events->InsertNextValue(vtkMRMLScene::EndBatchProcessEvent);
  this->SetAndObserveMRMLSceneEventsInternal(newScene, events.GetPointer());
}

//-----------------------------------------------------------------------------
void vtkSlicerSkeletalRepresentationInitializerLogic::RegisterNodes()
{
  assert(this->GetMRMLScene() != nullptr);
}

//---------------------------------------------------------------------------
void vtkSlicerSkeletalRepresentationInitializerLogic::UpdateFromMRMLScene()
{
  assert(this->GetMRMLScene() != nullptr);
}

//---------------------------------------------------------------------------
void vtkSlicerSkeletalRepresentationInitializerLogic
::OnMRMLSceneNodeAdded(vtkMRMLNode* vtkNotUsed(node))
{
}

//---------------------------------------------------------------------------
void vtkSlicerSkeletalRepresentationInitializerLogic
::OnMRMLSceneNodeRemoved(vtkMRMLNode* vtkNotUsed(node))
{
}

// flow surface in one step
// basic idea: When the user select a mesh file, make a copy of vtk file in the application path.
// In each step of flow, read in that copy, flow it and save it the same place with same name.
// TODO: cleanup the hard disk when the module exits
int vtkSlicerSkeletalRepresentationInitializerLogic::FlowSurfaceOneStep(const std::string &filename, double dt, double smooth_amount)
{
    vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();
    vtkSmartPointer<vtkPolyData> mesh = reader->GetOutput();
    if(mesh == nullptr)
    {
        vtkErrorMacro("No mesh has read in this module. Please select input mesh file first.");
        return -1;
    }
    vtkSmartPointer<vtkMassProperties> mass_filter =
        vtkSmartPointer<vtkMassProperties>::New();
    mass_filter->SetInputData(mesh);
    mass_filter->Update();
    double original_volume = mass_filter->GetVolume();
//    std::cout << "Original Volume: " << original_volume << std::endl;
    vtkSmartPointer<vtkWindowedSincPolyDataFilter> smooth_filter =
        vtkSmartPointer<vtkWindowedSincPolyDataFilter>::New();
    smooth_filter->SetPassBand(smooth_amount);
    smooth_filter->NonManifoldSmoothingOn();
    smooth_filter->NormalizeCoordinatesOn();
    smooth_filter->SetNumberOfIterations(20);
    smooth_filter->FeatureEdgeSmoothingOff();
    smooth_filter->BoundarySmoothingOff();
    smooth_filter->SetInputData(mesh);
    smooth_filter->Update();
    if(smooth_amount > 0) {
        mesh = smooth_filter->GetOutput();
    }

//    normal filter
    vtkSmartPointer<vtkPolyDataNormals> normal_filter =
        vtkSmartPointer<vtkPolyDataNormals>::New();
    normal_filter->SplittingOff();
    normal_filter->ComputeCellNormalsOff();
    normal_filter->ComputePointNormalsOn();
    normal_filter->SetInputData(mesh);
    normal_filter->Update();
    // curvature filter
    vtkSmartPointer<vtkCurvatures> curvature_filter =
        vtkSmartPointer<vtkCurvatures>::New();
    curvature_filter->SetCurvatureTypeToMean();
    curvature_filter->SetInputData(mesh);
    curvature_filter->Update();

    // compute the center of surface mesh
    vtkSmartPointer<vtkCenterOfMass> centerMassFilter =
        vtkSmartPointer<vtkCenterOfMass>::New();
    centerMassFilter->SetInputData(mesh);
    centerMassFilter->SetUseScalarsAsWeights(false);
    centerMassFilter->Update();
    double center[3];
    centerMassFilter->GetCenter(center);

    vtkSmartPointer<vtkDoubleArray> H =
        vtkDoubleArray::SafeDownCast(curvature_filter->GetOutput()->GetPointData()->GetArray("Mean_Curvature"));

    curvature_filter->SetCurvatureTypeToGaussian();
    curvature_filter->Update();
    vtkSmartPointer<vtkDoubleArray> K =
            vtkDoubleArray::SafeDownCast(curvature_filter->GetOutput()->GetPointData()->GetArray("Gauss_Curvature"));
    if(H == nullptr) {
        vtkErrorMacro("error in getting mean curvature");
        return -1;
    }
    vtkDataArray* N = normal_filter->GetOutput()->GetPointData()->GetNormals();
    if(N == nullptr) {
        vtkErrorMacro("error in getting normals");
        return -1;
    }

    // perform the flow
    vtkSmartPointer<vtkPoints> points = mesh->GetPoints();

    // Gaussian map
    vtkSmartPointer<vtkPolyData> spherePolys =
            vtkSmartPointer<vtkPolyData>::New();

    vtkSmartPointer<vtkPoints> spherePts = vtkSmartPointer<vtkPoints>::New();

    // concave pts
    vtkSmartPointer<vtkPolyData> concavePolyData = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints> concavePts = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> concavePolys = vtkSmartPointer<vtkCellArray>::New();

    // hyperbolic pts
    vtkSmartPointer<vtkPolyData> hyperPolyData = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints> hyperPts = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> hyperPolys = vtkSmartPointer<vtkCellArray>::New();
    for(int i = 0; i < points->GetNumberOfPoints(); ++i) {
        double p[3];
        points->GetPoint(i, p);
        double curr_N[3];
        N->GetTuple(i, curr_N);
        double curr_H = H->GetValue(i);
        double curr_K = K->GetValue(i);
        if(curr_K < 0)
        {
            vtkIdType newId = hyperPts->InsertNextPoint(p);
            GetNeighborCells(mesh, i, static_cast<int>(newId), hyperPolys, hyperPts);
        }
        else {
            // see the explanation on vtkCurvatures signs
            if(curr_H < 0)
            {
                vtkIdType newId = concavePts->InsertNextPoint(p);
                GetNeighborCells(mesh, i, static_cast<int>(newId), concavePolys, concavePts);
            }
        }

        double ptSphere[3];
        ptSphere[0] = curr_N[0] * std::pow( original_volume , 1.0 / 3.0 ) / vtkMath::Pi() + center[0];
        ptSphere[1] = curr_N[1] * std::pow( original_volume , 1.0 / 3.0 ) / vtkMath::Pi() + center[1];
        ptSphere[2] = curr_N[2] * std::pow( original_volume , 1.0 / 3.0 ) / vtkMath::Pi() + center[2];
        spherePts->InsertNextPoint(ptSphere);

        for(int idx = 0; idx < 3; ++idx) {
            p[idx] -= dt * curr_H * curr_N[idx];
        }
        points->SetPoint(i, p);
    }
    points->Modified();
    spherePts->Modified();
    hyperPts->Modified();
    hyperPolys->Modified();

    concavePts->Modified();
    concavePolys->Modified();
    concavePolyData->SetPoints(concavePts);
    concavePolyData->SetPolys(concavePolys);
    concavePolyData->Modified();

    spherePolys->SetPoints(spherePts);
    spherePolys->SetPolys(mesh->GetPolys());
    spherePolys->Modified();

    hyperPolyData->SetPoints(hyperPts);
    hyperPolyData->SetPolys(hyperPolys);
    hyperPolyData->Modified();

    const std::string hyperbolicRegionName("hyperbolic_points");
    const std::string concaveRegionName("concave_points");
    const std::string sphereName("gauss_sphere_map");
    const std::string modelName("curvature_flow_result");

    // firstly get other intermediate result invisible
    HideNodesByNameByClass(modelName.c_str(),"vtkMRMLModelNode");
    HideNodesByNameByClass(sphereName.c_str(),"vtkMRMLModelNode");
    HideNodesByNameByClass(hyperbolicRegionName.c_str(),"vtkMRMLModelNode");
    HideNodesByNameByClass(concaveRegionName.c_str(),"vtkMRMLModelNode");
    HideNodesByNameByClass("best_fitting_ellipsoid_polydata", "vtkMRMLModelNode");

    // then add this new intermediate result
    //AddModelNodeToScene(spherePolys, sphereName.c_str(), true, 1, 0,0);
    AddModelNodeToScene(mesh, modelName.c_str(), true);
    AddModelNodeToScene(hyperPolyData, hyperbolicRegionName.c_str(), true, 1, 0, 1);
    AddModelNodeToScene(concavePolyData, concaveRegionName.c_str(), true, 0,0,1);
    vtkSmartPointer<vtkPolyDataWriter> writer =
        vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetInputData(mesh);
    writer->SetFileName(filename.c_str());
    writer->Update();
    return 0;
}
void vtkSlicerSkeletalRepresentationInitializerLogic::SetInputFileName(const std::string &filename)
{
    vtkSmartPointer<vtkPolyDataReader> reader =
        vtkSmartPointer<vtkPolyDataReader>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();

    vtkSmartPointer<vtkPolyData> mesh;
    mesh = reader->GetOutput();
    // output the original mesh
    const std::string modelName("original");
    AddModelNodeToScene(mesh, modelName.c_str(), true, 0.88, 0.88, 0.88);

    // save
    const std::string tempDir(this->GetApplicationLogic()->GetTemporaryPath());
    std::string tempFileName;
    tempFileName = tempDir + "/initial_surface.vtk";

    vtkSmartPointer<vtkPolyDataWriter> writer =
        vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetInputData(mesh);
    writer->SetFileName(tempFileName.c_str());
    writer->Update();
}

// flow surface to the end: either it's ellipsoidal enough or reach max_iter
int vtkSlicerSkeletalRepresentationInitializerLogic::FlowSurfaceMesh(const std::string &filename, double dt, double smooth_amount, int max_iter, int freq_output)
{
    vtkSmartPointer<vtkPolyDataReader> reader =
        vtkSmartPointer<vtkPolyDataReader>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();

    vtkSmartPointer<vtkPolyData> mesh =
        vtkSmartPointer<vtkPolyData>::New();
    mesh = reader->GetOutput();

    vtkSmartPointer<vtkMassProperties> mass_filter =
        vtkSmartPointer<vtkMassProperties>::New();
    mass_filter->SetInputData(mesh);
    mass_filter->Update();
    double original_volume = mass_filter->GetVolume();
    // default parameters
    // double dt = 0.001;
    // double smooth_amount = 0.03;
    // int max_iter = 500;
    int iter = 0;
    double tolerance = 0.05;
    double q = 1.0;

    // create folder if not exist
    const std::string tempFolder(this->GetApplicationLogic()->GetTemporaryPath());
    std::string forwardFolder;
    forwardFolder = tempFolder + "/forward";
    std::cout << "forward folder" << forwardFolder << std::endl;
    if (!vtksys::SystemTools::FileExists(forwardFolder, false))
    {
      if (!vtksys::SystemTools::MakeDirectory(forwardFolder))
      {
              std::cout << "Failed to create folder : " << forwardFolder << std::endl;

      }
    }
    while(q > tolerance && iter < max_iter) {
        // smooth filter
        vtkSmartPointer<vtkWindowedSincPolyDataFilter> smooth_filter =
            vtkSmartPointer<vtkWindowedSincPolyDataFilter>::New();
        smooth_filter->SetPassBand(smooth_amount);
        smooth_filter->NonManifoldSmoothingOn();
        smooth_filter->NormalizeCoordinatesOn();
        smooth_filter->SetNumberOfIterations(20);
        smooth_filter->FeatureEdgeSmoothingOff();
        smooth_filter->BoundarySmoothingOff();
        smooth_filter->SetInputData(mesh);
        smooth_filter->Update();
        if(smooth_amount > 0) {
            mesh = smooth_filter->GetOutput();
        }

        // normal filter
        vtkSmartPointer<vtkPolyDataNormals> normal_filter =
            vtkSmartPointer<vtkPolyDataNormals>::New();
        normal_filter->SplittingOff();
        normal_filter->ComputeCellNormalsOff();
        normal_filter->ComputePointNormalsOn();
        normal_filter->SetInputData(mesh);
        normal_filter->Update();
        // curvature filter
        vtkSmartPointer<vtkCurvatures> curvature_filter =
            vtkSmartPointer<vtkCurvatures>::New();
        curvature_filter->SetCurvatureTypeToMean();
        curvature_filter->SetInputData(mesh);
        curvature_filter->Update();


        vtkSmartPointer<vtkDoubleArray> H =
            vtkDoubleArray::SafeDownCast(curvature_filter->GetOutput()->GetPointData()->GetArray("Mean_Curvature"));
        if(H == nullptr) {
            std::cerr << "error in getting mean curvature" << std::endl;
            return EXIT_FAILURE;
        }
        vtkDataArray* N = normal_filter->GetOutput()->GetPointData()->GetNormals();
        if(N == nullptr) {
            std::cerr << "error in getting normals" << std::endl;
            return EXIT_FAILURE;
        }

        // perform the flow
        vtkSmartPointer<vtkPoints> points = mesh->GetPoints();
        for(int i = 0; i < points->GetNumberOfPoints(); ++i) {
            double p[3];
            points->GetPoint(i, p);
            double curr_N[3];
            N->GetTuple(i, curr_N);
            double curr_H = H->GetValue(i);
            for(int idx = 0; idx < 3; ++idx) {
                p[idx] -= dt * curr_H * curr_N[idx];
            }
            points->SetPoint(i, p);
        }
        points->Modified();
        mass_filter->SetInputData(mesh);
        mass_filter->Update();
        double curr_volume = mass_filter->GetVolume();
        for(int i = 0; i < points->GetNumberOfPoints(); ++i) {
            double p[3];
            points->GetPoint(i, p);
            for(int j = 0; j < 3; ++j) {
                p[j] *= std::pow( original_volume / curr_volume , 1.0 / 3.0 );
            }
//            points->SetPoint(i, p);
        }
        points->Modified();
        // save the result for the purpose of backward flow
        std::string fileName;
        fileName = forwardFolder + "/" + std::to_string(iter+1) + ".vtk";
        vtkSmartPointer<vtkPolyDataWriter> writer =
            vtkSmartPointer<vtkPolyDataWriter>::New();
        writer->SetInputData(mesh);
        writer->SetFileName(fileName.c_str());
        writer->Update();

        if((iter +1) % freq_output == 0)
        {
            const std::string modelName = "output" + std::to_string(iter+1);
            AddModelNodeToScene(mesh, modelName.c_str(), false);
            vtkSmartPointer<vtkCenterOfMass> centerMassFilter =
                vtkSmartPointer<vtkCenterOfMass>::New();
            centerMassFilter->SetInputData(mesh);
            centerMassFilter->SetUseScalarsAsWeights(false);
            centerMassFilter->Update();
            double center[3];
            centerMassFilter->GetCenter(center);
        }
        q -= 0.0001;
        iter++;
    }
    forwardCount = iter;
    double rx, ry, rz;
    ShowFittingEllipsoid(mesh, rx, ry, rz);

    GenerateSrepForEllipsoid(mesh, mRows, mCols, forwardCount);
    return 1;
}

void vtkSlicerSkeletalRepresentationInitializerLogic::AddModelNodeToScene(vtkPolyData* mesh, const char* modelName, bool isModelVisible, double r, double g, double b)
{
    std::cout << "AddModelNodeToScene: parameters:" << modelName << std::endl;
    vtkMRMLScene *scene = this->GetMRMLScene();
    if(!scene)
    {
        vtkErrorMacro(" Invalid scene");
        return;
    }

    // model node
    vtkSmartPointer<vtkMRMLModelNode> modelNode;
    modelNode = vtkSmartPointer<vtkMRMLModelNode>::New();
    modelNode->SetScene(scene);

    modelNode->SetName(modelName);
    modelNode->SetAndObservePolyData(mesh);

    // display node
    vtkSmartPointer<vtkMRMLModelDisplayNode> displayModelNode;

    displayModelNode = vtkSmartPointer<vtkMRMLModelDisplayNode>::New();
    if(displayModelNode == nullptr)
    {
        vtkErrorMacro("displayModelNode is NULL");
        return;
    }
    displayModelNode->SetColor(r, g, b);
    displayModelNode->SetScene(scene);
    displayModelNode->SetLineWidth(2.0);
    displayModelNode->SetBackfaceCulling(0);
    displayModelNode->SetRepresentation(1);
    if(isModelVisible)
    {
        // make the 1st mesh after flow visible
        displayModelNode->SetVisibility(1);
    }
    else
    {
        displayModelNode->SetVisibility(0);
    }

    scene->AddNode(displayModelNode);
    modelNode->AddAndObserveDisplayNodeID(displayModelNode->GetID());

    scene->AddNode(modelNode);

}
void vtkSlicerSkeletalRepresentationInitializerLogic::ShowFittingEllipsoid(vtkPolyData* mesh, double &rx, double &ry, double &rz)
{
    vtkSmartPointer<vtkPoints> points = mesh->GetPoints();
    Eigen::MatrixXd point_matrix(points->GetNumberOfPoints(), 3);
    for(int i = 0; i < points->GetNumberOfPoints(); ++i)
    {
        double p[3];
        points->GetPoint(i, p);
        point_matrix.row(i) << p[0], p[1], p[2];
    }
    // compute best fitting ellipsoid
    // For now assume that the surface is centered and rotationally aligned
    // 1. compute the second moment after centering the data points
    Eigen::MatrixXd center = point_matrix.colwise().mean();
    Eigen::MatrixXd centered_point_mat = point_matrix - center.replicate(point_matrix.rows(), 1);
    Eigen::MatrixXd point_matrix_transposed = centered_point_mat.transpose();
    Eigen::Matrix3d second_moment = point_matrix_transposed * centered_point_mat;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(second_moment);
    Eigen::VectorXd radii = es.eigenvalues();
    radii(0) = sqrt(radii(0));
    radii(1) = sqrt(radii(1));
    radii(2) = sqrt(radii(2));

    double ellipsoid_volume = 4 / 3.0 * vtkMath::Pi() * radii(0) * radii(1) * radii(2);
    vtkSmartPointer<vtkMassProperties> mass =
        vtkSmartPointer<vtkMassProperties>::New();
    mass->SetInputData(mesh);
    mass->Update();

    double volume_factor = pow(mass->GetVolume() / ellipsoid_volume, 1.0 / 3.0);
    radii(0) *= volume_factor;
    radii(1) *= volume_factor;
    radii(2) *= volume_factor;
    // obtain the best fitting ellipsoid from the second moment matrix
    vtkSmartPointer<vtkParametricEllipsoid> ellipsoid =
        vtkSmartPointer<vtkParametricEllipsoid>::New();
    ellipsoid->SetXRadius(radii(0));
    ellipsoid->SetYRadius(radii(1));
    ellipsoid->SetZRadius(radii(2));

    vtkSmartPointer<vtkParametricFunctionSource> parametric_function =
        vtkSmartPointer<vtkParametricFunctionSource>::New();
    parametric_function->SetParametricFunction(ellipsoid);
    parametric_function->SetUResolution(30);
    parametric_function->SetVResolution(30);
    parametric_function->Update();
    vtkSmartPointer<vtkPolyData> ellipsoid_polydata = parametric_function->GetOutput();

    using namespace Eigen;
    // Get ellipsoid points into the matrix
    MatrixXd ellipsoid_points_matrix(ellipsoid_polydata->GetNumberOfPoints(), 3);
    for(int i = 0; i < ellipsoid_polydata->GetNumberOfPoints(); ++i) {
        double p[3];
        ellipsoid_polydata->GetPoint(i,p);
        ellipsoid_points_matrix(i,0) = p[0];
        ellipsoid_points_matrix(i,1) = p[1];
        ellipsoid_points_matrix(i,2) = p[2];
    }
    MatrixXd rotation;
    rotation = es.eigenvectors(); // 3 by 3 rotation matrix

    // rotate the points
    MatrixXd rotated_ellipsoid_points = rotation * (ellipsoid_points_matrix.transpose());
    rotated_ellipsoid_points.transposeInPlace(); // n x 3
    // translate the points
    MatrixXd translated_points = rotated_ellipsoid_points + center.replicate(rotated_ellipsoid_points.rows(),1);

    // convert eigen matrix to vtk polydata
    vtkSmartPointer<vtkPolyData> best_fitting_ellipsoid_polydata =
        vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints> best_fitting_ellipsoid_points =
        vtkSmartPointer<vtkPoints>::New();
    for(int i = 0; i < translated_points.rows(); ++i) {
        double p[3] = {translated_points(i,0), translated_points(i,1), translated_points(i,2)};
        best_fitting_ellipsoid_points->InsertNextPoint(p);
    }
    best_fitting_ellipsoid_polydata->SetPoints(best_fitting_ellipsoid_points);
    best_fitting_ellipsoid_polydata->SetPolys(ellipsoid_polydata->GetPolys());
    best_fitting_ellipsoid_polydata->Modified();

    AddModelNodeToScene(best_fitting_ellipsoid_polydata, "best_fitting_ellipsoid", true, 1, 1, 0);
    rx = radii(2); ry = radii(1); rz = radii(0);

    // output to file
    const std::string tempFolder(this->GetApplicationLogic()->GetTemporaryPath());

    const std::string ellipsoidFile = tempFolder + "/ellipsoid.vtk";
    vtkSmartPointer<vtkPolyDataWriter> writer =
        vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetInputData(best_fitting_ellipsoid_polydata);
    writer->SetFileName(ellipsoidFile.c_str());
    writer->Update();

}

const double ELLIPSE_SCALE = 0.9;
const double EPS = 1e-6;
void vtkSlicerSkeletalRepresentationInitializerLogic::GenerateSrepForEllipsoid(vtkPolyData *mesh,
                                                                              int nRows, int nCols, int totalNum,
                                                                               bool rotateX, bool rotateY,
                                                                               bool rotateZ)
{
    // create folder if not exist
    const std::string tempFolder(this->GetApplicationLogic()->GetTemporaryPath());
    const std::string modelFolder = tempFolder + "/model";
    std::cout << "s-reps folder" << modelFolder << std::endl;
    if (!vtksys::SystemTools::FileExists(modelFolder, false))
    {
      if (!vtksys::SystemTools::MakeDirectory(modelFolder))
      {
        std::cout << "Failed to create folder : " << modelFolder << std::endl;
      }
    }

    // copy surface of ellipsoid to model folder
    const std::string ellSurfaceFile = tempFolder + "/ellipsoid.vtk";
    const std::string newEllSurfaceFile = tempFolder + "/forward/" + std::to_string(totalNum + 1) + ".vtk";
    vtksys::SystemTools::CopyAFile(ellSurfaceFile, newEllSurfaceFile, true);

    using namespace Eigen;
    // the number of rows should be odd number
    double shift = 0.1; // shift fold curve off the inner spokes

    // 1. derive the best fitting ellipsoid from the deformed mesh
    vtkSmartPointer<vtkPoints> points = mesh->GetPoints();
    MatrixXd point_matrix(points->GetNumberOfPoints(), 3);
    for(int i = 0; i < points->GetNumberOfPoints(); ++i)
    {
        double p[3];
        points->GetPoint(i, p);
        point_matrix.row(i) << p[0], p[1], p[2];
    }
    MatrixXd center = point_matrix.colwise().mean();
    MatrixXd centered_point_mat = point_matrix - center.replicate(point_matrix.rows(), 1);
    MatrixXd point_matrix_transposed = centered_point_mat.transpose();
    Matrix3d second_moment = point_matrix_transposed * centered_point_mat;
    SelfAdjointEigenSolver<Eigen::MatrixXd> es_obj(second_moment);
    VectorXd radii = es_obj.eigenvalues();

    // notations consistent with wenqi's slides
    double rz = sqrt(radii(0));
    double ry = sqrt(radii(1));
    double rx = sqrt(radii(2));
    double ellipsoid_volume = 4 / 3.0 * vtkMath::Pi() * rx * ry * rz;
    vtkSmartPointer<vtkMassProperties> mass =
        vtkSmartPointer<vtkMassProperties>::New();
    mass->SetInputData(mesh);
    mass->Update();

    double volume_factor = pow(mass->GetVolume() / ellipsoid_volume, 1.0 / 3.0);
    rz *= volume_factor;
    ry *= volume_factor;
    rx *= volume_factor;

    double mrx_o = (rx*rx-rz*rz)/rx; // m_a
    double mry_o = (ry*ry-rz*rz)/ry; // m_b
    double mrb = mry_o * ELLIPSE_SCALE;
    double mra = mrx_o * ELLIPSE_SCALE;

    // 2. compute the skeletal points
    int nCrestPoints = nRows*2 + (nCols-2)*2;
    double deltaTheta = 2*vtkMath::Pi()/nCrestPoints;
    MatrixXd skeletal_points_x(nRows, nCols);
    MatrixXd skeletal_points_y(nRows, nCols);
    //MatrixXd skeletal_points_z(nRows, nCols);
    int r = 0, c = 0;
    for(int i = 0; i < nCrestPoints; ++i)
    {
        double theta = vtkMath::Pi() - deltaTheta * floor(nRows/2) - deltaTheta*i;
        double x = mra * cos(theta);
        double y = mrb * sin(theta);

        // these crest points have no inward points (side or corner of the s-rep)
        skeletal_points_x(r, c) = x;
        skeletal_points_y(r, c) = y;
        //skeletal_points_z(r, c) = z;
        if(i < nCols - 1)
        {
            // top row of crest points
            c += 1;
        }
        else if(i < nCols - 1 + nRows - 1)
        {
            // right side col of crest points ( if the top-left point is the origin)
            r = r + 1;
        }
        else if(i < nCols - 1 + nRows - 1 + nCols - 1)
        {
            // bottom row of crest points
            c = c - 1;
        }
        else
        {
            // left side col of crest points
            r = r - 1;
        }
        if((i < nCols - 1 && i > 0) || (i > nCols + nRows - 2 && i < 2*nCols + nRows - 3))
        {
            // compute skeletal points inward
            double mx_ = (mra * mra - mrb * mrb) * cos(theta) / mra; // this is the middle line
            double my_ = .0;
            double dx_ = x - mx_;
            double dy_ = y - my_;
            int numSteps = static_cast<int>(floor(nRows/2)); // steps from crest point to the skeletal point
            double stepSize = 1.0 / double(numSteps); // step size on the half side of srep
            for(int j = 0; j <= numSteps; ++j)
            {
                double tempX_ = mx_ + stepSize * j * dx_;
                double tempY_ = my_ + stepSize * j * dy_;
                if(i < nCols - 1)
                {
                    // step from medial to top at current iteration on the top line
                    int currR = numSteps - j;
                    skeletal_points_x(currR, c-1) = tempX_;
                    skeletal_points_y(currR, c-1) = tempY_;
                }
                else
                {
                    int currR = j + numSteps;
                    skeletal_points_x(currR, c+1) = tempX_;
                    skeletal_points_y(currR, c+1) = tempY_;
                }

            }

        }
    }

    // 3. compute the head points of spokes
    MatrixXd skeletal_points(nRows*nCols, 3);
    MatrixXd bdry_points_up(nRows*nCols, 3);
    MatrixXd bdry_points_down(nRows*nCols, 3);
    MatrixXd bdry_points_crest(nCrestPoints, 3);
    MatrixXd skeletal_points_crest(nCrestPoints, 3);
    int id_pt = 0; int id_crest = 0;
    MatrixXd shift_dir(nCrestPoints, 3); // shift direction for every crest point
    for(int i = 0; i < nRows; ++i)
    {
        for(int j = 0; j < nCols; ++j)
        {
            double mx = skeletal_points_x(i,j);
            double my = skeletal_points_y(i,j);
            double sB = my * mrx_o;
            double cB = mx * mry_o;
            double l = sqrt(sB*sB + cB*cB);
            double sB_n, cB_n; // sin(theta), cos(theta)
            if(l < EPS)
            {
                sB_n = sB;
                cB_n = cB;
            }
            else
            {
                sB_n = sB / l;
                cB_n = cB / l;
            }
            double cA = l / (mrx_o * mry_o); // cos(phi)
            double sA = sqrt(1 - cA*cA); // sin(phi)
            double sx = rx * cA * cB_n - mx;
            double sy = ry * cA * sB_n - my;
            double sz = rz * sA;

            double bx = (sx + mx);
            double by = (sy + my);
            double bz = (sz);

            skeletal_points.row(id_pt) << mx, my, 0.0;
            bdry_points_up.row(id_pt) << bx, by, bz;
            bdry_points_down.row(id_pt) << bx, by, -bz;
            id_pt++;
            // fold curve
            if(i == 0 || i == nRows - 1 || j == 0 || j == nCols - 1)
            {
                double cx = rx * cB_n - mx;
                double cy = ry * sB_n - my;
                double cz = 0;
                Vector3d v, v2, v3;
                v << cx, cy, cz;
                v2 << sx, sy, 0.0;
                double v_n = v.norm();
                v2.normalize(); // v2 is the unit vector pointing out to norm dir
                v3 = v_n * v2;
                double bx = (v3(0) + mx);
                double by = (v3(1) + my);
                double bz = v3(2);
                bdry_points_crest.row(id_crest) << bx, by, bz;
                skeletal_points_crest.row(id_crest) << mx, my, 0.0;
                shift_dir.row(id_crest) << v2(0), v2(1), v2(2);
                id_crest++;
            }
        }
    }

    // 4. transform the s-rep
    MatrixXd transpose_srep = skeletal_points.transpose(); // 3xn
    Matrix3d srep_secondMoment = transpose_srep * skeletal_points; // 3x3
    SelfAdjointEigenSolver<Eigen::MatrixXd> es_srep(srep_secondMoment);

    Matrix3d rotation;
    std::string arrowName = "Before rotation orientations";

    rotation = es_obj.eigenvectors(); // 3 by 3 rotation relative to deformed object

    if(rotateX) {
        rotation(0,1) *= -1;
        rotation(1,1) *= -1;
        rotation(2,1) *= -1;
        arrowName = "After rotation orientations";
    }
    if(rotateY) {
        rotation(0,2) *= -1;
        rotation(1,2) *= -1;
        rotation(2,2) *= -1;
        arrowName = "After rotation orientations";
    }
    if(rotateZ) {
        rotation(0,0) *= -1;
        rotation(1,0) *= -1;
        rotation(2,0) *= -1;
        arrowName = "After rotation orientations";
    }
    Matrix3d rot_srep;
    rot_srep = es_srep.eigenvectors().transpose();
    rotation = rotation * rot_srep;

    // all skeletal points
    MatrixXd trans_srep = (rotation * transpose_srep).transpose();
    MatrixXd transformed_skeletal_points = trans_srep+
            center.replicate(trans_srep.rows(), 1);

    // up spoke head point on the bdry
    MatrixXd transpose_up_pdm = bdry_points_up.transpose();
    MatrixXd trans_up_pdm = (rotation * transpose_up_pdm).transpose();
    MatrixXd transformed_up_pdm =  trans_up_pdm +
            center.replicate(trans_up_pdm.rows(), 1);

    // down spoke head point on the bdry
    MatrixXd transpose_down_pdm = bdry_points_down.transpose();
    MatrixXd trans_down_pdm = (rotation * transpose_down_pdm).transpose();
    MatrixXd transformed_down_pdm = trans_down_pdm +
            center.replicate(trans_down_pdm.rows(), 1);

    // crest head point on the bdry
    MatrixXd transpose_crest_pdm = bdry_points_crest.transpose();
    MatrixXd trans_crest_pdm = (rotation * transpose_crest_pdm).transpose();
    MatrixXd transformed_crest_pdm = trans_crest_pdm + center.replicate(trans_crest_pdm.rows(), 1);

    // crest base point on the skeletal sheet
    MatrixXd transpose_crest_base = skeletal_points_crest.transpose();
    MatrixXd trans_crest_base = (rotation * transpose_crest_base).transpose();
    MatrixXd transformed_crest_base = trans_crest_base + center.replicate(trans_crest_base.rows(), 1);

    // 5. transfer points to polydata
    // srep_poly is supposed to form a mesh grid connecting skeletal points
    vtkSmartPointer<vtkPolyData>  srep_poly       = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints>    skeletal_sheet  = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> skeletal_mesh   = vtkSmartPointer<vtkCellArray>::New();

    vtkSmartPointer<vtkPolyData>  upSpokes_poly      = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints>    upSpokes_pts       = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> upSpokes_lines     = vtkSmartPointer<vtkCellArray>::New();

    vtkSmartPointer<vtkPolyData>  downSpokes_poly      = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints>    downSpokes_pts       = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> downSpokes_lines     = vtkSmartPointer<vtkCellArray>::New();

    // TODO:crest spokes should be a little off the inner spokes
    vtkSmartPointer<vtkPolyData>  crestSpokes_poly      = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints>    crestSpokes_pts       = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> crestSpokes_lines     = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPolyData> foldCurve_poly         = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints>    foldCurve_pts         = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> fold_curve            = vtkSmartPointer<vtkCellArray>::New();

    skeletal_sheet->SetDataTypeToDouble();
    upSpokes_pts->SetDataTypeToDouble();
    downSpokes_pts->SetDataTypeToDouble();
    crestSpokes_pts->SetDataTypeToDouble();
    foldCurve_pts->SetDataTypeToDouble();

    vtkSmartPointer<vtkDoubleArray> upSpokeLengths = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> downSpokeLengths = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> crestSpokeLengths = vtkSmartPointer<vtkDoubleArray>::New();
    upSpokeLengths->SetNumberOfComponents(1);
    downSpokeLengths->SetNumberOfComponents(1);
    crestSpokeLengths->SetNumberOfComponents(1);

    upSpokeLengths->SetName("spokeLength");
    downSpokeLengths->SetName("spokeLength");
    crestSpokeLengths->SetName("spokeLength");

    vtkSmartPointer<vtkDoubleArray> upSpokeDirs = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> downSpokeDirs = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> crestSpokeDirs = vtkSmartPointer<vtkDoubleArray>::New();

    upSpokeDirs->SetNumberOfComponents(3);
    downSpokeDirs->SetNumberOfComponents(3);
    crestSpokeDirs->SetNumberOfComponents(3);

    upSpokeDirs->SetName("spokeDirection");
    downSpokeDirs->SetName("spokeDirection");
    crestSpokeDirs->SetName("spokeDirection");


    for(int i = 0; i < nRows * nCols; ++i)
    {
        // skeletal points
        double mx = transformed_skeletal_points(i,0);
        double my = transformed_skeletal_points(i,1);
        double mz = transformed_skeletal_points(i,2);
        int id0 = static_cast<int>(upSpokes_pts->InsertNextPoint(mx,my, mz));

        double bx_up = transformed_up_pdm(i, 0);
        double by_up = transformed_up_pdm(i, 1);
        double bz_up = transformed_up_pdm(i, 2);
        int id1 = static_cast<int>(upSpokes_pts->InsertNextPoint(bx_up, by_up, bz_up));
        if(i == nRows * nCols / 2) {
            vtkSmartPointer<vtkConeSource> cone = vtkSmartPointer<vtkConeSource>::New();
            cone->SetCenter(bx_up, by_up, bz_up);
            cone->SetDirection(bx_up-mx, by_up-my, bz_up-mz);
            cone->Update();
            AddModelNodeToScene(cone->GetOutput(), arrowName.c_str(), true, 0, 0, 1);
        }

        // spoke length and dir
        vtkVector3d upSpoke(bx_up-mx, by_up-my, bz_up-mz);
        double upSpokeLength = upSpoke.Normalize();
        upSpokeLengths->InsertNextTuple1(upSpokeLength);
        upSpokeDirs->InsertNextTuple3(upSpoke.GetX(), upSpoke.GetY(), upSpoke.GetZ());

        // form up spokes
        vtkSmartPointer<vtkLine> up_arrow = vtkSmartPointer<vtkLine>::New();
        up_arrow->GetPointIds()->SetId(0, id0);
        up_arrow->GetPointIds()->SetId(1, id1);
        upSpokes_lines->InsertNextCell(up_arrow);

        // form down spokes
        int id2 = static_cast<int>(downSpokes_pts->InsertNextPoint(mx, my, mz));
        double bx_down = transformed_down_pdm(i,0);
        double by_down = transformed_down_pdm(i,1);
        double bz_down = transformed_down_pdm(i,2);
        int id3 = static_cast<int>(downSpokes_pts->InsertNextPoint(bx_down,by_down,bz_down));

        // spoke length and dir
        vtkVector3d downSpoke(bx_down-mx, by_down-my, bz_down-mz);
        double downSpokeLength = downSpoke.Normalize();
        downSpokeLengths->InsertNextTuple1(downSpokeLength);
        downSpokeDirs->InsertNextTuple3(downSpoke.GetX(), downSpoke.GetY(), downSpoke.GetZ());

        vtkSmartPointer<vtkLine> down_arrow = vtkSmartPointer<vtkLine>::New();
        down_arrow->GetPointIds()->SetId(0, id2);
        down_arrow->GetPointIds()->SetId(1, id3);
        downSpokes_lines->InsertNextCell(down_arrow);

    }
    // display up spokes
    upSpokes_poly->SetPoints(upSpokes_pts);
    upSpokes_poly->SetLines(upSpokes_lines);

    upSpokes_poly->GetPointData()->AddArray(upSpokeDirs);
    upSpokes_poly->GetPointData()->SetActiveVectors("spokeDirection");
    upSpokes_poly->GetPointData()->AddArray(upSpokeLengths);
    upSpokes_poly->GetPointData()->SetActiveScalars("spokeLength");
    AddModelNodeToScene(upSpokes_poly, "up spokes for ellipsoid", true, 0, 1, 1);

    // write to file
    vtkSmartPointer<vtkPolyDataWriter> upSpokeWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
    const std::string upFileName = modelFolder + "/up" + std::to_string(totalNum) + ".vtk";
    const std::string downFileName = modelFolder + "/down" + std::to_string(totalNum) + ".vtk";
    const std::string crestFileName = modelFolder + "/crest" + std::to_string(totalNum) + ".vtk";
    upSpokeWriter->SetFileName(upFileName.c_str());
    upSpokeWriter->SetInputData(upSpokes_poly);
    upSpokeWriter->Update();

    // display down spokes
    downSpokes_poly->SetPoints(downSpokes_pts);
    downSpokes_poly->SetLines(downSpokes_lines);

    downSpokes_poly->GetPointData()->AddArray(downSpokeDirs);
    downSpokes_poly->GetPointData()->SetActiveVectors("spokeDirection");
    downSpokes_poly->GetPointData()->AddArray(downSpokeLengths);
    downSpokes_poly->GetPointData()->SetActiveScalars("spokeLength");

    AddModelNodeToScene(downSpokes_poly, "down spokes for ellipsoid", true, 1, 0, 1);

    vtkSmartPointer<vtkPolyDataWriter> downSpokeWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
    downSpokeWriter->SetFileName(downFileName.c_str());
    downSpokeWriter->SetInputData(downSpokes_poly);
    downSpokeWriter->Update();

    // deal with skeletal mesh
    for(int i = 0; i < nRows * nCols; ++i)
    {
        double mx = transformed_skeletal_points(i, 0);
        double my = transformed_skeletal_points(i, 1);
        double mz = transformed_skeletal_points(i, 2);
        int current_id = static_cast<int>(skeletal_sheet->InsertNextPoint(mx, my, mz));
        int current_row = static_cast<int>(floor(i / nCols));
        int current_col = i - current_row * nCols;
        if(current_col >= 0 && current_row >= 0
                && current_row < nRows-1 && current_col < nCols - 1)
        {
            vtkSmartPointer<vtkQuad> quad = vtkSmartPointer<vtkQuad>::New();
            quad->GetPointIds()->SetId(0, current_id);
            quad->GetPointIds()->SetId(1, current_id + nCols);
            quad->GetPointIds()->SetId(2, current_id + nCols + 1);
            quad->GetPointIds()->SetId(3, current_id + 1);
            skeletal_mesh->InsertNextCell(quad);
        }
    }
    srep_poly->SetPoints(skeletal_sheet);
    srep_poly->SetPolys(skeletal_mesh);
    AddModelNodeToScene(srep_poly, "skeletal mesh for ellipsoid", true, 0, 0, 0);

    const std::string meshFileName = modelFolder + "/mesh" + std::to_string(totalNum) + ".vtk";
    vtkSmartPointer<vtkPolyDataWriter> meshWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
    meshWriter->SetFileName(meshFileName.c_str());
    meshWriter->SetInputData(srep_poly);
    meshWriter->Update();
    // deal with crest spokes
    for(int i = 0; i < nCrestPoints; ++i)
    {
        // tail point
        double cx_t = transformed_crest_base(i, 0);
        double cy_t = transformed_crest_base(i, 1);
        double cz_t = transformed_crest_base(i, 2);
        // head point (_b means boundary)
        double cx_b = transformed_crest_pdm(i, 0);
        double cy_b = transformed_crest_pdm(i, 1);
        double cz_b = transformed_crest_pdm(i, 2);

        if(shift > 0)
        {
            double shift_x = (cx_b - cx_t) * shift;
            double shift_y = (cy_b - cy_t) * shift;
            double shift_z = (cz_b - cz_t) * shift;

            cx_t += shift_x;
            cy_t += shift_y;
            cz_t += shift_z;
        }
        if(i == nCols + 2 * (nRows / 2 - 1)) {
            vtkSmartPointer<vtkConeSource> cone = vtkSmartPointer<vtkConeSource>::New();
            cone->SetCenter(cx_b, cy_b, cz_b);
            cone->SetDirection(cx_b-cx_t, cy_b-cy_t, cz_b-cz_t);
            cone->Update();
            AddModelNodeToScene(cone->GetOutput(), arrowName.c_str(), true, 1, 0, 0);

        }
        if(i == nCols / 2) {
            vtkSmartPointer<vtkConeSource> cone = vtkSmartPointer<vtkConeSource>::New();
            cone->SetCenter(cx_b, cy_b, cz_b);
            cone->SetDirection(cx_b-cx_t, cy_b-cy_t, cz_b-cz_t);
            cone->Update();
            AddModelNodeToScene(cone->GetOutput(), arrowName.c_str(), true, 0, 1, 0);

        }

        int id0 = static_cast<int>(crestSpokes_pts->InsertNextPoint(cx_t, cy_t, cz_t));
        int id1 = static_cast<int>(crestSpokes_pts->InsertNextPoint(cx_b, cy_b, cz_b));

        vtkSmartPointer<vtkLine> crest_arrow = vtkSmartPointer<vtkLine>::New();
        crest_arrow->GetPointIds()->SetId(0, id0);
        crest_arrow->GetPointIds()->SetId(1, id1);
        crestSpokes_lines->InsertNextCell(crest_arrow);

        vtkVector3d crestSpoke(cx_b-cx_t, cy_b-cy_t, cz_b-cz_t);
        double crestSpokeLength = crestSpoke.Normalize();

        crestSpokeLengths->InsertNextTuple1(crestSpokeLength);
        crestSpokeDirs->InsertNextTuple3(crestSpoke.GetX(), crestSpoke.GetY(), crestSpoke.GetZ());
    }
    crestSpokes_poly->SetPoints(crestSpokes_pts);
    crestSpokes_poly->SetLines(crestSpokes_lines);

    crestSpokes_poly->GetPointData()->AddArray(crestSpokeDirs);
    crestSpokes_poly->GetPointData()->SetActiveVectors("spokeDirection");
    crestSpokes_poly->GetPointData()->AddArray(crestSpokeLengths);
    crestSpokes_poly->GetPointData()->SetActiveScalars("spokeLength");

    AddModelNodeToScene(crestSpokes_poly, "crest spokes for ellipsoid", true, 1, 0, 0);
    vtkSmartPointer<vtkPolyDataWriter> crestSpokeWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
    crestSpokeWriter->SetFileName(crestFileName.c_str());
    crestSpokeWriter->SetInputData(crestSpokes_poly);
    crestSpokeWriter->Update();

    // deal with fold curve
    for(int i = 0; i < nCrestPoints; ++i)
    {
        double cx_t = transformed_crest_base(i, 0);
        double cy_t = transformed_crest_base(i, 1);
        double cz_t = transformed_crest_base(i, 2);
        double cx_b = transformed_crest_pdm(i, 0);
        double cy_b = transformed_crest_pdm(i, 1);
        double cz_b = transformed_crest_pdm(i, 2);

        if(shift > 0)
        {
            double shift_x = (cx_b - cx_t) * shift;
            double shift_y = (cy_b - cy_t) * shift;
            double shift_z = (cz_b - cz_t) * shift;

            cx_t += shift_x;
            cy_t += shift_y;
            cz_t += shift_z;
        }
        int id0 = static_cast<int>(foldCurve_pts->InsertNextPoint(cx_t, cy_t, cz_t));

        if(id0 > 0 && i < nCols)
        {
            // first row
            vtkSmartPointer<vtkLine> fold_seg = vtkSmartPointer<vtkLine>::New();
            fold_seg->GetPointIds()->SetId(0, id0-1);
            fold_seg->GetPointIds()->SetId(1, id0);
            fold_curve->InsertNextCell(fold_seg);
        }

        if(i > nCols && i < nCols + 2*(nRows-2) + 1 && (i-nCols) % 2 == 1)
        {
            // right side of crest
            vtkSmartPointer<vtkLine> fold_seg = vtkSmartPointer<vtkLine>::New();
            fold_seg->GetPointIds()->SetId(0, id0-2);
            fold_seg->GetPointIds()->SetId(1, id0);
            fold_curve->InsertNextCell(fold_seg);
        }
        if(i > nCols && i < nCols + 2*(nRows-2) + 1 && (i-nCols) % 2 == 0)
        {
            // part of left side
            vtkSmartPointer<vtkLine> fold_seg = vtkSmartPointer<vtkLine>::New();
            fold_seg->GetPointIds()->SetId(0, id0-2);
            fold_seg->GetPointIds()->SetId(1, id0);
            fold_curve->InsertNextCell(fold_seg);
        }

        if(i == nCols)
        {
            // remaining part of left side
            vtkSmartPointer<vtkLine> fold_seg = vtkSmartPointer<vtkLine>::New();
            fold_seg->GetPointIds()->SetId(0, 0);
            fold_seg->GetPointIds()->SetId(1, id0);
            fold_curve->InsertNextCell(fold_seg);
        }
        if(i > nCols + 2*(nRows-2))
        {
            //bottom side
            vtkSmartPointer<vtkLine> fold_seg = vtkSmartPointer<vtkLine>::New();
            fold_seg->GetPointIds()->SetId(0, id0-1);
            fold_seg->GetPointIds()->SetId(1, id0);
            fold_curve->InsertNextCell(fold_seg);
        }
        if(i == nCrestPoints - 1)
        {
            // bottome right
            vtkSmartPointer<vtkLine> fold_seg = vtkSmartPointer<vtkLine>::New();
            fold_seg->GetPointIds()->SetId(0, id0-nCols);
            fold_seg->GetPointIds()->SetId(1, id0);
            fold_curve->InsertNextCell(fold_seg);
        }
    }
    foldCurve_poly->SetPoints(foldCurve_pts);
    foldCurve_poly->SetLines(fold_curve);
    AddModelNodeToScene(foldCurve_poly, "fold curve of ellipsoid", true, 1, 1, 0);

    const std::string curveFileName = modelFolder + "/curve" + std::to_string(totalNum) + ".vtk";
    vtkSmartPointer<vtkPolyDataWriter> curveWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
    curveWriter->SetFileName(curveFileName.c_str());
    curveWriter->SetInputData(foldCurve_poly);
    curveWriter->Update();
}

void vtkSlicerSkeletalRepresentationInitializerLogic::HideNodesByNameByClass(const std::string & nodeName, const std::string &className)
{
    std::cout << "node name:" << nodeName << std::endl;
    std::cout << "class name:" << className << std::endl;
    std::vector<vtkMRMLNode*> vectModelNodes;
    vtkSmartPointer<vtkCollection> modelNodes = this->GetMRMLScene()->GetNodesByClassByName(className.c_str(), nodeName.c_str());
    modelNodes->InitTraversal();
    for(int i = 0; i < modelNodes->GetNumberOfItems(); i++)
    {
        vtkSmartPointer<vtkMRMLModelNode> thisModelNode = vtkMRMLModelNode::SafeDownCast(modelNodes->GetNextItemAsObject());
        vtkSmartPointer<vtkMRMLModelDisplayNode> displayNode;
        displayNode = thisModelNode->GetModelDisplayNode();
        if(displayNode == nullptr)
        {
            continue;
        }

        displayNode->SetVisibility(0);

    }

}

void vtkSlicerSkeletalRepresentationInitializerLogic::HideNodesByClass(const std::string &className)
{
    vtkSmartPointer<vtkCollection> modelNodes = this->GetMRMLScene()->GetNodesByClass(className.c_str());
    modelNodes->InitTraversal();
    for(int i = 0; i < modelNodes->GetNumberOfItems(); i++)
    {
        vtkSmartPointer<vtkMRMLModelNode> thisModelNode = vtkMRMLModelNode::SafeDownCast(modelNodes->GetNextItemAsObject());
        vtkSmartPointer<vtkMRMLModelDisplayNode> displayNode;
        displayNode = thisModelNode->GetModelDisplayNode();
        if(displayNode == nullptr)
        {
            continue;
        }

        displayNode->SetVisibility(0);

    }

}

int vtkSlicerSkeletalRepresentationInitializerLogic::InklingFlow(const std::string &filename, double dt, double smooth_amount, int max_iter, int freq_output, double /*threshold*/)
{
    //std::cout << threshold << std::endl;
    std::cout << filename << std::endl;
    std::cout << dt << std::endl;
    std::cout << smooth_amount << std::endl;
    std::cout << max_iter << std::endl;
    std::cout << freq_output << std::endl;

    std::string tempFileName(this->GetApplicationLogic()->GetTemporaryPath());
    tempFileName += "/temp_output.vtk";

    vtkSmartPointer<vtkPolyDataReader> reader =
        vtkSmartPointer<vtkPolyDataReader>::New();
    reader->SetFileName(tempFileName.c_str());
    reader->Update();

    vtkSmartPointer<vtkPolyData> mesh =
        vtkSmartPointer<vtkPolyData>::New();
    mesh = reader->GetOutput();

    vtkSmartPointer<vtkMassProperties> mass_filter =
        vtkSmartPointer<vtkMassProperties>::New();
    mass_filter->SetInputData(mesh);
    mass_filter->Update();
    double original_volume = mass_filter->GetVolume();

    int iter = 0;
    double tolerance = 0.05;
    double q = 1.0;

    while(q > tolerance && iter < max_iter)
    {
        // smooth filter
        vtkSmartPointer<vtkWindowedSincPolyDataFilter> smooth_filter =
            vtkSmartPointer<vtkWindowedSincPolyDataFilter>::New();
        smooth_filter->SetPassBand(smooth_amount);
        smooth_filter->NonManifoldSmoothingOn();
        smooth_filter->NormalizeCoordinatesOn();
        smooth_filter->SetNumberOfIterations(20);
        smooth_filter->FeatureEdgeSmoothingOff();
        smooth_filter->BoundarySmoothingOff();
        smooth_filter->SetInputData(mesh);
        smooth_filter->Update();
        if(smooth_amount > 0) {
            mesh = smooth_filter->GetOutput();
        }

        // normal filter
        vtkSmartPointer<vtkPolyDataNormals> normal_filter =
            vtkSmartPointer<vtkPolyDataNormals>::New();
        normal_filter->SplittingOff();
        normal_filter->ComputeCellNormalsOff();
        normal_filter->ComputePointNormalsOn();
        normal_filter->SetInputData(mesh);
        normal_filter->Update();
        vtkDataArray* N = normal_filter->GetOutput()->GetPointData()->GetNormals();
        if(N == nullptr) {
            std::cerr << "error in getting normals" << std::endl;
            return EXIT_FAILURE;
        }

        // mean curvature filter
        vtkSmartPointer<vtkCurvatures> curvature_filter =
            vtkSmartPointer<vtkCurvatures>::New();
        curvature_filter->SetCurvatureTypeToMean();
        curvature_filter->SetInputData(mesh);
        curvature_filter->Update();

        vtkSmartPointer<vtkDoubleArray> H =
            vtkDoubleArray::SafeDownCast(curvature_filter->GetOutput()->GetPointData()->GetArray("Mean_Curvature"));


        if(H == nullptr) {
            std::cerr << "error in getting mean curvature" << std::endl;
            return EXIT_FAILURE;
        }

        curvature_filter->SetCurvatureTypeToGaussian();
        curvature_filter->Update();
        vtkSmartPointer<vtkDoubleArray> K =
            vtkDoubleArray::SafeDownCast(curvature_filter->GetOutput()->GetPointData()->GetArray("Gauss_Curvature"));

        if(K == nullptr) {
            std::cerr << "error in getting Gaussian curvature" << std::endl;
            return EXIT_FAILURE;
        }

        curvature_filter->SetCurvatureTypeToMaximum();
        curvature_filter->Update();
        vtkSmartPointer<vtkDoubleArray> MC =
            vtkDoubleArray::SafeDownCast(curvature_filter->GetOutput()->GetPointData()->GetArray("Maximum_Curvature"));

        if(MC == nullptr) {
            std::cerr << "error in getting max curvature" << std::endl;
            return EXIT_FAILURE;
        }

        curvature_filter->SetCurvatureTypeToMinimum();
        curvature_filter->Update();
        vtkSmartPointer<vtkDoubleArray> MinC =
            vtkDoubleArray::SafeDownCast(curvature_filter->GetOutput()->GetPointData()->GetArray("Minimum_Curvature"));
        if(MinC == nullptr)
        {
            std::cout << "error in getting min curvature" << std::endl;
            return -1;
        }

        // perform the flow
        vtkSmartPointer<vtkPoints> points = mesh->GetPoints();
        double maxVal = -10000.0;
        vtkIdType maxIndex = -1;
        for(int i = 0; i < points->GetNumberOfPoints(); ++i) {
            double p[3];
            points->GetPoint(i, p);
            double curr_N[3];
            N->GetTuple(i, curr_N);
            double curr_H = H->GetValue(i);
            //double curr_K = K->GetValue(i);
            double curr_max = MC->GetValue(i);
            double curr_min = MinC->GetValue(i);

            double delta = 0.0;
            double diffK = curr_max - curr_min;
            delta = dt * curr_H; // * 1 / (diffK * diffK);
            diffK = std::abs(diffK);

            if(diffK > maxVal)
            {
                maxVal = diffK;
                maxIndex = i;
            }
            if(curr_max >= 0 && curr_min >= 0)
            {
                delta = dt * curr_max;
            }
            else if(curr_max < 0 && curr_min < 0)
            {
                delta = dt * curr_min;
            }
            else
            {
                delta = dt * curr_H;
            }
            for(int idx = 0; idx < 3; ++idx)
            {
                p[idx] -= delta * curr_N[idx]; //dt * curr_H * curr_N[idx];
            }
            points->SetPoint(i, p);
        }
        points->Modified();

        mass_filter->SetInputData(mesh);
        mass_filter->Update();
        double curr_volume = mass_filter->GetVolume();

        for(int i = 0; i < points->GetNumberOfPoints(); ++i) {
            double p[3];
            points->GetPoint(i, p);
            for(int j = 0; j < 3; ++j) {
                p[j] *= std::pow( original_volume / curr_volume , 1.0 / 3.0 );
            }
        }
        points->Modified();

        double testRender[3];
        points->GetPoint(maxIndex, testRender);
//        AddPointToScene(testRender[0], testRender[1], testRender[2], 13); // sphere3D

        HideNodesByNameByClass("output_inkling","vtkMRMLModelNode");

         if((iter +1) % freq_output == 0)
         {
             std::string modelName("output_inkling");
             modelName += std::to_string(iter+1);
             AddModelNodeToScene(mesh, modelName.c_str(), false);
         }

        q -= 0.0001;
        iter++;
    }

    return 1;

}

void vtkSlicerSkeletalRepresentationInitializerLogic::AddPointToScene(double x, double y, double z, int glyphType, double r, double g, double b)
{
    std::cout << "AddPointToScene: parameters:" << x << ", y:" << y << ", z:" << z << std::endl;
    vtkMRMLScene *scene = this->GetMRMLScene();
    if(!scene)
    {
        vtkErrorMacro(" Invalid scene");
        return;
    }

    // node which controls display properties
    vtkSmartPointer<vtkMRMLMarkupsDisplayNode> displayNode;
    displayNode = vtkSmartPointer<vtkMRMLMarkupsDisplayNode>::New();
//    this->SetDisplayNodeToDefaults(displayNode);
    displayNode->SetGlyphScale(1.10);
    displayNode->SetTextScale(3.40);
    displayNode->SetSelectedColor(r, g, b);

    displayNode->SetGlyphType(glyphType); // 13: sphere3D
    scene->AddNode(displayNode);

    // model node
    vtkSmartPointer<vtkMRMLMarkupsFiducialNode> fidNode;

    fidNode = vtkSmartPointer<vtkMRMLMarkupsFiducialNode>::New();
    if(fidNode == nullptr)
    {
        vtkErrorMacro("fidNode is NULL");
        return;
    }
    fidNode->SetAndObserveDisplayNodeID(displayNode->GetID());
    fidNode->SetLocked(true);
    fidNode->SetName("Hi");
    scene->AddNode(fidNode);

    fidNode->AddFiducial(x, y, z);

}

// compute and apply tps transformation matrix
void vtkSlicerSkeletalRepresentationInitializerLogic::ComputePairwiseTps(int totalNum)
{

    typedef itkThinPlateSplineExtended TransformType;
    typedef itk::Point< CoordinateRepType, 3 > PointType;
    typedef TransformType::PointSetType PointSetType;
    typedef PointSetType::PointIdentifier PointIdType;

    // create folder if not exist

    const std::string tempFolder(this->GetApplicationLogic()->GetTemporaryPath());
    const std::string backwardFolder = tempFolder + "/backward";
    std::cout << "backward folder" << backwardFolder << std::endl;
    if (!vtksys::SystemTools::FileExists(backwardFolder, false))
    {
      if (!vtksys::SystemTools::MakeDirectory(backwardFolder))
      {
            std::cout << "Failed to create folder : " << backwardFolder << std::endl;

      }
    }

    // compute transformation matrix by current surface and backward surface

    vtkSmartPointer<vtkPolyDataReader> sourceSurfaceReader = vtkSmartPointer<vtkPolyDataReader>::New();
    vtkSmartPointer<vtkPolyDataReader> targetSurfaceReader = vtkSmartPointer<vtkPolyDataReader>::New();

    std::cout << "Computing pairwise transformation matrix via TPS for " << totalNum << " cases..." << std::endl;
    for(int stepNum = totalNum; stepNum > 1; --stepNum)
    {
        const std::string inputMeshFile = tempFolder + "/forward/" + std::to_string(stepNum) + ".vtk";
        const std::string nextMeshFile = tempFolder + "/forward/" + std::to_string(stepNum - 1) + ".vtk";
        sourceSurfaceReader->SetFileName(inputMeshFile.c_str());
        sourceSurfaceReader->Update();

        // current surface mesh
        vtkSmartPointer<vtkPolyData> polyData_source = sourceSurfaceReader->GetOutput();

        targetSurfaceReader->SetFileName(nextMeshFile.c_str());
        targetSurfaceReader->Update();
        // next surface mesh which back flow to
        vtkSmartPointer<vtkPolyData> polyData_target = targetSurfaceReader->GetOutput();

        PointSetType::Pointer sourceLandMarks = PointSetType::New();
        PointSetType::Pointer targetLandMarks = PointSetType::New();
        PointType p1; PointType p2; // same as double p1[3];
        PointSetType::PointsContainer::Pointer sourceLandMarkContainer
                = sourceLandMarks->GetPoints();
        PointSetType::PointsContainer::Pointer targetLandMarkContainer
                = targetLandMarks->GetPoints();

        PointIdType id_s = itk::NumericTraits< PointIdType >::Zero;
        PointIdType id_t = itk::NumericTraits< PointIdType >::Zero;

        // Read in the source points set
        for(unsigned int i = 0; i < polyData_source->GetNumberOfPoints(); i += 10){
            double p[3];
            polyData_source->GetPoint(i,p);
            p1[0] = p[0];
            p1[1] = p[1];
            p1[2] = p[2];
            sourceLandMarkContainer->InsertElement(id_s, p1);
            id_s++;
        }

        // Read in the target points set
        for(unsigned int i = 0; i < polyData_target->GetNumberOfPoints(); i += 10){
            double p[3];
            polyData_target->GetPoint(i,p);
            p2[0] = p[0];
            p2[1] = p[1];
            p2[2] = p[2];
            targetLandMarkContainer->InsertElement(id_t, p2);
            id_t++;
        }

        TransformType::Pointer tps = TransformType::New();
        tps->SetSourceLandmarks(sourceLandMarks);
        tps->SetTargetLandmarks(targetLandMarks);
        tps->ComputeWMatrix();

        // Apply tps on the srep in future (from the ellipsoid)
        // srep of totalNum == srep of ellipsoid
//        std::cout << "Applying transformation matrix on from step:" << stepNum << " to " << stepNum-1 << std::endl;
        std::string upFileName = tempFolder + "/model/up" + std::to_string(stepNum) + ".vtk";
        std::string downFileName = tempFolder + "/model/down" + std::to_string(stepNum) + ".vtk";
        std::string crestFileName = tempFolder + "/model/crest" + std::to_string(stepNum) + ".vtk";
        // read source srep
        vtkSmartPointer<vtkPolyDataReader> upSpokeReader = vtkSmartPointer<vtkPolyDataReader>::New();
        upSpokeReader->SetFileName(upFileName.c_str());
        upSpokeReader->Update();
        vtkSmartPointer<vtkPolyData> upSpokes = upSpokeReader->GetOutput();

        vtkSmartPointer<vtkPolyDataReader> downSpokeReader = vtkSmartPointer<vtkPolyDataReader>::New();
        downSpokeReader->SetFileName(downFileName.c_str());
        downSpokeReader->Update();
        vtkSmartPointer<vtkPolyData> downSpokes = downSpokeReader->GetOutput();

        vtkSmartPointer<vtkPolyDataReader> crestSpokeReader = vtkSmartPointer<vtkPolyDataReader>::New();
        crestSpokeReader->SetFileName(crestFileName.c_str());
        crestSpokeReader->Update();
        vtkSmartPointer<vtkPolyData> crestSpokes = crestSpokeReader->GetOutput();

        upFileName = tempFolder + "/model/up" + std::to_string(stepNum - 1) + ".vtk";
        std::string upOutputFile(upFileName.c_str());
        TransformNOutput(tps, upSpokes, upOutputFile);
        // display upspokes

        downFileName = tempFolder + "/model/down" + std::to_string(stepNum - 1) + ".vtk";
        std::string downOutputFile(downFileName.c_str());
        TransformNOutput(tps, downSpokes, downOutputFile);

        crestFileName = tempFolder + "/model/crest" + std::to_string(stepNum - 1) + ".vtk";
        std::string crestOutputFile(crestFileName.c_str());
        TransformNOutput(tps, crestSpokes, crestOutputFile);

        // Apply transformation matrix on the skeletal sheet and fold curve
        std::string meshFileName(tempFolder);
        meshFileName = tempFolder + "/model/mesh" + std::to_string(stepNum) + ".vtk";// Transform input
        vtkSmartPointer<vtkPolyDataReader> meshReader = vtkSmartPointer<vtkPolyDataReader>::New();
        meshReader->SetFileName(meshFileName.c_str());
        meshReader->Update();
        vtkSmartPointer<vtkPolyData> meshPoly = meshReader->GetOutput();
        meshFileName = tempFolder + "/model/mesh" + std::to_string(stepNum-1) + ".vtk";// output
        TransformPoints(tps, meshPoly, meshFileName);

        std::string curveFileName(tempFolder);
        curveFileName = tempFolder + "/model/curve" + std::to_string(stepNum) + ".vtk"; // Transform input
        vtkSmartPointer<vtkPolyDataReader> curveReader = vtkSmartPointer<vtkPolyDataReader>::New();
        curveReader->SetFileName(curveFileName.c_str());
        curveReader->Update();
        vtkSmartPointer<vtkPolyData> curvePoly = curveReader->GetOutput();
        curveFileName = tempFolder + "/model/curve" + std::to_string(stepNum-1) + ".vtk"; // output
        TransformPoints(tps, curvePoly, curveFileName);
    }
}

void vtkSlicerSkeletalRepresentationInitializerLogic::BackwardFlow(int totalNum)
{

    // 1. compute pairwise TPS
    ComputePairwiseTps(totalNum);
    std::cout << "Finished computing transformation matrix." << std::endl;

    std::cout << "Finished applying transformation matrix." << std::endl;

    // 2. display the srep for the initial object
    DisplayResultSrep();
}

void vtkSlicerSkeletalRepresentationInitializerLogic::SetOutputPath(const std::string &outputPath)
{
    mOutputPath = outputPath;
}

void vtkSlicerSkeletalRepresentationInitializerLogic::SetRows(int r)
{

    mRows = r;
}

void vtkSlicerSkeletalRepresentationInitializerLogic::SetCols(int c)
{
    mCols = c;
}

void vtkSlicerSkeletalRepresentationInitializerLogic::DisplayResultSrep(bool flip)
{
    // Hide other nodes.
    HideNodesByClass("vtkMRMLModelNode");
    const std::string tempFolder(this->GetApplicationLogic()->GetTemporaryPath());
    std::string upFileName = tempFolder + "/model/up1.vtk";
    std::string downFileName = tempFolder + "/model/down1.vtk";
    const std::string crestFileName = tempFolder +"/model/crest1.vtk";
    const std::string meshFileName = tempFolder + "/model/mesh1.vtk";;
    const std::string curveFileName = tempFolder + "/model/curve1.vtk";
    if(flip) {
        upFileName = downFileName;
        downFileName = tempFolder + "/model/up1.vtk";
    }

    vtkSmartPointer<vtkPolyDataReader> upSpokeReader = vtkSmartPointer<vtkPolyDataReader>::New();
    upSpokeReader->SetFileName(upFileName.c_str());
    upSpokeReader->Update();
    vtkSmartPointer<vtkPolyData> upSpoke_poly = upSpokeReader->GetOutput();

    vtkSmartPointer<vtkPolyDataReader> downSpokeReader = vtkSmartPointer<vtkPolyDataReader>::New();
    downSpokeReader->SetFileName(downFileName.c_str());
    downSpokeReader->Update();
    vtkSmartPointer<vtkPolyData> downSpoke_poly = downSpokeReader->GetOutput();

    vtkSmartPointer<vtkPolyDataReader> crestSpokeReader = vtkSmartPointer<vtkPolyDataReader>::New();
    crestSpokeReader->SetFileName(crestFileName.c_str());
    crestSpokeReader->Update();
    vtkSmartPointer<vtkPolyData> crestSpoke_poly = crestSpokeReader->GetOutput();

    vtkSmartPointer<vtkPolyDataReader> meshReader = vtkSmartPointer<vtkPolyDataReader>::New();
    meshReader->SetFileName(meshFileName.c_str());
    meshReader->Update();
    vtkSmartPointer<vtkPolyData> meshPoly = meshReader->GetOutput();

    vtkSmartPointer<vtkPolyDataReader> curveReader = vtkSmartPointer<vtkPolyDataReader>::New();
    curveReader->SetFileName(curveFileName.c_str());
    curveReader->Update();
    vtkSmartPointer<vtkPolyData> curvePoly = curveReader->GetOutput();

    AddModelNodeToScene(upSpoke_poly, "up spokes for initial object", true, 0, 1, 1);
    AddModelNodeToScene(downSpoke_poly, "down spokes for initial object", true, 1, 0, 1);
    AddModelNodeToScene(crestSpoke_poly, "crest spokes for initial object", true, 1, 0, 0);
    AddModelNodeToScene(meshPoly, "skeletal mesh for initial object", true, 0, 0, 0);
    AddModelNodeToScene(curvePoly, "fold curve for initial object", true, 1, 1, 0);

    // add spokeLength and spokeDirection to polyData
    vtkSmartPointer<vtkPolyData> vtpUpSpoke = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPolyData> vtpDownSpoke = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPolyData> vtpCrestSpoke = vtkSmartPointer<vtkPolyData>::New();
    CompletePolyData(upSpoke_poly, vtpUpSpoke);
    CompletePolyData(downSpoke_poly, vtpDownSpoke);
    CompletePolyData(crestSpoke_poly, vtpCrestSpoke, true);
    // save to files
    std::string outputUpFileName = mOutputPath + "/up.vtp";
    std::string outputDownFileName = mOutputPath + "/down.vtp";
    std::string outputCrestFileName = mOutputPath + "/crest.vtp";

    vtkSmartPointer<vtkXMLPolyDataWriter> vtpWriter = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    vtpWriter->SetDataModeToAscii();
    vtpWriter->SetFileName(outputUpFileName.c_str());
    vtpWriter->SetInputData(vtpUpSpoke);
    vtpWriter->Update();

    vtpWriter->SetFileName(outputDownFileName.c_str());
    vtpWriter->SetInputData(vtpDownSpoke);
    vtpWriter->Update();

    vtpWriter->SetFileName(outputCrestFileName.c_str());
    vtpWriter->SetInputData(vtpCrestSpoke);
    vtpWriter->Update();

    // compose header file
    std::stringstream output;

    output<<"<s-rep>"<<std::endl;
    output<<"  <nRows>"<<mRows<<"</nRows>"<<std::endl;
    output<<"  <nCols>"<<mCols<<"</nCols>"<<std::endl;
    output<<"  <meshType>Quad</meshType>"<< std::endl;
    output<<"  <color>"<<std::endl;
    output<<"    <red>0</red>"<<std::endl;
    output<<"    <green>0.5</green>"<<std::endl;
    output<<"    <blue>0</blue>"<<std::endl;
    output<<"  </color>"<<std::endl;
    output<<"  <isMean>False</isMean>"<<std::endl;
    output<<"  <meanStatPath/>"<<std::endl;
    output<<"  <upSpoke>"<< outputUpFileName<<"</upSpoke>"<<std::endl;
    output<<"  <downSpoke>"<< outputDownFileName << "</downSpoke>"<<std::endl;
    output<<"  <crestSpoke>"<< outputCrestFileName << "</crestSpoke>"<<std::endl;
    output<<"</s-rep>"<<std::endl;

    std::string header_file = mOutputPath + "/header.xml"; ;
    std::ofstream out_file;
    out_file.open(header_file);
    out_file << output.rdbuf();
    out_file.close();
}

void vtkSlicerSkeletalRepresentationInitializerLogic::RotateSkeleton(bool rotateX, bool rotateY, bool rotateZ)
{
    const std::string tempFolder(this->GetApplicationLogic()->GetTemporaryPath());
    const std::string newEllSurfaceFile = tempFolder + "/forward/" + std::to_string(forwardCount) + ".vtk";

    vtkSmartPointer<vtkPolyDataReader> reader =
        vtkSmartPointer<vtkPolyDataReader>::New();
    reader->SetFileName(newEllSurfaceFile.c_str());
    reader->Update();

    vtkSmartPointer<vtkPolyData> mesh =
        vtkSmartPointer<vtkPolyData>::New();
    mesh = reader->GetOutput();
    GenerateSrepForEllipsoid(mesh, mRows, mCols, forwardCount, rotateX, rotateY, rotateZ);
}

void vtkSlicerSkeletalRepresentationInitializerLogic::ReorderSpokes(vtkPolyData *input, vtkPoints* outputPts, vtkCellArray* outputPolys)
{
    if(input->GetNumberOfPoints() == 0) return;
    outputPts->SetDataTypeToDouble();

    vtkSmartPointer<vtkPolyData> testPoly = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints> testPts = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> testPolys = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPolyData> test2Poly = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints> test2Pts = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> test2Polys = vtkSmartPointer<vtkCellArray>::New();


    for(int i = 0; i < mRows; ++i) {
        double mx[mCols];
        double my[mCols];
        double mz[mCols];
        double bx[mCols];
        double by[mCols];
        double bz[mCols];
        for(int j = 0; j < mCols; ++j) {
            double skeletalPt[3];
            double bdryPt[3];
            int idx = i * mCols + j;
            input->GetPoint(idx*2, skeletalPt);
            input->GetPoint(idx*2+1, bdryPt);
            mx[j] = skeletalPt[0];
            my[j] = skeletalPt[1];
            mz[j] = skeletalPt[2];
            bx[j] = bdryPt[0];
            by[j] = bdryPt[1];
            bz[j] = bdryPt[2];
        }

        for(int j = mCols - 1; j >=0; --j) {
            double newSkeletalPt[3], newBdryPt[3];
            newSkeletalPt[0] = mx[j];
            newSkeletalPt[1] = my[j];
            newSkeletalPt[2] = mz[j];
            vtkIdType id0 = outputPts->InsertNextPoint(newSkeletalPt);

            newBdryPt[0] = bx[j];
            newBdryPt[1] = by[j];
            newBdryPt[2] = bz[j];
            vtkIdType id1 = outputPts->InsertNextPoint(newBdryPt);

            vtkSmartPointer<vtkLine> arrow = vtkSmartPointer<vtkLine>::New();
            arrow->GetPointIds()->SetId(0, id0);
            arrow->GetPointIds()->SetId(1, id1);
            outputPolys->InsertNextCell(arrow);

            if(i == 0 || i == mRows - 1) {
                double skeletalPt[3];
                double bdryPt[3];
                int idx = i * mCols + j;
                input->GetPoint(idx*2, skeletalPt);
                input->GetPoint(idx*2+1, bdryPt);
                vtkIdType id0Test = testPts->InsertNextPoint(skeletalPt);
                vtkIdType id1Test = testPts->InsertNextPoint(bdryPt);
                vtkSmartPointer<vtkLine> arrowTest = vtkSmartPointer<vtkLine>::New();
                arrowTest->GetPointIds()->SetId(0, id0Test);
                arrowTest->GetPointIds()->SetId(1, id1Test);
                testPolys->InsertNextCell(arrowTest);

                vtkIdType id0Test2 = test2Pts->InsertNextPoint(newSkeletalPt);
                vtkIdType id1Test2 = test2Pts->InsertNextPoint(newBdryPt);

                vtkSmartPointer<vtkLine> arrowTest2 = vtkSmartPointer<vtkLine>::New();
                arrowTest2->GetPointIds()->SetId(0, id0Test2);
                arrowTest2->GetPointIds()->SetId(1, id1Test2);
                test2Polys->InsertNextCell(arrowTest2);

            }
        }
    }
    testPts->Modified();
    testPoly->Modified();
    test2Pts->Modified();
    test2Polys->Modified();
    testPoly->SetPoints(testPts);
    testPoly->SetPolys(testPolys);
    test2Poly->SetPoints(test2Pts);
    test2Poly->SetPolys(test2Polys);

    outputPts->Modified();
    outputPolys->Modified();

//    AddModelNodeToScene(testPoly, "up spokes for initial object", true, 0, 1, 1);
//    AddModelNodeToScene(test2Poly, "reordered up spokes for initial object", true, 1, 0, 1);

}

void vtkSlicerSkeletalRepresentationInitializerLogic::TransformNOutput(itkThinPlateSplineExtended::Pointer tps, vtkPolyData* spokes, const std::string& outputFileName)
{
    vtkPoints* newPoints = vtkPoints::New();
    newPoints->SetDataTypeToDouble();
    vtkSmartPointer<vtkPolyData> newSpokePoly = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkCellArray> spokeLines     = vtkSmartPointer<vtkCellArray>::New();
    for(int i = 0; i < spokes->GetNumberOfPoints(); i += 2)
    {
        double p[3];
        spokes->GetPoint(i,p);

        // transform medial point by tps
        PointType hub;
        hub[0] = p[0];
        hub[1] = p[1];
        hub[2] = p[2];
        PointType transHub = tps->TransformPoint(hub);
        double newP[3];
        newP[0] = transHub[0];
        newP[1] = transHub[1];
        newP[2] = transHub[2];
        int id0 = static_cast<int>(newPoints->InsertNextPoint(newP));

        // transform implied boundary point by tps
        double p_bdry[3];
        spokes->GetPoint(i+1, p_bdry);
        PointType bdry;
        bdry[0] = p_bdry[0];
        bdry[1] = p_bdry[1];
        bdry[2] = p_bdry[2];
        PointType transB = tps->TransformPoint(bdry);
        double newB[3];
        newB[0] = transB[0];
        newB[1] = transB[1];
        newB[2] = transB[2];
        int id1 = static_cast<int>(newPoints->InsertNextPoint(newB));

        vtkSmartPointer<vtkLine> arrow = vtkSmartPointer<vtkLine>::New();
        arrow->GetPointIds()->SetId(0, id0);
        arrow->GetPointIds()->SetId(1, id1);
        spokeLines->InsertNextCell(arrow);

    }
    newSpokePoly->SetPoints(newPoints);
    newSpokePoly->SetLines(spokeLines);

    vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetFileName(outputFileName.c_str());
    writer->SetInputData(newSpokePoly);
    writer->Update();
}

void vtkSlicerSkeletalRepresentationInitializerLogic::TransformPoints(itkThinPlateSplineExtended::Pointer tps,
                                                                     vtkPolyData *poly, const std::string &outputFileName)
{
    vtkSmartPointer<vtkPoints> pts = poly->GetPoints();
    for(int i = 0; i < poly->GetNumberOfPoints(); ++i)
    {
        double pt[3];
        poly->GetPoint(i, pt);

        PointType itkPt;
        itkPt[0] = pt[0];
        itkPt[1] = pt[1];
        itkPt[2] = pt[2];
        PointType transPt = tps->TransformPoint(itkPt);
        pts->SetPoint(i, transPt[0], transPt[1], transPt[2]);
    }
    pts->Modified();
    vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetFileName(outputFileName.c_str());
    writer->SetInputData(poly);
    writer->Update();
}
double vtkSlicerSkeletalRepresentationInitializerLogic::CalculateSpokeLength(PointType tail, PointType tip){
    PointType spokeVector;
    for(unsigned int dim = 0; dim < 3; dim++){
        spokeVector[dim] = tip[dim] - tail[dim];
    }

    // Compute the spoke length
    double spokeRadiu=0;
    for(unsigned int dim=0; dim<3; dim++){
        spokeRadiu += spokeVector[dim] * spokeVector[dim];
    }

    return sqrt(spokeRadiu);
}

void vtkSlicerSkeletalRepresentationInitializerLogic::GetNeighborCells(vtkPolyData *mesh, int ptId, int newId, vtkCellArray *output, vtkPoints* morePts)
{
    vtkNew<vtkIdList> cellIdList;
    mesh->GetPointCells(ptId, cellIdList);
    for(vtkIdType i = 0; i < cellIdList->GetNumberOfIds(); i++)
    {

        vtkCell* cell = mesh->GetCell(cellIdList->GetId(i));
        //cout << "The cell has " << cell->GetNumberOfEdges() << " edges." << endl;

        //if the cell doesn't have any edges, it is a line
        if(cell->GetNumberOfEdges() <= 0)
        {
            continue;
        }

        for(int j = 0; j < cell->GetNumberOfEdges(); ++j)
        {
            vtkCell* edge = cell->GetEdge(j);
            vtkIdList* pointIdList = edge->GetPointIds();
            double pt0[3], pt1[3];

            mesh->GetPoint(pointIdList->GetId(0), pt0);
            mesh->GetPoint(pointIdList->GetId(1), pt1);

            vtkIdType id2;
            vtkSmartPointer<vtkLine> newEdge = vtkSmartPointer<vtkLine>::New();
            if(pointIdList->GetId(0) == ptId)
            {
                id2 = morePts->InsertNextPoint(pt1);
                newEdge->GetPointIds()->SetId(0, newId);
                newEdge->GetPointIds()->SetId(1, id2);
            }
            else {
                id2 = morePts->InsertNextPoint(pt0);
                newEdge->GetPointIds()->SetId(1, newId);
                newEdge->GetPointIds()->SetId(0, id2);
            }

            output->InsertNextCell(newEdge);
        }
    }
}

void vtkSlicerSkeletalRepresentationInitializerLogic::CompletePolyData(vtkPolyData *poly, vtkPolyData *output, bool isCrest)
{
    poly->GetLines()->InitTraversal();
    vtkNew<vtkIdList> idList;
    //vtkSmartPointer<vtkIdList> idList = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkDoubleArray> spokeDirection = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> spokeLengths = vtkSmartPointer<vtkDoubleArray>::New();

    spokeLengths->SetNumberOfComponents(1);
    spokeLengths->SetName("spokeLength");

    spokeDirection->SetNumberOfComponents(3);
    spokeDirection->SetName("spokeDirection");
    vtkSmartPointer<vtkPoints> skeletalPts = vtkSmartPointer<vtkPoints>::New();

    while(poly->GetLines()->GetNextCell(idList))
    {
        if(idList->GetNumberOfIds() != 2)
        {
            continue;
        }
        double basePt[3], bdryPt[3], spoke[3];
        poly->GetPoint(idList->GetId(0), basePt);
        poly->GetPoint(idList->GetId(1), bdryPt);
        spoke[0] = bdryPt[0] - basePt[0];
        spoke[1] = bdryPt[1] - basePt[1];
        spoke[2] = bdryPt[2] - basePt[2];

        double r = vtkMath::Normalize(spoke);
        spokeLengths->InsertNextTuple1(r);
        spokeDirection->InsertNextTuple(spoke);

        skeletalPts->InsertNextPoint(basePt);
    }

    output->SetPoints(skeletalPts);

    output->GetPointData()->AddArray(spokeDirection);
    output->GetPointData()->SetActiveVectors("spokeDirection");
    output->GetPointData()->AddArray(spokeLengths);
    output->GetPointData()->SetActiveScalars("spokeLength");

    // connection of skeletal points
    if(isCrest)
    {
        vtkSmartPointer<vtkCellArray> curve = vtkSmartPointer<vtkCellArray>::New();
        for(int i = 0; i < skeletalPts->GetNumberOfPoints() - 1; ++i)
        {
            if(i < mCols-1)
            {
                // Horizontal connection for top row
                vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
                line->GetPointIds()->SetId(0, i);
                line->GetPointIds()->SetId(1, i+1);
                curve->InsertNextCell(line);
            }
            else if(i > mCols + 2 * (mRows-2)-1){
                // Backward horizontal connection for bot row
                vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
                line->GetPointIds()->SetId(0, i+1);
                line->GetPointIds()->SetId(1, i);
                curve->InsertNextCell(line);
            }
            else if((i - mCols) % 2 == 0){
                // vertical connection for left side points
                vtkSmartPointer<vtkLine> line1 = vtkSmartPointer<vtkLine>::New();
                line1->GetPointIds()->SetId(0, i+2);
                line1->GetPointIds()->SetId(1, i);
                curve->InsertNextCell(line1);
            }
            else if(i > mCols-1 && (i - mCols - 1) % 2 == 0){
                // vertical connection for right side points
                vtkSmartPointer<vtkLine> line2 = vtkSmartPointer<vtkLine>::New();
                line2->GetPointIds()->SetId(0, i);
                line2->GetPointIds()->SetId(1, i - 2);
                curve->InsertNextCell(line2);
            }

        }
        vtkSmartPointer<vtkLine> lineEnd = vtkSmartPointer<vtkLine>::New();
        lineEnd->GetPointIds()->SetId(0, mCols);
        lineEnd->GetPointIds()->SetId(1, 0);

        curve->InsertNextCell(lineEnd);
        lineEnd->GetPointIds()->SetId(0, skeletalPts->GetNumberOfPoints() - 1);
        lineEnd->GetPointIds()->SetId(1, skeletalPts->GetNumberOfPoints() - 1 - mCols);

        curve->InsertNextCell(lineEnd);
        output->SetLines(curve);
    }
    else {
        vtkSmartPointer<vtkCellArray> quads = vtkSmartPointer<vtkCellArray>::New();

        for(int i = 0; i < mRows - 1; ++i)
        {
            for(int j = 0; j < mCols-1; ++j)
            {
                int id0 = i * mCols + j;
                int id1 = id0 + 1;
                int id2 = id0 + mCols;
                int id3 = id2 + 1;
                vtkSmartPointer<vtkQuad> quad = vtkSmartPointer<vtkQuad>::New();
                quad->GetPointIds()->SetId(0, id0);
                quad->GetPointIds()->SetId(1, id2);
                quad->GetPointIds()->SetId(2, id3);
                quad->GetPointIds()->SetId(3, id1);
                quads->InsertNextCell(quad);
            }
        }
        output->SetPolys(quads);
    }


}
