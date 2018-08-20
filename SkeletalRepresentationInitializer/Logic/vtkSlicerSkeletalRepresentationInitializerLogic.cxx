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
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkCenterOfMass.h>
#include <vtkObjectFactory.h>
#include <vtkCurvatures.h>
#include <vtkWindowedSincPolyDataFilter.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkParametricEllipsoid.h>
#include <vtkParametricFunctionSource.h>
#include <vtkImplicitPolyDataDistance.h>
#include <vtkPolyDataNormals.h>

#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkMassProperties.h>

// Eigen includes
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

// STD includes
#include <cassert>
#include <iostream>


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
  assert(this->GetMRMLScene() != 0);
}

//---------------------------------------------------------------------------
void vtkSlicerSkeletalRepresentationInitializerLogic::UpdateFromMRMLScene()
{
  assert(this->GetMRMLScene() != 0);
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
int vtkSlicerSkeletalRepresentationInitializerLogic::FlowSurfaceOneStep(double dt, double smooth_amount)
{
    std::cout << "flow one step : dt-" << dt << std::endl;
    std::cout << "flow one step : smooth amount-" << smooth_amount << std::endl;
    std::string name = "temp_output.vtk";
    vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
    reader->SetFileName(name.c_str());
    reader->Update();
    vtkSmartPointer<vtkPolyData> mesh = reader->GetOutput();
    if(mesh == NULL)
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

    vtkSmartPointer<vtkDoubleArray> H =
        vtkDoubleArray::SafeDownCast(curvature_filter->GetOutput()->GetPointData()->GetArray("Mean_Curvature"));
    if(H == NULL) {
        vtkErrorMacro("error in getting mean curvature");
        return -1;
    }
    vtkDataArray* N = normal_filter->GetOutput()->GetPointData()->GetNormals();
    if(N == NULL) {
        vtkErrorMacro("error in getting normals");
        return -1;
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
//        points->SetPoint(i, p);
    }
    points->Modified();

    // firstly get other intermediate result invisible
    HideNodesByNameByClass("curvate_flow_result","vtkMRMLModelNode");
    HideNodesByNameByClass("best_fitting_ellipsoid_polydata", "vtkMRMLModelNode");

    // then add this new intermediate result
    std::string modelName("curvate_flow_result");
    AddModelNodeToScene(mesh, modelName.c_str(), true);

    vtkSmartPointer<vtkPolyDataWriter> writer =
        vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetInputData(mesh);
    writer->SetFileName(name.c_str());
    writer->Update();

    // compute the fitting ellipsoid
    vtkSmartPointer<vtkCenterOfMass> centerMassFilter =
        vtkSmartPointer<vtkCenterOfMass>::New();
    centerMassFilter->SetInputData(mesh);
    centerMassFilter->SetUseScalarsAsWeights(false);
    centerMassFilter->Update();
    double center[3];
    centerMassFilter->GetCenter(center);

    ShowFittingEllipsoid(points, curr_volume, center);
    return 0;
}
int vtkSlicerSkeletalRepresentationInitializerLogic::SetInputFileName(const std::string &filename)
{
    vtkSmartPointer<vtkPolyDataReader> reader =
        vtkSmartPointer<vtkPolyDataReader>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();

    vtkSmartPointer<vtkPolyData> mesh;
    mesh = reader->GetOutput();
    // output the original mesh
    std::string modelName("original");
    AddModelNodeToScene(mesh, modelName.c_str(), true, 0.88, 0.88, 0.88);

    // save
    std::string name = "temp_output.vtk";
    vtkSmartPointer<vtkPolyDataWriter> writer =
        vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetInputData(mesh);
    writer->SetFileName(name.c_str());
    writer->Update();

    return 0;
}

// flow surface to the end: either it's ellipsoidal enough or reach max_iter
int vtkSlicerSkeletalRepresentationInitializerLogic::FlowSurfaceMesh(const std::string &filename, double dt, double smooth_amount, int max_iter, int freq_output)
{
    // std::cout << filename << std::endl;
    // std::cout << dt << std::endl;
    // std::cout << smooth_amount << std::endl;
    // std::cout << max_iter << std::endl;
    // std::cout << freq_output << std::endl;
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
//    std::cout << "Original Volume: " << original_volume << std::endl;

    // default parameters
    // double dt = 0.001;
    // double smooth_amount = 0.03;
    // int max_iter = 300;
    int iter = 0;
    double tolerance = 0.05;
    double q = 1.0;


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
        if(H == NULL) {
            std::cerr << "error in getting mean curvature" << std::endl;
            return EXIT_FAILURE;
        }
        vtkDataArray* N = normal_filter->GetOutput()->GetPointData()->GetNormals();
        if(N == NULL) {
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
        // double test_point[3];
        // points->GetPoint(10, test_point);
        // std::cout << test_point[0] << " , " << test_point[1] << " , " << test_point[2] << std::endl;

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

        if((iter +1) % freq_output == 0)
        {
            char modelName[128];
            sprintf(modelName, "output#%04d", iter+1);
            AddModelNodeToScene(mesh, modelName, false);
            vtkSmartPointer<vtkCenterOfMass> centerMassFilter =
                vtkSmartPointer<vtkCenterOfMass>::New();
            centerMassFilter->SetInputData(mesh);
            centerMassFilter->SetUseScalarsAsWeights(false);
            centerMassFilter->Update();
            double center[3];
            centerMassFilter->GetCenter(center);

            ShowFittingEllipsoid(points, curr_volume, center);

        }
        q -= 0.0001;
        iter++;
    }

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
    if(displayModelNode == NULL)
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
int vtkSlicerSkeletalRepresentationInitializerLogic::ShowFittingEllipsoid(vtkPoints* points, double curr_volume, double center[3])
{
    Eigen::MatrixXd point_matrix(points->GetNumberOfPoints(), 3);
    for(int i = 0; i < points->GetNumberOfPoints(); ++i)
    {
        double p[3];
        points->GetPoint(i, p);
        point_matrix.row(i) << p[0], p[1], p[2];
    }
    // compute best fitting ellipsoid
    // For now assume that the surface is centered and rotationally aligned
    // 1. compute the second moment
    Eigen::MatrixXd point_matrix_transposed = point_matrix.transpose();
    Eigen::Matrix3d second_moment = point_matrix_transposed * point_matrix;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(second_moment);
    Eigen::VectorXd radii = es.eigenvalues();
    radii(0) = sqrt(radii(0));
    radii(1) = sqrt(radii(1));
    radii(2) = sqrt(radii(2));

    double ellipsoid_volume = 4 / 3.0 * M_PI * radii(0) * radii(1) * radii(2);
    double volume_factor = pow(curr_volume / ellipsoid_volume, 1.0 / 3.0); 
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
    MatrixXd cog(1, 3); // center of gravity
    cog << center[0], center[1], center[2];
    MatrixXd translated_points = rotated_ellipsoid_points + cog.replicate(rotated_ellipsoid_points.rows(),1);

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

    // save ellipsoid mesh
    // std::string ellFileName("temp_ellipsoid.vtk");
    // writer->SetInputData(best_fitting_ellipsoid_polydata);
    // writer->SetFileName(ellFileName.c_str());
    // writer->Update();

    // std::string filename("temp_ellipsoid.vtk");
    // vtkSmartPointer<vtkPolyDataReader> reader =
    //     vtkSmartPointer<vtkPolyDataReader>::New();
    // reader->SetFileName(filename.c_str());
    // reader->Update();

    // vtkSmartPointer<vtkPolyData> mesh =
    //     vtkSmartPointer<vtkPolyData>::New();
    // mesh = reader->GetOutput();

    AddModelNodeToScene(best_fitting_ellipsoid_polydata, "best_fitting_ellipsoid", false, 1, 1, 0);
    return 0;
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
        if(displayNode == NULL)
        {
            continue;
        }

        displayNode->SetVisibility(0);

    }

}

int vtkSlicerSkeletalRepresentationInitializerLogic::InklingFlow(const std::string &filename, double dt, double smooth_amount, int max_iter, int freq_output, double threshold)
{
    std::cout << filename << std::endl;
    std::cout << dt << std::endl;
    std::cout << smooth_amount << std::endl;
    std::cout << max_iter << std::endl;
    std::cout << freq_output << std::endl;
    vtkSmartPointer<vtkPolyDataReader> reader =
        vtkSmartPointer<vtkPolyDataReader>::New();
    reader->SetFileName("temp_output.vtk");
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

//    while(q > tolerance && iter < max_iter)
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
        if(N == NULL) {
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


        if(H == NULL) {
            std::cerr << "error in getting mean curvature" << std::endl;
            return EXIT_FAILURE;
        }

        curvature_filter->SetCurvatureTypeToGaussian();
        curvature_filter->Update();
        vtkSmartPointer<vtkDoubleArray> K =
            vtkDoubleArray::SafeDownCast(curvature_filter->GetOutput()->GetPointData()->GetArray("Gauss_Curvature"));

        if(K == NULL) {
            std::cerr << "error in getting Gaussian curvature" << std::endl;
            return EXIT_FAILURE;
        }

        curvature_filter->SetCurvatureTypeToMaximum();
        curvature_filter->Update();
        vtkSmartPointer<vtkDoubleArray> MC =
            vtkDoubleArray::SafeDownCast(curvature_filter->GetOutput()->GetPointData()->GetArray("Maximum_Curvature"));

        if(MC == NULL) {
            std::cerr << "error in getting max curvature" << std::endl;
            return EXIT_FAILURE;
        }

        curvature_filter->SetCurvatureTypeToMinimum();
        curvature_filter->Update();
        vtkSmartPointer<vtkDoubleArray> MinC =
            vtkDoubleArray::SafeDownCast(curvature_filter->GetOutput()->GetPointData()->GetArray("Minimum_Curvature"));
        if(MinC == NULL)
        {
            std::cout << "error in getting min curvature" << std::endl;
            return -1;
        }

        // perform the flow
        vtkSmartPointer<vtkPoints> points = mesh->GetPoints();
        double maxVal = -10000.0;
        double maxIndex = -1;
        for(int i = 0; i < points->GetNumberOfPoints(); ++i) {
            double p[3];
            points->GetPoint(i, p);
            double curr_N[3];
            N->GetTuple(i, curr_N);
            double curr_H = H->GetValue(i);
            double curr_K = K->GetValue(i);
            double curr_max = MC->GetValue(i);
            double curr_min = MinC->GetValue(i);

            double delta = 0.0;
            double diffK = curr_max - curr_min;
            delta = dt * curr_H; // * 1 / (diffK * diffK);
            diffK = abs(diffK);

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
        // std::cout << "min diff^2:" << maxVal << std::endl;
        // std::cout << "min position:" << maxIndex << std::endl;
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

        // TODO: best fitting ellipsoid
        HideNodesByNameByClass("output_inkling","vtkMRMLModelNode");

        // then add this new intermediate result

        vtkSmartPointer<vtkPolyDataWriter> writer =
            vtkSmartPointer<vtkPolyDataWriter>::New();
        writer->SetInputData(mesh);
        writer->SetFileName("temp_output.vtk");
        writer->Update();

        char modelName[128];
        sprintf(modelName, "output_inkling");
        AddModelNodeToScene(mesh, modelName, true);
        // if((iter +1) % freq_output == 0)
        // {
        //     char modelName[128];
        //     sprintf(modelName, "output#%04d", iter+1);
        //     AddModelNodeToScene(mesh, modelName, false);
        // }

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
    if(fidNode == NULL)
    {
        vtkErrorMacro("fidNode is NULL");
        return;
    }
    fidNode->SetAndObserveDisplayNodeID(displayNode->GetID());
    fidNode->SetLocked(true);
    fidNode->SetName("Highlight");
    scene->AddNode(fidNode);

    fidNode->AddFiducial(x, y, z);

}