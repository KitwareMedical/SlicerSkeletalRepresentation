// This implementation is based on conventional formulation
// Author: Zhiyuan Liu
// Date: Sept, 2018
#include "vtkForwardFlowLogic.h"
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


// flow surface to the end: either it's ellipsoidal enough or reach max_iter
int vtkForwardFlowLogic::FlowSurfaceMesh(const std::string &filename, double dt, double smooth_amount, int max_iter, int freq_output)
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

// flow surface in one step
// basic idea: When the user select a mesh file, make a copy of vtk file in the application path.
// In each step of flow, read in that copy, flow it and save it the same place with same name.
// TODO: cleanup the hard disk when the module exits
int vtkForwardFlowLogic::FlowSurfaceOneStep(double dt, double smooth_amount)
{
    std::cout << "flow one step : dt-" << dt << std::endl;
    std::cout << "flow one step : smooth amount-" << smooth_amount << std::endl;
    std::string name = "temp_output.vtk";
    vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
    reader->SetFileName(name.c_str());
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

    vtkSmartPointer<vtkDoubleArray> H =
        vtkDoubleArray::SafeDownCast(curvature_filter->GetOutput()->GetPointData()->GetArray("Mean_Curvature"));
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
void vtkForwardFlowLogic::AddModelNodeToScene(vtkPolyData* mesh, const char* modelName, bool isModelVisible, double r, double g, double b)
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
