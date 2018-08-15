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

// VTK includes
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkCurvatures.h>
#include <vtkWindowedSincPolyDataFilter.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkParametricEllipsoid.h>
#include <vtkParametricFunction.h>
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
int vtkSlicerSkeletalRepresentationInitializerLogic::FlowSurfaceMesh(const std::string &filename)
{
    // TODO: set parameters: dt, smooth_amount, max_iter
    vtkSmartPointer<vtkPolyDataReader> reader =
        vtkSmartPointer<vtkPolyDataReader>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();


    vtkSmartPointer<vtkPolyData> mesh =
        vtkSmartPointer<vtkPolyData>::New();
    mesh = reader->GetOutput();
    AddModelNodeToScene(reader->GetOutput(), "original", true);

    vtkSmartPointer<vtkMassProperties> mass_filter =
        vtkSmartPointer<vtkMassProperties>::New();
    mass_filter->SetInputData(mesh);
    mass_filter->Update();
    double original_volume = mass_filter->GetVolume();
    std::cout << "Original Volume: " << original_volume << std::endl;

    // default parameters
    double dt = 0.001;
    double smooth_amount = 100.0;
    int max_iter = 100;
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
        }
        points->Modified();
	    Eigen::MatrixXd point_matrix(points->GetNumberOfPoints(), 3);
        // compute best fitting ellipsoid
	    // For now assume that the surface is centered and rotationally aligned
	    // 1. compute the second moment
	    Eigen::MatrixXd point_matrix_transposed = point_matrix.transpose();
	    Eigen::Matrix3d second_moment = point_matrix_transposed * point_matrix;
	    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(second_moment);
        Eigen::VectorXd radii = es.eigenvalues();


	    // save the intermediate surface & output result
        // char name[128];
        // sprintf(name, "temp_%04d.vtk", iter);
        // vtkSmartPointer<vtkPolyDataWriter> writer =
        //     vtkSmartPointer<vtkPolyDataWriter>::New();
        // writer->SetInputData(mesh);
        // writer->SetFileName(name);
        // writer->Update();
        char modelName[128];
        sprintf(modelName, "output#%04d", iter);

        AddModelNodeToScene(mesh, modelName, false);

        q -= 0.0001;
        iter++;
    }

    return 1;
}

void vtkSlicerSkeletalRepresentationInitializerLogic::AddModelNodeToScene(vtkPolyData* mesh, char* modelName, bool isModelVisible)
{
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
    displayModelNode->SetColor(0, 1, 1);
    displayModelNode->SetScene(scene);
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