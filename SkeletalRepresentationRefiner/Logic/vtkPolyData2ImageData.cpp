/*==============================================================================

  This class creates an image of the poly data of surface mesh.
  The image is mapped to the unit cube centered at (0.5, 0.5, 0.5).
  By doing so, the image can be both shown in Pablo and Slicer.

==============================================================================*/
#include "vtkPolyData2ImageData.h"
#include <vtkVersion.h>
#include <vtkPolyData.h>
#include <vtkImageData.h>
#include <vtkSphereSource.h>
#include <vtkMetaImageWriter.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkImageStencil.h>
#include <vtkPointData.h>
#include <vtkGenericDataObjectReader.h>
#include <vtkPolyDataWriter.h>
#include <math.h>

vtkPolyData2ImageData::vtkPolyData2ImageData()
{

}

void vtkPolyData2ImageData::Convert(const std::string &inputFileName, vtkSmartPointer<vtkImageData> output)
{
    // Get all data from the file
    vtkSmartPointer<vtkGenericDataObjectReader> reader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
    reader->SetFileName(inputFileName.c_str());
    reader->Update();

    vtkPolyData* inputData = reader->GetPolyDataOutput();
    double bounds[6], range[3];

    // 1. transform the mesh into unit cube
    vtkSmartPointer<vtkPolyData> transMesh = vtkSmartPointer<vtkPolyData>::New();
    inputData->GetBounds(bounds);
    double newCenter[3] = {0.5, 0.5, 0.5};
    range[0] = bounds[1] - bounds[0]; // The range of x coordinate
    range[1] = bounds[3] - bounds[2];
    range[2] = bounds[5] - bounds[4];
    double ratioYX, ratioZX, ratioZY;
    ratioYX = range[1] / range[0];
    ratioZX = range[2] / range[0];
    ratioZY = range[2] / range[1];

    double newBounds[6] = {0.};
    // put the longest axis to [0,1], scale other coordinates accordingly
    if(range[0] >= range[1] && range[0] >= range[2])
    {
        newBounds[0] = 0.0;
        newBounds[1] = 1.0;
        newBounds[2] = newCenter[1] - 0.5 * ratioYX;
        newBounds[3] = newCenter[1] + 0.5 * ratioYX;
        newBounds[4] = newCenter[2] - 0.5 * ratioZX;
        newBounds[5] = newCenter[2] + 0.5 * ratioZX;
    }
    else if(range[1] >= range[0] && range[1] >= range[2])
    {
        newBounds[0] = newCenter[0] - 0.5 / ratioYX;
        newBounds[1] = newCenter[0] + 0.5 / ratioYX;
        newBounds[2] = 0.0;
        newBounds[3] = 1.0;
        newBounds[4] = newCenter[2] - 0.5 * ratioZY;
        newBounds[5] = newCenter[2] + 0.5 * ratioZY;
    }
    else if(range[2] >= range[0] && range[2] >= range[1])
    {
        newBounds[0] = newCenter[0] - 0.5 / ratioZX;
        newBounds[1] = newCenter[0] + 0.5 / ratioZX;
        newBounds[2] = newCenter[1] - 0.5 / ratioZY;
        newBounds[3] = newCenter[1] + 0.5 / ratioZY;
        newBounds[4] = 0.0;
        newBounds[5] = 1.0;
    }

    double rangeTransMesh[3];
    rangeTransMesh[0] = newBounds[1] - newBounds[0];
    rangeTransMesh[1] = newBounds[3] - newBounds[2];
    rangeTransMesh[2] = newBounds[5] - newBounds[4];
    vtkSmartPointer<vtkPoints> newPts = vtkSmartPointer<vtkPoints>::New();
    for(int i = 0; i < inputData->GetNumberOfPoints(); ++i)
    {
        double oldPt[3], newPt[3];
        inputData->GetPoint(i, oldPt);
        newPt[0] = rangeTransMesh[0] * (oldPt[0] - bounds[0]) / range[0] + newBounds[0];
        newPt[1] = rangeTransMesh[1] * (oldPt[1] - bounds[2]) / range[1] + newBounds[2];
        newPt[2] = rangeTransMesh[2] * (oldPt[2] - bounds[4]) / range[2] + newBounds[4];
        newPts->InsertPoint(i, newPt);
    }
    newPts->Modified();
    transMesh->SetPoints(newPts);
    transMesh->SetPolys(inputData->GetPolys());

    vtkSmartPointer<vtkImageData> whiteImage = vtkSmartPointer<vtkImageData>::New();

    double spacing[3]; // desired volume spacing
    double voxel_spacing = 0.005;
    spacing[0] = voxel_spacing;
    spacing[1] = voxel_spacing;
    spacing[2] = voxel_spacing;
    whiteImage->SetSpacing(spacing);

    // compute dimensions
    int dim[3];
    for (int i = 0; i < 3; i++)
    {
        dim[i] = static_cast<int>(1 / spacing[i]);
    }
    whiteImage->SetDimensions(dim);
    whiteImage->SetExtent(0, dim[0] - 1, 0, dim[1] - 1, 0, dim[2] - 1);

    double origin[3];
    origin[0] = newBounds[0]; //bounds[0] + spacing[0]/2;//0.5 * (bounds[0] + bounds[1]);
    origin[1] = newBounds[2]; //bounds[2] + spacing[1]/2;//0.5 * (bounds[3] + bounds[2]);
    origin[2] = newBounds[4]; //bounds[4] + spacing[2]/2;//0.5 * (bounds[5] + bounds[4]);

    whiteImage->SetOrigin(origin);

#if VTK_MAJOR_VERSION <= 5
    whiteImage->SetScalarTypeToUnsignedChar();
    whiteImage->AllocateScalars();
#else
    whiteImage->AllocateScalars(VTK_UNSIGNED_CHAR,1);
#endif
    // fill the image with foreground voxels:
    unsigned char inval = 255;
    unsigned char outval = 0;
    vtkIdType count = whiteImage->GetNumberOfPoints();
    for (vtkIdType i = 0; i < count; ++i)
    {
        whiteImage->GetPointData()->GetScalars()->SetTuple1(i, inval);
    }

    // polygonal data --> image stencil:
    vtkSmartPointer<vtkPolyDataToImageStencil> pol2stenc = vtkSmartPointer<vtkPolyDataToImageStencil>::New();
#if VTK_MAJOR_VERSION <= 5
    pol2stenc->SetInput(transMesh);
#else
    pol2stenc->SetInputData(transMesh);
#endif
    pol2stenc->SetTolerance(0);
    double polOrigin[3] = {0.0, 0.0, 0.0};
    pol2stenc->SetOutputOrigin(polOrigin);
    pol2stenc->SetOutputSpacing(spacing);
    pol2stenc->SetOutputWholeExtent(whiteImage->GetExtent());
    pol2stenc->Update();

    // cut the corresponding white image and set the background:
    vtkSmartPointer<vtkImageStencil> imgstenc =
        vtkSmartPointer<vtkImageStencil>::New();
#if VTK_MAJOR_VERSION <= 5
    imgstenc->SetInput(whiteImage);
    imgstenc->SetStencil(pol2stenc->GetOutput());
#else
    imgstenc->SetInputData(whiteImage);
    imgstenc->SetStencilConnection(pol2stenc->GetOutputPort());
#endif
    imgstenc->ReverseStencilOff();
    imgstenc->SetBackgroundValue(outval);
    imgstenc->Update();

    output->DeepCopy(imgstenc->GetOutput());
}
