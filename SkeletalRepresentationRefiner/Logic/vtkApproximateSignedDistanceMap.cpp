#include "vtkApproximateSignedDistanceMap.h"

#include "itkVTKImageToImageFilter.h"
#include "vtkImageMagnitude.h"
#include "itkImageToVTKImageFilter.h"
#include "itkApproximateSignedDistanceMapImageFilter.h"
#include "itkImageFileWriter.h"
#include "itksys/SystemTools.hxx"
vtkApproximateSignedDistanceMap::vtkApproximateSignedDistanceMap()
{
    
}

void vtkApproximateSignedDistanceMap::Convert(vtkImageData *input, vtkSmartPointer<vtkImageData> output)
{
    std::cout << "Started to convert image to signed distance map..." << std::endl;
    using FilterType = itk::VTKImageToImageFilter< ImageType >;
   
    vtkSmartPointer< vtkImageMagnitude > magnitude =
       vtkSmartPointer< vtkImageMagnitude >::New();
    magnitude->SetInputData(input );
    magnitude->Update();
    FilterType::Pointer filter = FilterType::New();
    filter->SetInput(magnitude->GetOutput());
    
    try
    {
        filter->Update();
        typedef  itk::ApproximateSignedDistanceMapImageFilter< ImageType, RealImage  > ApproximateSignedDistanceMapImageFilterType;
          ApproximateSignedDistanceMapImageFilterType::Pointer approximateSignedDistanceMapImageFilter =
            ApproximateSignedDistanceMapImageFilterType::New();
          approximateSignedDistanceMapImageFilter->SetInput(filter->GetOutput());
          approximateSignedDistanceMapImageFilter->SetInsideValue(255);
          approximateSignedDistanceMapImageFilter->SetOutsideValue(0);
        
        typedef itk::ImageFileWriter< RealImage > WriterType;
        WriterType::Pointer writer = WriterType::New();
        writer->SetFileName("/playpen/workspace/newuoa/testsdm.mhd");
        writer->SetInput(approximateSignedDistanceMapImageFilter->GetOutput());
        writer->Update();
        // Convert itk image back to vtk image data
        using BackFilterType = itk::ImageToVTKImageFilter< RealImage >;
        BackFilterType::Pointer backFilter = BackFilterType::New();
        backFilter->SetInput(approximateSignedDistanceMapImageFilter->GetOutput());
        backFilter->Update();
        
        output->DeepCopy(backFilter->GetOutput());
        std::cout << "Finished generating distance map." << std::endl;
        
    }
    catch( itk::ExceptionObject & error )
    {
        std::cerr << "Error: " << error << std::endl;
        return;
    }
}

