#include "ApproximateSignedDistanceMap.h"

#include "itkVTKImageToImageFilter.h"
#include "vtkImageMagnitude.h"
#include "itkImageToVTKImageFilter.h"
#include "itkApproximateSignedDistanceMapImageFilter.h"
#include "itkImageFileWriter.h"
#include "itksys/SystemTools.hxx"
ApproximateSignedDistanceMap::ApproximateSignedDistanceMap()
{

}

void ApproximateSignedDistanceMap::Convert(vtkImageData *input, RealImage::Pointer outputITK)
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
          approximateSignedDistanceMapImageFilter->Update();

        RealImage::Pointer image = approximateSignedDistanceMapImageFilter->GetOutput();

        DeepCopy<RealImage>(image, outputITK);
        std::cout << "Finished generating distance map." << std::endl;

    }
    catch( itk::ExceptionObject & error )
    {
        std::cerr << "Error: " << error << std::endl;
        return;
    }
}

