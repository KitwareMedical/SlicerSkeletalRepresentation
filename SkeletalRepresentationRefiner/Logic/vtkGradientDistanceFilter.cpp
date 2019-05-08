#include "vtkGradientDistanceFilter.h"
#include "itkGradientImageFilter.h"
#include "itkImageFileWriter.h"
vtkGradientDistanceFilter::vtkGradientDistanceFilter()
{
    
}

bool vtkGradientDistanceFilter::CompareImages(VectorImage::Pointer input, VectorImage::Pointer output)
{
    VectorImage::RegionType region = input->GetLargestPossibleRegion();
    
    VectorImage::SizeType imageSize = region.GetSize();
    for(int i = 0; i < imageSize[0]; ++i)
    {
        for (int j = 0; j < imageSize[1]; ++j) {
            for (int k = 0; k < imageSize[2]; ++k) {
                VectorImage::IndexType index = {{i,j,k}};
                VectorImage::PixelType pixel = input->GetPixel(index);
                VectorImage::PixelType pixel2 = output->GetPixel(index);
                if(pixel[0] != pixel2[0] || pixel[1] != pixel2[1] || pixel[2] != pixel2[2])
                    return false;
            }
        }
    }
    return true;
}
void vtkGradientDistanceFilter::Filter(RealImage::Pointer input, VectorImage::Pointer output)
{
    typedef itk::GradientImageFilter<RealImage, float>  GradientFilterType;
    
    GradientFilterType::Pointer gradientFilter = GradientFilterType::New();
    gradientFilter->SetInput(input);
    
    try {
        gradientFilter->Update();
        VectorImage::Pointer gradImage = gradientFilter->GetOutput();
        DeepCopy(gradImage, output);
        //std::cout << "The two images are same:" << CompareImages(gradImage, output) << std::endl;

    } catch (itk::ExceptionObject & excep) {
        std::cerr << "Exception caught !" << std::endl;
        std::cerr << excep << std::endl;
    }
    
    
    
}

void vtkGradientDistanceFilter::DeepCopy(VectorImage::Pointer input, VectorImage::Pointer output)
{
    VectorImage::RegionType region = input->GetLargestPossibleRegion();
    output->SetRegions(region);
    output->Allocate();
    
    VectorImage::SizeType imageSize = region.GetSize();
    for(int i = 0; i < imageSize[0]; ++i)
    {
        for (int j = 0; j < imageSize[1]; ++j) {
            for (int k = 0; k < imageSize[2]; ++k) {
                VectorImage::IndexType index = {{i,j,k}};
                VectorImage::PixelType pixel = input->GetPixel(index);
                output->SetPixel(index, pixel);
            }
        }
    }
}
