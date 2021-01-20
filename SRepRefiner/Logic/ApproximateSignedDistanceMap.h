#ifndef ApproximateSignedDistanceMap_H
#define ApproximateSignedDistanceMap_H
#include <itkImage.h>
#include <string>
#include "itkImageRegionIterator.h"
class vtkImageData;
template<typename TImage>
void DeepCopy(typename TImage::Pointer input, typename TImage::Pointer output)
{
  output->SetRegions(input->GetLargestPossibleRegion());
  output->Allocate();

  itk::ImageRegionConstIterator<TImage> inputIterator(input, input->GetLargestPossibleRegion());
  itk::ImageRegionIterator<TImage> outputIterator(output, output->GetLargestPossibleRegion());

  while(!inputIterator.IsAtEnd())
    {
    outputIterator.Set(inputIterator.Get());
    ++inputIterator;
    ++outputIterator;
    }
}

class ApproximateSignedDistanceMap
{
public:
    typedef unsigned char Pixel;
    typedef itk::Image<Pixel, 3> ImageType;
    typedef ImageType::Pointer ImagePointer;

    typedef itk::Image<float, 3> RealImage;
    ApproximateSignedDistanceMap();

    void Convert(vtkImageData *input, RealImage::Pointer outputITK);
};

#endif // ApproximateSignedDistanceMap_H
