#ifndef VTKAPPROXIMATESIGNEDDISTANCEMAP_H
#define VTKAPPROXIMATESIGNEDDISTANCEMAP_H
#include <itkImage.h>
#include <string>
#include <vtkSmartPointer.h>
class vtkImageData;
class vtkApproximateSignedDistanceMap
{
public:
    typedef unsigned char Pixel;
    typedef itk::Image<Pixel, 3> ImageType;
    typedef ImageType::Pointer ImagePointer;
    
    typedef itk::Image<float, 3> RealImage;
    vtkApproximateSignedDistanceMap();
    
    void Convert(vtkImageData *input, vtkSmartPointer<vtkImageData> output);
};

#endif // VTKAPPROXIMATESIGNEDDISTANCEMAP_H
