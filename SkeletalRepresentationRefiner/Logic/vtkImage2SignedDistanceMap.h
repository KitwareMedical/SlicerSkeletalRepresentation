#ifndef VTKIMAGE2SIGNEDDISTANCEMAP_H
#define VTKIMAGE2SIGNEDDISTANCEMAP_H
#include <itkImage.h>
#include "vtkSmartPointer.h"

class vtkImageData;
class vtkImage2SignedDistanceMap
{
public:
    typedef unsigned char Pixel;
    typedef itk::Image<Pixel, 3> ImageType;
    typedef ImageType::Pointer ImagePointer;
    
    typedef itk::Image<float, 3> RealImage;
    
    vtkImage2SignedDistanceMap();
    
    void Convert(vtkImageData *input, vtkSmartPointer<vtkImageData> output);
    
private:
    ImagePointer CloseImage(ImagePointer image, double radius);
};

#endif // VTKIMAGE2SIGNEDDISTANCEMAP_H
