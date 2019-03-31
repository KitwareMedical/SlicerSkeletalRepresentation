#ifndef VTKGRADIENTDISTANCEFILTER_H
#define VTKGRADIENTDISTANCEFILTER_H

#include "itkImage.h"
#include "itkCovariantVector.h"
// This class compute gradient of signed distance map, resulting in normals of the image everywhere
class vtkGradientDistanceFilter
{
public:
    typedef itk::Image<float, 3> RealImage;
    typedef itk::Image<itk::CovariantVector<float, 3>, 3> VectorImage;
    
    vtkGradientDistanceFilter();
    
    void Filter(RealImage::Pointer input, VectorImage::Pointer output);
    
private:
    void DeepCopy(VectorImage::Pointer input, VectorImage::Pointer output);
    bool CompareImages(VectorImage::Pointer input, VectorImage::Pointer output);
};

#endif // VTKGRADIENTDISTANCEFILTER_H
