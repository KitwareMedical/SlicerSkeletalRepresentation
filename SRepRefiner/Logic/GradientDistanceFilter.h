#ifndef GradientDistanceFilter_H
#define GradientDistanceFilter_H

#include "itkImage.h"
#include "itkCovariantVector.h"
// This class compute gradient of signed distance map, resulting in normals of the image everywhere
class GradientDistanceFilter
{
public:
    typedef itk::Image<float, 3> RealImage;
    typedef itk::Image<itk::CovariantVector<float, 3>, 3> VectorImage;

    GradientDistanceFilter();

    void Filter(RealImage::Pointer input, VectorImage::Pointer output);

private:
    void DeepCopy(VectorImage::Pointer input, VectorImage::Pointer output);
    bool CompareImages(VectorImage::Pointer input, VectorImage::Pointer output);
};

#endif // GradientDistanceFilter_H
