#ifndef VTKANTIALIAS_H
#define VTKANTIALIAS_H

#include "vtkSmartPointer.h"

class vtkImageData;
class vtkAntiAlias
{
public:
    vtkAntiAlias();
    
    void Filter(vtkImageData *input, vtkSmartPointer<vtkImageData> output);
};

#endif // VTKANTIALIAS_H
