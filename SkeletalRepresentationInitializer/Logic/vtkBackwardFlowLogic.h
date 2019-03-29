// This class provides logic of backward flow
// author: Zhiyuan Liu
// Date: Sept. 4, 2018
#ifndef __vtkBackwardFlowLogic_h
#define __vtkBackwardFlowLogic_h

#include "vtkObjectBase.h"

class vtkPolyData;
class vtkBackwardFlowLogic : public vtkObjectBase{
public:
    vtkBackwardFlowLogic(){}
    ~vtkBackwardFlowLogic(){}

    void runApplyTPS();
    void computePairwiseTPS(vtkPolyData* afterFlow, vtkPolyData* beforeFlow, const char* outputFileName);
    // void generateEllipsoidSrep(int numRow, int numCol, double ra, double rb, double rc, const char* outputPath);
};
#endif
