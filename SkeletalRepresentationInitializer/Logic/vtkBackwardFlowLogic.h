// This class provides logic of backward flow
// author: Zhiyuan Liu
// Date: Sept. 4, 2018
#ifndef __vtkBackwardFlowLogic_h
#define __vtkBackwardFlowLogic_h

class vtkPolyData;
class vtkBackwardFlowLogic {
public:
    vtkBackwardFlowLogic(){}
    ~vtkBackwardFlowLogic(){}

    void runApplyTPS();
    void computePairwiseTPS(vtkPolyData* afterFlow, vtkPolyData* beforeFlow, const char* outputFileName);

};
#endif