#ifndef VTKSREP_H
#define VTKSREP_H
#include <vector>

class vtkSpoke;
class vtkSrep
{
public:
    vtkSrep();
    
    // new points
    vtkSrep(int r, int c, std::vector<double> &radii, std::vector<double> &dirs, std::vector<double> &skeletalPoints);
    
    // delete points
    ~vtkSrep(); 
    
    // look for spoke at row = r, col = c
    vtkSpoke *GetSpoke(int r, int c) const;
    
    // False if there is no spoke in this srep
    bool IsEmpty() const;
    
    // Get the vector of spokes
    std::vector<vtkSpoke *> &GetAllSpokes();
    
private:
    int nRows;
    int nCols;
    std::vector<vtkSpoke*> spokes;
};

#endif // VTKSREP_H
