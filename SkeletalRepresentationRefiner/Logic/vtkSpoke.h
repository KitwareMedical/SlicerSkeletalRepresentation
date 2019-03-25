#ifndef VTKSPOKE_H
#define VTKSPOKE_H

class vtkSpoke
{
public:
    vtkSpoke();
    vtkSpoke(double radius, double px, double py, double pz, double ux, double uy, double uz);
    // deep copy in assignment 
    vtkSpoke& operator=(const vtkSpoke& other);
    
    void SetRadius(double r);
    
    void SetSkeletalPoint(double px, double py, double pz);
    
    void SetDirection(double *u);
    
    void GetDirection(double *output) const;
    
    double GetRadius() const;
    
    void GetSkeletalPoint(double *output) const;
    
    void Add(vtkSpoke* another, double* output) const;
    
    void GetBoundaryPoint(double *output) const;
private:
    double mR;
    double mPx;
    double mPy;
    double mPz;
    double mUx;
    double mUy;
    double mUz;
};

#endif // VTKSPOKE_H

