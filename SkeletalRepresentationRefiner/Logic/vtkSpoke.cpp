#include "vtkSpoke.h"
vtkSpoke::vtkSpoke(){}

vtkSpoke::vtkSpoke(double radius, double px, double py, double pz, double ux, double uy, double uz)
{
    mR = radius;
    mPx = px;
    mPy = py;
    mPz = pz;
    mUx = ux;
    mUy = uy;
    mUz = uz;
}

vtkSpoke &vtkSpoke::operator=(const vtkSpoke &other)
{
    double u[3], p[3], r;
    other.GetDirection(u);
    other.GetSkeletalPoint(p);
    r = other.GetRadius();
    
    mR = r;
    mPx = p[0];
    mPy = p[1];
    mPz = p[2];
    
    mUx = u[0];
    mUy = u[1];
    mUz = u[2];
    return *this;
}

void vtkSpoke::SetRadius(double r)
{
    mR = r;
}

void vtkSpoke::SetSkeletalPoint(double px, double py, double pz)
{
    mPx = px;
    mPy = py;
    mPz = pz;
}

void vtkSpoke::SetDirection(double *u)
{
    mUx = u[0];
    mUy = u[1];
    mUz = u[2];
}

void vtkSpoke::GetDirection(double *output) const
{
    output[0] = this->mUx;
    output[1] = this->mUy;
    output[2] = this->mUz;
}

double vtkSpoke::GetRadius() const
{
    return mR;
}

void vtkSpoke::GetSkeletalPoint(double *output) const
{
    output[0] = mPx;
    output[1] = mPy;
    output[2] = mPz;
}

void vtkSpoke::Add(vtkSpoke *another, double *output) const
{
    double eX = this->mR * this->mUx;
    double eY = this->mR * this->mUy;
    double eZ = this->mR * this->mUz;
    
    double aX = another->mR * another->mUx;
    double aY = another->mR * another->mUy;
    double aZ = another->mR * another->mUz;
    
    output[0] = eX + aX;
    output[1] = eY + aY;
    output[2] = eZ + aZ;
    
}

void vtkSpoke::GetBoundaryPoint(double *output) const
{
    output[0] = mPx + mR * mUx;
    output[1] = mPy + mR * mUy;
    output[2] = mPz + mR * mUz;
}
