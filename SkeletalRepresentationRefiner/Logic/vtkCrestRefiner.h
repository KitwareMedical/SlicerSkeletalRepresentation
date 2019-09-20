#ifndef VTKCRESTREFINER_H
#define VTKCRESTREFINER_H
#include <vector>
#include "vtkSlicerSkeletalRepresentationRefinerLogic.h"
class vtkSpoke;

class vtkCrestRefiner : public vtkSlicerSkeletalRepresentationRefinerLogic
{
public:
    vtkCrestRefiner();
    // Description: Override operator (). Required by min_newuoa.
    // Parameter: @coeff: the pointer to coefficients
    double operator () (double *coeff);

    // This function returns the cost value result from defined objective function
    // given the current coeff array
    double EvaluateObjectiveFunction(double *coeff);

    void ChangeSpokes(double *coeff, std::vector<vtkSpoke*> &orig, std::vector<vtkSpoke*> &newSpokes);

    void SetSpokes(std::vector<vtkSpoke*> &spokes);

    void SetInterpolationLevel(int interpLevel);

    void SetInteriorSpokes(std::vector<vtkSpoke*> &upSpokes, std::vector<vtkSpoke*> &downSpokes);

    void SetRowsCols(int r, int c);

    void SetTargetMeshFilePath(const std::string &meshFile);

    void SetWeights(double imageMatch, double normalMatch, double srad);

    void SetTransformationMat(double mat[][4]);

    // optimize crest spokes directions to be the dir pointing to upS + downS
    void AdjustAngles();

private:
    void CommensuratePenalties(double* imageMatch, double* normalMatch, double* srad);
private:
    RealImage::Pointer mAntiAliasedImage = RealImage::New();
    VectorImage::Pointer mGradDistImage = VectorImage::New();
    std::string mTargetMeshFilePath;
    std::vector<vtkSpoke*> mSpokes;
    int mInterpolationLevel;
    std::vector<vtkSpoke*> mUpSpokes, mDownSpokes;
    int mNumRows, mNumCols;
    // when apply this transformation: [x, y, z, 1] * mTransformationMat
    double mTransformationMat[4][4]; // homogeneous matrix transfrom from srep to unit cube cs.
    // weights in optimization algorithm
    double mWtImageMatch;
    double mWtNormalMatch;
    double mWtSrad;
};

#endif // VTKCRESTREFINER_H
