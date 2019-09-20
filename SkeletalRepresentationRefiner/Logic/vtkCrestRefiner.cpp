#include "vtkCrestRefiner.h"
#include <iostream>
#include <vtkParametricSpline.h>
#include <vtkParametricFunctionSource.h>
#include "vtkSpoke.h"
vtkCrestRefiner::vtkCrestRefiner()
{
}

double vtkCrestRefiner::operator ()(double *coeff)
{
    double cost = 0.0;
    cost = EvaluateObjectiveFunction(coeff);
    return cost;
}

double vtkCrestRefiner::EvaluateObjectiveFunction(double *coeff)
{
    // this temporary srep is constructed to compute the cost function value
    // The original srep should not be changed by each iteration
//    for (int i = 0; i < mSpokes.size(); ++i) {
//        std::cout << "coeff[" << i << "]=" << coeff[i] << std::endl;
//    }
    std::vector<vtkSpoke*> newSpokes;
    ChangeSpokes(coeff, mSpokes, newSpokes);
    double imageDist = 0.0, normal = 0.0, srad = 0.0;
    // interpolate all crest spokes according to interpolation level
//    std::vector<vtkSpoke *> upInterpSpokes, downInterpSpokes, crestInterpSpokes, tempInterp;

//    InterpolateCrest(newSpokes, mUpSpokes, mInterpolationLevel, mNumRows, mNumCols, crestInterpSpokes, upInterpSpokes);

//    InterpolateCrest(newSpokes, mDownSpokes, mInterpolationLevel, mNumRows, mNumCols, tempInterp, downInterpSpokes);
//    double ptCrest[3], ptUp[3], ptDown[3], ptInterp[3], ptSkeletal[3], du[9];
    // shares from up spoke to resp. down spoke
//    int shares = 2 * static_cast<int>(pow(2, mInterpolationLevel));
//    double interval = static_cast<double>((1.0/ shares));

    for (size_t i = 0; i < newSpokes.size(); ++i) {
        vtkSpoke *thisSpoke = newSpokes[i];

        // compute distance for this spoke
        imageDist += ComputeDistance(thisSpoke, &normal);

        // interpolate from up to down and compute each distance to boundary
//        vtkSmartPointer<vtkPoints> radialCurve = vtkSmartPointer<vtkPoints>::New();
//        upInterpSpokes[i]->GetBoundaryPoint(ptUp);
//        downInterpSpokes[i]->GetBoundaryPoint(ptDown);
//        crestInterpSpokes[i]->GetBoundaryPoint(ptCrest);

//        radialCurve->InsertNextPoint(ptUp);
//        radialCurve->InsertNextPoint(ptCrest);
//        radialCurve->InsertNextPoint(ptDown);
//        vtkSmartPointer<vtkParametricSpline> splineRadial =
//                vtkSmartPointer<vtkParametricSpline>::New();
//        splineRadial->SetPoints(radialCurve);
//        vtkSmartPointer<vtkParametricFunctionSource> functionSourceRadial =
//                vtkSmartPointer<vtkParametricFunctionSource>::New();
//        functionSourceRadial->SetParametricFunction(splineRadial);
//        functionSourceRadial->Update();
//        // share the base point among all other interpolated spokes
//        crestInterpSpokes[i]->GetSkeletalPoint(ptSkeletal);
//        for(int j = 0; j <= shares; ++j)
//        {
//            double uInterp = j * interval;
//            double u[3] = {uInterp, uInterp, uInterp};
//            splineRadial->Evaluate(u, ptInterp, du);
//            vtkSpoke* interpSpoke = new vtkSpoke(ptSkeletal, ptInterp);
//            imageDist += ComputeDistance(interpSpoke, &normal);
//        }
    }
    std::cout << "image match weight is:" << mWtImageMatch << " while the total lost is:" << mWtImageMatch * imageDist << std::endl;
    for (int i = 0; i < newSpokes.size(); ++i) {
        delete newSpokes[i];
        newSpokes[i] = nullptr;
    }
    return mWtImageMatch * imageDist/* + mWtNormalMatch * normal + mWtSrad * srad*/;
}

void vtkCrestRefiner::ChangeSpokes(double *coeff, std::vector<vtkSpoke *> &orig, std::vector<vtkSpoke *> &newSpokes)
{
    if(orig.empty())
    {
        return;
    }
    for(size_t i = 0; i < orig.size(); ++i)
    {
//        size_t idx = i * 4;
//        double newU[3], newR, oldR;
//        newU[0] = coeff[idx];
//        newU[1] = coeff[idx+1];
//        newU[2] = coeff[idx+2];

        vtkSpoke* thisSpoke = orig[i];
        double oldR = thisSpoke->GetRadius();
        double newR = exp(coeff[i]) * oldR;
        double oldU[3];
        thisSpoke->GetDirection(oldU);

        double ptSkeletal[3];
        thisSpoke->GetSkeletalPoint(ptSkeletal);
        vtkSpoke* newSpoke = new vtkSpoke(newR, ptSkeletal[0], ptSkeletal[1], ptSkeletal[2], oldU[0], oldU[1], oldU[2]);
        newSpokes.push_back(newSpoke);
    }
}

void vtkCrestRefiner::SetSpokes(std::vector<vtkSpoke *> &spokes)
{
    mSpokes = spokes;
}

void vtkCrestRefiner::SetInterpolationLevel(int interpLevel)
{
    mInterpolationLevel = interpLevel;
}

void vtkCrestRefiner::SetInteriorSpokes(std::vector<vtkSpoke *> &upSpokes,
                                        std::vector<vtkSpoke *> &downSpokes)
{
    mUpSpokes = upSpokes;
    mDownSpokes = downSpokes;
}

void vtkCrestRefiner::SetRowsCols(int r, int c)
{
    mNumCols = c;
    mNumRows = r;
}

void vtkCrestRefiner::SetTargetMeshFilePath(const std::string &meshFile)
{
    mTargetMeshFilePath = meshFile;
    // Prepare signed distance image
    AntiAliasSignedDistanceMap(mTargetMeshFilePath);

}

void vtkCrestRefiner::SetWeights(double imageMatch, double normalMatch, double srad)
{
    mWtImageMatch = imageMatch;
    mWtNormalMatch = normalMatch;
    mWtSrad = srad;
}

void vtkCrestRefiner::SetTransformationMat(double mat[][4])
{
    mTransformationMat[0][0] = mat[0][0];
    mTransformationMat[0][1] = mat[0][1];
    mTransformationMat[0][2] = mat[0][2];
    mTransformationMat[0][3] = mat[0][3];

    mTransformationMat[1][0] = mat[1][0];
    mTransformationMat[1][1] = mat[1][1];
    mTransformationMat[1][2] = mat[1][2];
    mTransformationMat[1][3] = mat[1][3];

    mTransformationMat[2][0] = mat[2][0];
    mTransformationMat[2][1] = mat[2][1];
    mTransformationMat[2][2] = mat[2][2];
    mTransformationMat[2][3] = mat[2][3];

    mTransformationMat[3][0] = mat[3][0];
    mTransformationMat[3][1] = mat[3][1];
    mTransformationMat[3][2] = mat[3][2];
    mTransformationMat[3][3] = mat[3][3];
}

void vtkCrestRefiner::AdjustAngles()
{
    int idCrest = 0;
    for(int i = 0; i < mUpSpokes.size(); ++i)
    {
        int r = i / mNumCols;
        int c = i - (r * mNumCols);
        if(r == 0 || r == mNumRows - 1 || c == 0 || c == mNumCols - 1) {
            // update related crest spokes
            double upU[3], downU[3];
            mUpSpokes[i]->GetDirection(upU);
            mDownSpokes[i]->GetDirection(downU);
            double upR, downR;
            upR = mUpSpokes[i]->GetRadius();
            downR = mDownSpokes[i]->GetRadius();
            // crest spoke dir = upS + downS
            double crestU[3];
            crestU[0] = upR * upU[0] + downR * downU[0];
            crestU[1] = upR * upU[1] + downR * downU[1];
            crestU[2] = upR * upU[2] + downR * downU[2];
            vtkMath::Normalize(crestU);
            mSpokes[idCrest]->SetDirection(crestU);
            idCrest++;
        }
    }
}

void vtkCrestRefiner::CommensuratePenalties(double* imageMatch, double* normalMatch, double* srad)
{
    // map them between 10 and 100
    while(*imageMatch < 10)
    {
        *imageMatch *= 10;
    }
    while(*imageMatch > 100) {
        *imageMatch /= 10;
    }
    while(*normalMatch < 10)
    {
        *normalMatch *= 10;
    }
    while(*normalMatch > 100) {
        *normalMatch /= 10;
    }
    while(*srad < 10)
    {
        *srad *= 10;
    }
    while(*srad > 100) {
        *srad /= 10;
    }
}

