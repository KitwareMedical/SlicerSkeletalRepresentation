/*==============================================================================

  Program: 3D Slicer

  Portions (c) Copyright Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

==============================================================================*/
#include "vtkSpoke.h"
#include <math.h>
#include <vtkMath.h>
// Eigen includes
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
vtkSpoke::vtkSpoke(){}

vtkSpoke::vtkSpoke(double radius, double px, double py, double pz, double ux, double uy, double uz)
{
    mR = radius;
    mPx = px;
    mPy = py;
    mPz = pz;

    double r = sqrt(ux * ux + uy * uy + uz * uz);
    mUx = ux / r;
    mUy = uy / r;
    mUz = uz / r;
}

vtkSpoke::vtkSpoke(double *ptSkeletal, double *ptBoundary)
{
    mPx = ptSkeletal[0];
    mPy = ptSkeletal[1];
    mPz = ptSkeletal[2];

    double diffX = ptBoundary[0] - ptSkeletal[0];
    double diffY = ptBoundary[1] - ptSkeletal[1];
    double diffZ = ptBoundary[2] - ptSkeletal[2];
    mR = sqrt(diffX * diffX + diffY * diffY + diffZ * diffZ);

    mUx = diffX / mR;
    mUy = diffY / mR;
    mUz = diffZ / mR;
}

vtkSpoke::vtkSpoke(const vtkSpoke &other)
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
    vtkMath::Normalize(u);
    mUx = u[0];
    mUy = u[1];
    mUz = u[2];
}

void vtkSpoke::GetDirection(double *output) const
{
    output[0] = this->mUx;
    output[1] = this->mUy;
    output[2] = this->mUz;

    vtkMath::Normalize(output);
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

void vtkSpoke::Diff(vtkSpoke *another, double *output) const
{
    double eX = this->mR * this->mUx;
    double eY = this->mR * this->mUy;
    double eZ = this->mR * this->mUz;

    double aX = another->mR * another->mUx;
    double aY = another->mR * another->mUy;
    double aZ = another->mR * another->mUz;

    output[0] = eX - aX;
    output[1] = eY - aY;
    output[2] = eZ - aZ;
}

void vtkSpoke::GetBoundaryPoint(double *output) const
{
    output[0] = mPx + mR * mUx;
    output[1] = mPy + mR * mUy;
    output[2] = mPz + mR * mUz;
}

void vtkSpoke::SetNeighborU(const std::vector<vtkSpoke *> &neighbors, bool isForward)
{
    mNeighborsU.clear();
    for (size_t i = 0; i < neighbors.size(); ++i) {
        mNeighborsU.push_back(neighbors[i]);
    }
    mIsForwardU = isForward;
}

void vtkSpoke::SetNeighborV(const std::vector<vtkSpoke *> &neighbors, bool isForward)
{
    mNeighborsV.clear();
    for (size_t i = 0; i < neighbors.size(); ++i) {
        mNeighborsV.push_back(neighbors[i]);
    }
    mIsForwardV = isForward;
}

double vtkSpoke::GetRSradPenalty(double stepSize)
{
    // 1. compute derivatives
    double dxdu[3], dxdv[3], dSdu[3], dSdv[3], U[3], drdu, drdv; // variables that will be used

    // neighbors in each direction could be either 1 or 2
    ComputeDerivatives(mNeighborsU, mIsForwardU, stepSize, dxdu, dSdu, &drdu);
    ComputeDerivatives(mNeighborsV, mIsForwardV, stepSize, dxdv, dSdv, &drdv);
    this->GetDirection(U);

    // 2. construct rSrad Matrix
    double UTU[3][3]; // UT*U - I
    UTU[0][0] = U[0] * U[0] - 1;
    UTU[0][1] = U[0] * U[1];
    UTU[0][2] = U[0] * U[2];
    UTU[1][0] = U[1] * U[0];
    UTU[1][1] = U[1] * U[1] -1;
    UTU[1][2] = U[1] * U[2];
    UTU[2][0] = U[2] * U[0];
    UTU[2][1] = U[2] * U[1];
    UTU[2][2] = U[2] * U[2] -1;

    // Notation in Han, Qiong's dissertation
    Eigen::MatrixXd Q(2,3);
    Q(0,0) = dxdu[0] * UTU[0][0] + dxdu[1] * UTU[1][0] + dxdu[2] * UTU[2][0];
    Q(0,1) = dxdu[0] * UTU[0][1] + dxdu[1] * UTU[1][1] + dxdu[2] * UTU[2][1];
    Q(0,2) = dxdu[0] * UTU[0][2] + dxdu[1] * UTU[1][2] + dxdu[2] * UTU[2][2];

    Q(1,0) = dxdv[0] * UTU[0][0] + dxdv[1] * UTU[1][0] + dxdv[2] * UTU[2][0];
    Q(1,1) = dxdv[0] * UTU[0][1] + dxdv[1] * UTU[1][1] + dxdv[2] * UTU[2][1];
    Q(1,2) = dxdv[0] * UTU[0][2] + dxdv[1] * UTU[1][2] + dxdv[2] * UTU[2][2];

    Eigen::MatrixXd leftSide(2,3), rightSide(3, 2);
    leftSide(0,0) = dSdu[0] - drdu * U[0];
    leftSide(0,1) = dSdu[1] - drdu * U[1];
    leftSide(0,2) = dSdu[2] - drdu * U[2];

    leftSide(1,0) = dSdv[0] - drdv * U[0];
    leftSide(1,1) = dSdv[1] - drdv * U[1];
    leftSide(1,2) = dSdv[2] - drdv * U[2];

    Eigen::Matrix2d QQT, QQT_inv;
    QQT = Q * Q.transpose();
    QQT_inv = QQT.inverse();

    rightSide = Q.transpose() * QQT_inv;

    Eigen::Matrix2d rSradMat;
    rSradMat = leftSide * rightSide;
    rSradMat.transposeInPlace();
    // 3. compute rSrad penalty
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigensolver(rSradMat);
    double maxEigen = eigensolver.eigenvalues()[1];

    if(maxEigen > 1) return maxEigen-1;
    else return 0.0;
//    double detRSrad = rSradMat.determinant();
//    if(detRSrad < 0) return 0.0;
//    else if(detRSrad < 1) return 0.0;
//    else return (abs(detRSrad)-1);
}

bool vtkSpoke::IsValid()
{
    return !(std::isnan(mR) ||
            std::isnan(mUx) || std::isnan(mUy) || std::isnan(mUz)/*
            || std::isnan(mPx) || std::isnan(mPy) || std::isnan(mPz)*/);
}

void vtkSpoke::ComputeDerivatives(std::vector<vtkSpoke*> neibors, bool isForward, double stepSize, // input
                                  double *dxdu, double *dSdu, double *drdu) // output
{
    if(neibors.size() == 1)
    {
        double neiborX[3], neiborR;
        neibors[0]->GetSkeletalPoint(neiborX);
        neiborR = neibors[0]->GetRadius();
        if(isForward)
        {
            dxdu[0] = neiborX[0] - mPx;
            dxdu[1] = neiborX[1] - mPy;
            dxdu[2] = neiborX[2] - mPz;
            *drdu = neiborR - mR;
            neibors[0]->Diff(this, dSdu);
        }
        else {
            dxdu[0] = mPx - neiborX[0];
            dxdu[1] = mPy - neiborX[1];
            dxdu[2] = mPz - neiborX[2];
            *drdu = mR - neiborR;
            this->Diff(neibors[0], dSdu);
        }

        dxdu[0] /= stepSize;
        dxdu[1] /= stepSize;
        dxdu[2] /= stepSize;
        *drdu /= stepSize;

        dSdu[0]  /= stepSize;
        dSdu[1]  /= stepSize;
        dSdu[2]  /= stepSize;
    }
    else if(neibors.size() == 2)
    {
        double neiborX1[3], neiborX0[3], neiborR1, neiborR0;
        neibors[1]->Diff(neibors[0], dSdu);
        neibors[1]->GetSkeletalPoint(neiborX1);
        neibors[0]->GetSkeletalPoint(neiborX0);

        neiborR1 = neibors[1]->GetRadius();
        neiborR0 = neibors[0]->GetRadius();

        *drdu = (neiborR1 - neiborR0) / stepSize / 2;
        dxdu[0] = (neiborX1[0] - neiborX0[0]) / stepSize / 2;
        dxdu[1] = (neiborX1[1] - neiborX0[1]) / stepSize / 2;
        dxdu[2] = (neiborX1[2] - neiborX0[2]) / stepSize / 2;
    }
}
