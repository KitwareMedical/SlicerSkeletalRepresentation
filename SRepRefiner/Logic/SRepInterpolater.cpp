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
#include "SRepInterpolater.h"
#include "Spoke.h"

// STD includes
#include <cmath>

SRepInterpolater::SRepInterpolater()
{

}

void SRepInterpolater::Interpolate(double u, double v,
                                                              Spoke **cornerSpokes, Spoke *interpolatedSpoke)
{
    // 1. interpolate spoke length and direction
    double lambda = 1.0;
    InterpolateQuad(cornerSpokes, u, v, lambda, interpolatedSpoke);

    double dir[3];
    interpolatedSpoke->GetDirection(dir);

    // 2. interpolate skeletal point
    double pt[3];
    InterpolateSkeletalPoint(cornerSpokes, u, v, pt);
    interpolatedSpoke->SetSkeletalPoint(pt[0], pt[1], pt[2]);
}

void SRepInterpolater::InterpolateQuad(Spoke **cornerSpokes,
                                                                  double u, double v,
                                                                  double lambda, Spoke *interpolatedSpoke)
{
    Spoke* Sp11 = cornerSpokes[0];
    Spoke* Sp21 = cornerSpokes[1];
    Spoke* Sp22 = cornerSpokes[2];
    Spoke* Sp12 = cornerSpokes[3];

    Spoke topMiddle, leftMiddle, rightMiddle, botMiddle;
    InterpolateMiddleSpoke(Sp11, Sp12, lambda, &topMiddle);
    InterpolateMiddleSpoke(Sp11, Sp21, lambda, &leftMiddle);
    InterpolateMiddleSpoke(Sp21, Sp22, lambda, &botMiddle);
    InterpolateMiddleSpoke(Sp22, Sp12, lambda, &rightMiddle);

    // interpolate center spoke in this quad
    Spoke centerA, centerB, center;
    bool retA = InterpolateMiddleSpoke(&topMiddle, &botMiddle, lambda, &centerA);
    bool retB = InterpolateMiddleSpoke(&leftMiddle, &rightMiddle, lambda, &centerB);

    double rCenter = 0.0;

    double uCenter[3], uCenterA[3], uCenterB[3];
    centerA.GetDirection(uCenterA);
    centerB.GetDirection(uCenterB);
    if(retA && retB)
    {
        rCenter = 0.5 * (centerA.GetRadius() + centerB.GetRadius());
        uCenter[0] = 0.5 * (uCenterA[0] + uCenterB[0]);
        uCenter[1] = 0.5 * (uCenterA[1] + uCenterB[1]);
        uCenter[2] = 0.5 * (uCenterA[2] + uCenterB[2]);
    }
    else if(retA)
    {
        rCenter = centerA.GetRadius();
        uCenter[0] = uCenterA[0];
        uCenter[1] = uCenterA[1];
        uCenter[2] = uCenterA[2];
    }
    else if(retB)
    {
        rCenter = centerB.GetRadius();
        uCenter[0] = uCenterB[0];
        uCenter[1] = uCenterB[1];
        uCenter[2] = uCenterB[2];
    }
    else {
        return;
    }

    center.SetDirection(uCenter);
    center.SetRadius(rCenter);
    double halfDist = lambda/2;

    double interpolatedU[3]; double interpolatedR;
    if(std::abs(u - halfDist ) <= tolerance && std::abs(v - halfDist) <= tolerance)
    {
        // if the target spoke locates in the center of this quad, return center spoke
        interpolatedSpoke->SetDirection(uCenter);
        interpolatedSpoke->SetRadius(rCenter);
    }
    else if(std::abs(u) < tolerance && std::abs(v) < tolerance)
    {
        // left top corner
        Sp11->GetDirection(interpolatedU);
        interpolatedR = Sp11->GetRadius();
        interpolatedSpoke->SetDirection(interpolatedU);
        interpolatedSpoke->SetRadius(interpolatedR);
    }
    else if(std::abs(u-lambda) < tolerance && std::abs(v) < tolerance)
    {
        // left bottom corner
        Sp21->GetDirection(interpolatedU);
        interpolatedR = Sp21->GetRadius();
        interpolatedSpoke->SetDirection(interpolatedU);
        interpolatedSpoke->SetRadius(interpolatedR);
    }
    else if(std::abs(v-lambda) < tolerance && std::abs(u) < tolerance)
    {
        // right top corner
        Sp12->GetDirection(interpolatedU);
        interpolatedR = Sp12->GetRadius();
        interpolatedSpoke->SetDirection(interpolatedU);
        interpolatedSpoke->SetRadius(interpolatedR);
    }
    else if(std::abs(v-lambda) < tolerance && std::abs(u-lambda) < tolerance)
    {
        // right bottom corner
        Sp22->GetDirection(interpolatedU);
        interpolatedR = Sp22->GetRadius();
        interpolatedSpoke->SetDirection(interpolatedU);
        interpolatedSpoke->SetRadius(interpolatedR);
    }
    else if(std::abs(u - halfDist) <= tolerance && std::abs(v) < tolerance)
    {
        // on the left edge
        leftMiddle.GetDirection(interpolatedU);
        interpolatedR = leftMiddle.GetRadius();
        interpolatedSpoke->SetDirection(interpolatedU);
        interpolatedSpoke->SetRadius(interpolatedR);
    }
    else if(std::abs(u) < tolerance && std::abs(v - halfDist) < tolerance)
    {
        // on the top edge
        topMiddle.GetDirection(interpolatedU);
        interpolatedR = topMiddle.GetRadius();
        interpolatedSpoke->SetDirection(interpolatedU);
        interpolatedSpoke->SetRadius(interpolatedR);
    }
    else if(std::abs(u - halfDist) < tolerance && std::abs(v-lambda) < tolerance)
    {
        // right edge
        rightMiddle.GetDirection(interpolatedU);
        interpolatedR = rightMiddle.GetRadius();
        interpolatedSpoke->SetDirection(interpolatedU);
        interpolatedSpoke->SetRadius(interpolatedR);
    }
    else if(std::abs(u-lambda) < tolerance && std::abs(v-halfDist) < tolerance)
    {
        // bot edge
        botMiddle.GetDirection(interpolatedU);
        interpolatedR = botMiddle.GetRadius();
        interpolatedSpoke->SetDirection(interpolatedU);
        interpolatedSpoke->SetRadius(interpolatedR);
    }
    else
    {
        // construct corner spokes
        Spoke *newCorner[4];

        // if the target spoke locates in the 1st quadrant
        if(u < halfDist && v > halfDist)
        {
            newCorner[0] = &topMiddle;
            newCorner[1] = &center;
            newCorner[2] = &rightMiddle;
            newCorner[3] = Sp12;
            InterpolateQuad(newCorner, u, v-halfDist, lambda/2, interpolatedSpoke);
        }

        // if the target spoke locates in the 2nd quadrant
        else if(u < halfDist && v < halfDist)
        {
            newCorner[0] = Sp11;
            newCorner[1] = &leftMiddle;
            newCorner[2] = &center;
            newCorner[3] = &topMiddle;
            InterpolateQuad(newCorner, u, v, lambda/2, interpolatedSpoke);
        }

        // if the target spoke locates in the 3rd quadrant
        else if(u > halfDist && v < halfDist)
        {
            newCorner[0] = &leftMiddle;
            newCorner[1] = Sp21;
            newCorner[2] = &botMiddle;
            newCorner[3] = &center;
            InterpolateQuad(newCorner, u-halfDist, v, lambda/2, interpolatedSpoke);

        }
        // if the target spoke locates in the 4th quadrant
        else if(u > halfDist && v > halfDist)
        {
            newCorner[0] = &center;
            newCorner[1] = &botMiddle;
            newCorner[2] = Sp22;
            newCorner[3] = &rightMiddle;
            InterpolateQuad(newCorner, u-halfDist, v-halfDist, lambda/2, interpolatedSpoke);

        }
        // fall on the vertical axis, interpolate it in a degenerate quad: line segment
        else if(std::abs(v-halfDist) < tolerance)
        {
            newCorner[0] = &topMiddle;
            newCorner[1] = &botMiddle;
            newCorner[2] = nullptr;
            newCorner[3] = nullptr;
            InterpolateSegment(newCorner, u, 1, interpolatedSpoke);
        }
        // fall on the horizontal axis
        else if(std::abs(u-halfDist) < tolerance)
        {
            newCorner[0] = &leftMiddle;
            newCorner[1] = &rightMiddle;
            newCorner[2] = nullptr;
            newCorner[3] = nullptr;
            InterpolateSegment(newCorner, v, 1, interpolatedSpoke);
        }
    }

}

void SRepInterpolater::InterpolateSegment(Spoke **endSpokes, double dist, double lambda, Spoke *interpolatedSpoke)
{
    Spoke* start = endSpokes[0];
    Spoke* end = endSpokes[1];
    Spoke middleSpoke;
    InterpolateMiddleSpoke(start, end, lambda, &middleSpoke);
    double halfDist = lambda/2;
    if(std::abs(dist-halfDist) < tolerance)
    {
        *interpolatedSpoke = middleSpoke;
    }
    else if(dist < halfDist)
    {
        Spoke *newEnds[2];
        newEnds[0] = start;
        newEnds[1] = &middleSpoke;
        InterpolateSegment(newEnds, dist, halfDist, interpolatedSpoke);
    }
    else if(dist > halfDist)
    {
        Spoke *newEnds[2];
        newEnds[0] = &middleSpoke;
        newEnds[1] = end;
        InterpolateSegment(newEnds, dist - halfDist, halfDist, interpolatedSpoke);
    }
}

bool SRepInterpolater::InterpolateMiddleSpoke(Spoke *startS, Spoke *endS, double d, Spoke *interpolatedSpoke)
{
    // 1. compute 2nd derivative
    double startU[3], endU[3];
    startS->GetDirection(startU);
    endS->GetDirection(endU);
    if(!startS->IsValid() || !endS->IsValid())
    {
        return false;
    }

    double Uvv_start[3],Uvv_end[3];
    compute2ndDerivative(startU, endU, endU, d, Uvv_end);
    compute2ndDerivative(startU, endU, startU, 0, Uvv_start);

    // 2. compute the interpolated spoke
    double avg[3];
    startS->Add(endS, avg);
    avg[0] /= 2;
    avg[1] /= 2;
    avg[2] /= 2;

    double halfDist = d / 2;
    double uMiddle[3];
    // fix the bug if startU == endU, then uMiddle should be same with startU or endU
    if(std::abs(startU[0] - endU[0] ) <= tolerance
            && std::abs(startU[1] - endU[1] ) <= tolerance
            && std::abs(startU[2] - endU[2] ) <= tolerance)
    {
        uMiddle[0] = startU[0];
        uMiddle[1] = startU[1];
        uMiddle[2] = startU[2];
        Uvv_start[0] = 0.0;
        Uvv_start[1] = 0.0;
        Uvv_start[2] = 0.0;
        Uvv_end[0] = 0.0;
        Uvv_end[1] = 0.0;
        Uvv_end[2] = 0.0;
    }
    else {
        slerp(startU, endU, halfDist, uMiddle);
    }
    double innerProd1 = uMiddle[0] * avg[0] +
                        uMiddle[1] * avg[1] +
                        uMiddle[2] * avg[2];

    double innerProd2 = startU[0] * Uvv_start[0] +
                        startU[1] * Uvv_start[1] +
                        startU[2] * Uvv_start[2];

    double innerProd3 = endU[0] * Uvv_end[0] +
                        endU[1] * Uvv_end[1] +
                        endU[2] * Uvv_end[2];

    double interpolatedR =
            innerProd1 - halfDist*halfDist*0.25*(innerProd2 + innerProd3);

    // 3. return the interpolated
    interpolatedSpoke->SetRadius(interpolatedR);
    interpolatedSpoke->SetDirection(uMiddle);
    return true;
}
double h1(double s) { return 2*(s * s * s) - 3*(s * s) + 1; }
double h2(double s) { return -2*(s * s * s) + 3*(s * s); }
double h3(double s) { return (s * s * s) - 2*(s * s) + s; }
double h4(double s) { return (s * s * s) - (s * s); }
void SRepInterpolater::InterpolateSkeletalPoint(Spoke **cornerSpokes,double u, double v, double *output)
{
    Spoke* Sp11 = cornerSpokes[0];
    Spoke* Sp21 = cornerSpokes[1];
    Spoke* Sp22 = cornerSpokes[2];
    Spoke* Sp12 = cornerSpokes[3];

    double x11[3], x12[3], x21[3], x22[3];
    Sp11->GetSkeletalPoint(x11);
    Sp12->GetSkeletalPoint(x12);
    Sp21->GetSkeletalPoint(x21);
    Sp22->GetSkeletalPoint(x22);

    double hx[4][4];
    double hy[4][4];
    double hz[4][4];

    hx[0][0] = x11[0];          hx[0][1] = x12[0];
    hx[1][0] = x21[0];          hx[1][1] = x22[0];
    hx[2][0] = dxdu11[0];       hx[2][1] = dxdu12[0];
    hx[3][0] = dxdu21[0];       hx[3][1] = dxdu22[0];
    hx[0][2] = dxdv11[0];       hx[0][3] = dxdv12[0];
    hx[1][2] = dxdv21[0];       hx[1][3] = dxdv22[0];
    hx[2][2] = 0;               hx[2][3] = 0;
    hx[3][2] = 0;               hx[3][3] = 0;


    hy[0][0] = x11[1];          hy[0][1] = x12[1];
    hy[1][0] = x21[1];          hy[1][1] = x22[1];
    hy[2][0] = dxdu11[1];       hy[2][1] = dxdu12[1];
    hy[3][0] = dxdu21[1];       hy[3][1] = dxdu22[1];
    hy[0][2] = dxdv11[1];       hy[0][3] = dxdv12[1];
    hy[1][2] = dxdv21[1];       hy[1][3] = dxdv22[1];
    hy[2][2] = 0;               hy[2][3] = 0;
    hy[3][2] = 0;               hy[3][3] = 0;

    hz[0][0] = x11[2];       hz[0][1] = x12[2];
    hz[1][0] = x21[2];       hz[1][1] = x22[2];
    hz[2][0] = dxdu11[2];    hz[2][1] = dxdu12[2];
    hz[3][0] = dxdu21[2];    hz[3][1] = dxdu22[2];
    hz[0][2] = dxdv11[2];    hz[0][3] = dxdv12[2];
    hz[1][2] = dxdv21[2];    hz[1][3] = dxdv22[2];
    hz[2][2] = 0;            hz[2][3] = 0;
    hz[3][2] = 0;            hz[3][3] = 0;

    double hu[4], hv[4];
    hu[0] = h1(u);
    hu[1] = h2(u);
    hu[2] = h3(u);
    hu[3] = h4(u);
    hv[0] = h1(v);
    hv[1] = h2(v);
    hv[2] = h3(v);
    hv[3] = h4(v);

    // supposed computation is these
    //    vnl_double_1x1 xn = hu.transpose() * hx * hv;
    //    vnl_double_1x1 yn = hu.transpose() * hy * hv;
    //    vnl_double_1x1 zn = hu.transpose() * hz * hv;
    double huThx[4], huThy[4], huThz[4];
    huThx[0] = hu[0] * hx[0][0] + hu[1] * hx[1][0] + hu[2] * hx[2][0] + hu[3] * hx[3][0];
    huThx[1] = hu[0] * hx[0][1] + hu[1] * hx[1][1] + hu[2] * hx[2][1] + hu[3] * hx[3][1];
    huThx[2] = hu[0] * hx[0][2] + hu[1] * hx[1][2] + hu[2] * hx[2][2] + hu[3] * hx[3][2];
    huThx[3] = hu[0] * hx[0][3] + hu[1] * hx[1][3] + hu[2] * hx[2][3] + hu[3] * hx[3][3];

    huThy[0] = hu[0] * hy[0][0] + hu[1] * hy[1][0] + hu[2] * hy[2][0] + hu[3] * hy[3][0];
    huThy[1] = hu[0] * hy[0][1] + hu[1] * hy[1][1] + hu[2] * hy[2][1] + hu[3] * hy[3][1];
    huThy[2] = hu[0] * hy[0][2] + hu[1] * hy[1][2] + hu[2] * hy[2][2] + hu[3] * hy[3][2];
    huThy[3] = hu[0] * hy[0][3] + hu[1] * hy[1][3] + hu[2] * hy[2][3] + hu[3] * hy[3][3];

    huThz[0] = hu[0] * hz[0][0] + hu[1] * hz[1][0] + hu[2] * hz[2][0] + hu[3] * hz[3][0];
    huThz[1] = hu[0] * hz[0][1] + hu[1] * hz[1][1] + hu[2] * hz[2][1] + hu[3] * hz[3][1];
    huThz[2] = hu[0] * hz[0][2] + hu[1] * hz[1][2] + hu[2] * hz[2][2] + hu[3] * hz[3][2];
    huThz[3] = hu[0] * hz[0][3] + hu[1] * hz[1][3] + hu[2] * hz[2][3] + hu[3] * hz[3][3];

    output[0] = huThx[0] * hv[0] + huThx[1] * hv[1] + huThx[2] * hv[2];
    output[1] = huThy[0] * hv[0] + huThy[1] * hv[1] + huThy[2] * hv[2];
    output[2] = huThz[0] * hv[0] + huThz[1] * hv[1] + huThz[2] * hv[2];
}

void SRepInterpolater::SetCornerDxdu(double *u11, double *u21, double *u22, double *u12)
{
    dxdu11[0] = u11[0]; dxdu11[1] = u11[1]; dxdu11[2] = u11[2];
    dxdu21[0] = u21[0]; dxdu21[1] = u21[1]; dxdu21[2] = u21[2];
    dxdu22[0] = u22[0]; dxdu22[1] = u22[1]; dxdu22[2] = u22[2];
    dxdu12[0] = u12[0]; dxdu12[1] = u12[1]; dxdu12[2] = u12[2];
}

void SRepInterpolater::SetCornerDxdv(double *v11, double *v21, double *v22, double *v12)
{
    dxdv11[0] = v11[0]; dxdv11[1] = v11[1]; dxdv11[2] = v11[2];
    dxdv21[0] = v21[0]; dxdv21[1] = v21[1]; dxdv21[2] = v21[2];
    dxdv22[0] = v22[0]; dxdv22[1] = v22[1]; dxdv22[2] = v22[2];
    dxdv12[0] = v12[0]; dxdv12[1] = v12[1]; dxdv12[2] = v12[2];
}

void SRepInterpolater::compute2ndDerivative(double *startU, double *endU, double *targetU, double d, double *output)
{
    double del = 1e-5;
    double Upv1[3], Upv2[3], /*Upv3[3],*/ Upv4[3], Upv5[3];
    slerp(startU, endU, d+2*del, Upv1);
    slerp(startU, endU, d+del, Upv2);
    //slerp(startU, endU, d, Upv3);
    slerp(startU, endU, d-del, Upv4);
    slerp(startU, endU, d-2*del, Upv5);

    output[0] = 0.25 * (Upv5[0] + Upv1[0] - 2.0 * targetU[0]);
    output[1] = 0.25 * (Upv5[1] + Upv1[1] - 2.0 * targetU[1]);
    output[2] = 0.25 * (Upv5[2] + Upv1[2] - 2.0 * targetU[2]);
}

void SRepInterpolater::slerp(double *U1, double *U2, double u, double *output)
{
    double u1Tu2 = U1[0] * U2[0] + U1[1] * U2[1] + U1[2] * U2[2];
    if(u1Tu2 > 1.0)
    {
        u1Tu2 = 1.0;
    }
    else if(u1Tu2 < -1.0)
    {
        u1Tu2 = -1.0;
    }
    double phi = acos(u1Tu2);

    output[0] = ( sin((1-u)*phi)/sin(phi) )*U1[0] + ( sin(u*phi)/sin(phi) )*U2[0];
    output[1] = ( sin((1-u)*phi)/sin(phi) )*U1[1] + ( sin(u*phi)/sin(phi) )*U2[1];
    output[2] = ( sin((1-u)*phi)/sin(phi) )*U1[2] + ( sin(u*phi)/sin(phi) )*U2[2];
}
