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

// .NAME vtkSlicerSkeletalRepresentationRefinerLogic - slicer logic class for volumes manipulation
// .SECTION Description
// This class implements the logic of interpolating sreps

#ifndef VTKSLICERSKELETALREPRESENTATIONINTERPOLATER_H
#define VTKSLICERSKELETALREPRESENTATIONINTERPOLATER_H

class vtkSpoke;
class vtkSlicerSkeletalRepresentationInterpolater
{
public:
    vtkSlicerSkeletalRepresentationInterpolater();
    // main entry of interpolate a target spoke at (u,v)
    // Input: uBase, vBase are the range of u, v coords
    // Input: u, v coordinates in [0,1] which is the base position of the target spoke
    // Input: the array of 4 corner spokes
    // Output: the interpolated spoke
    void Interpolate(double u, double v, vtkSpoke** conerSpokes, vtkSpoke* interpolatedSpoke);
    // Interpolate the spoke length and dir for a target spoke in a quadrent given by cornerPoints of a 3d S-rep
    // Input: cornerPoints from top-left to bottom-left to bottom-right to top-right, ccw direction
    // Input: the u,v coord of the target spoke
    // Input: current interpolate level lambda
    // Output: the interpolated spoke containing direction and length
    // Return: error code
    void InterpolateQuad(vtkSpoke** cornerSpokes, double u, double v, double lambda, vtkSpoke* interpolatedSpoke);
    // Interpolate the spoke length and dir in a line segment (could be vertical or hor)
    // Input: 2 spokes designated two ends of this line segment
    // Input: the distance (dist) from the target spoke to the first spoke
    // Input: lambda is the total length between two end spokes
    // Output: the target spoke
    void InterpolateSegment(vtkSpoke** endSpokes, double dist, double lambda, vtkSpoke* interpolatedSpoke);
    // Interpolate the spoke between start spoke (startS) and end spoke(endS) given the distance (d) from
    // the start spoke to the end.
    void InterpolateMiddleSpoke(vtkSpoke* startS, vtkSpoke* endS, double d, vtkSpoke* interpolatedSpoke);
    // interpolate skeletal point
    // Input corner spokes with radius, direction and base point
    void InterpolateSkeletalPoint(vtkSpoke** cornerSpokes, double u, double v, double *output);
    void SetCornerDxdu(double *u11, double *u21, double *u22, double *u12);
    void SetCornerDxdv(double *v11, double *v21, double *v22, double *v12);
private:
    void compute2ndDerivative(double *startU, double *endU, double *targetU, double d, double *output);
    void slerp(double *U1, double *U2, double u, double* output);
    void computeDxdu(double *output);
    void computeDxdv(double *output);
private:
    double dxdu11[3];
    double dxdu12[3];
    double dxdu21[3];
    double dxdu22[3];
    double dxdv11[3];
    double dxdv12[3];
    double dxdv21[3];
    double dxdv22[3];
    const double tolerance = 1e-6;
};

#endif // VTKSLICERSKELETALREPRESENTATIONINTERPOLATER_H
