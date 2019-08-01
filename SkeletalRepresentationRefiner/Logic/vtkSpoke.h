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
#ifndef VTKSPOKE_H
#define VTKSPOKE_H

#include <vector>
#include <vtkObject.h>

/**
 * @brief The vtkSpoke class
 * Each spoke consists of skeletal point, radius and direction.
 * The boundary point implied by a spoke is computed every time according to those three components above.
 */
class vtkSpoke : public vtkObject
{
public:
    static vtkSpoke* New();

    vtkSpoke(double radius, double px, double py, double pz, double ux, double uy, double uz);

    // copy constructor implements deep copy
    vtkSpoke(const vtkSpoke& src);

    // deep copy in assignment
    vtkSpoke& operator=(const vtkSpoke& other);

    void SetRadius(double r);

    void SetSkeletalPoint(double px, double py, double pz);

    void SetDirection(double *u);

    void GetDirection(double *output) const;

    // S = rU
    //void GetS(double *output) const;

    double GetRadius() const;

    void GetSkeletalPoint(double *output) const;

    // Addition between two spokes S1+S2
    void Add(vtkSpoke* another, double* output) const;

    // S_this - S_another
    // S = rU
    void Diff(vtkSpoke* another, double* output) const;

    void GetBoundaryPoint(double *output) const;

    // neighbors in u direction
    // if this neighbor is forward, use neighbor - this as difference when computing finite diff
    void SetNeighborU(const std::vector<vtkSpoke*> &neighbors, bool isForward);

    // neighbors in v direction
    void SetNeighborV(const std::vector<vtkSpoke*> &neighbors, bool isForward);

    // input delta: step size in finite difference
    // output rSrad penalty
    double GetRSradPenalty(double delta);

protected:
    vtkSpoke();
    ~vtkSpoke();
private:
    void ComputeDerivatives(std::vector<vtkSpoke*> neibors, bool isForward, double stepSize, // input
                            double *dxdu, double *dSdu, double *drdu);
private:
    double mR;
    double mPx;
    double mPy;
    double mPz;
    double mUx;
    double mUy;
    double mUz;
    bool mIsForwardU, mIsForwardV;
    std::vector<vtkSpoke*> mNeighborsU;
    std::vector<vtkSpoke*> mNeighborsV;
};

#endif // VTKSPOKE_H

