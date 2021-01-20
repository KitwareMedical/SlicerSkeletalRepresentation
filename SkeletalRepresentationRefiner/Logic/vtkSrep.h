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

#ifndef VTKSREP_H
#define VTKSREP_H
#include <vector>

class vtkSpoke;

/**
 * @brief The vtkSrep class
 * Srep is a collection of spokes and provides utilities to conveniently compute quantities needed.
 * E.g., convert from vectors of direction, skeletal points and radii into spokes
 * OR shift spokes along their own directions.
 */
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

    std::vector<vtkSpoke *> & copyFrom(std::vector<vtkSpoke *> &source);

    // Get all skeletal points
    std::vector<double> &GetAllSkeletalPoints();

    // Update spoke lengths and dirs
    // The input array is from NEWUOA, formed in order of (ux, uy, uz, x_r)
    // where x_r is the logarithm of the ratio
    void Refine(const double *coeff);

    // Add spokes. It's ok that spokes share skeletal points with existing spokes
    void AddSpokes(std::vector<double> &radii, std::vector<double> &dirs, std::vector<double> &skeletalPoints);

    // Shift spokes along the spoke direction. Useful in transformation between legacy and new format sreps
    void ShiftSpokes(double shift);

    void DeepCopy(vtkSrep& src);

    int GetNumRows() const;

    int GetNumCols() const;

private:
    int nRows;
    int nCols;
    std::vector<vtkSpoke*> spokes;
    std::vector<double> skeletalPts;
};

#endif // VTKSREP_H
