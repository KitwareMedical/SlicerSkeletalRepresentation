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
    
    // Update spoke lengths and dirs
    // The input array is from NEWUOA, formed in order of (ux, uy, uz, x_r)
    // where x_r is the logarithm of the ratio
    void Refine(const double *coeff);
    
private:
    int nRows;
    int nCols;
    std::vector<vtkSpoke*> spokes;
};

#endif // VTKSREP_H
