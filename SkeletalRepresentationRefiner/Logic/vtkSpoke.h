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

