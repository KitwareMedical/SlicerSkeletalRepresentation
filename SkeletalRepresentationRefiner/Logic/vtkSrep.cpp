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

#include "vtkSrep.h"
#include <math.h>
#include "vtkSpoke.h"

vtkSrep::vtkSrep()
{
    
}

vtkSrep::vtkSrep(int r, int c,  std::vector<double> &radii, std::vector<double> &dirs, std::vector<double> &skeletalPoints)
{
    nRows = r;
    nCols = c;
    for(int i = 0; i < r * c; ++i)
    {
        int idTuple = i * 3;
        vtkSpoke *s = new vtkSpoke(radii[i], skeletalPoints[idTuple], skeletalPoints[idTuple + 1], skeletalPoints[idTuple + 2],
                dirs[idTuple], dirs[idTuple + 1], dirs[idTuple + 2]);
        spokes.push_back(s);
    }
}

vtkSrep::~vtkSrep()
{
    for(int i = 0; i < spokes.size(); ++i)
    {
        if(spokes[i] == NULL)
        {
            continue;
        }
        delete spokes[i];
        spokes[i] = NULL;
    }
}

vtkSpoke *vtkSrep::GetSpoke(int r, int c) const
{
    if(spokes.empty())
    {
        return NULL;
    }
    int id = r * nCols + c;
    return spokes[id];
}

bool vtkSrep::IsEmpty() const
{
    return spokes.empty();
}

std::vector<vtkSpoke *> &vtkSrep::GetAllSpokes()
{
    return spokes;
}

void vtkSrep::Refine(const double *coeff)
{
    if(spokes.empty())
    {
        return;
    }
    for(int i = 0; i < spokes.size(); ++i)
    {
        int idx = i * 4;
        double newU[3], newR, oldR;
        newU[0] = coeff[idx];
        newU[1] = coeff[idx+1];
        newU[2] = coeff[idx+2];
        
        vtkSpoke* thisSpoke = spokes[i];
        oldR = thisSpoke->GetRadius();
        newR = exp(coeff[idx+3]) * oldR;
        
        thisSpoke->SetDirection(newU);
        thisSpoke->SetRadius(newR);
    }
}
