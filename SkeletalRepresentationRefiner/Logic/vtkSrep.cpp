#include "vtkSrep.h"
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
