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
// This class provides functions for forward mean curvature flow
// author: Zhiyuan Liu
// Date: Sept. 2018
#ifndef __vtkForwardFlowLogic_h
#define __vtkForwardFlowLogic_h

class vtkMRMLScene;
class vtkPolyData;
class vtkForwardFlowLogic{
public:
    vtkForwardFlowLogic(){}
    ~vtkForwardFlowLogic(){}

    inline void SetScene(vtkMRMLScene* _scene){scene = _scene;}

  // flow surface to the end: either it's ellipsoidal enough or reach max_itre
  // input[dt]: delta t in each move
  // input[smooth_amount]: 0-2 double value for smooth filter
  // input[max_iter]: maximum of iteration number
  // input[freq_output]: the frequence of output model(intermediate surface mesh of flow) node to scene
  int FlowSurfaceMesh(const std::string &filename, double dt, double smooth_amount, int max_iter, int freq_output);

  // flow one step only
  // input[dt]: delta t in each move
  // input[smooth_amount]: 0-2 double value for smooth filter
  int FlowSurfaceOneStep(double dt, double smooth_amount);

    // Show fitting ellipsoid in 3D window
  // Can be called after one step flow or overall flow
  // if called after one step flow,  render the ellipsoid generated just now
  // otherwise render the ellipsoid at the end
  //int ShowFittingEllipsoid(vtkPoints* points, double curr_volume, double center[3]);

private:
  void AddModelNodeToScene(vtkPolyData* mesh, const char* modelName, bool isModelVisible, double r, double g, double b);
  
private:
    vtkMRMLScene* scene;
};
#endif
