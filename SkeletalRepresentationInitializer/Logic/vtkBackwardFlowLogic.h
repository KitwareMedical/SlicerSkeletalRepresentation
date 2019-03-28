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
// This class provides logic of backward flow
// author: Zhiyuan Liu
// Date: Sept. 4, 2018
#ifndef __vtkBackwardFlowLogic_h
#define __vtkBackwardFlowLogic_h

class vtkPolyData;
class vtkBackwardFlowLogic {
public:
    vtkBackwardFlowLogic(){}
    ~vtkBackwardFlowLogic(){}

    void runApplyTPS();
    void computePairwiseTPS(vtkPolyData* afterFlow, vtkPolyData* beforeFlow, const char* outputFileName);
    void generateEllipsoidSrep(int numRow, int numCol, double ra, double rb, double rc, const char* outputPath);
};
#endif
