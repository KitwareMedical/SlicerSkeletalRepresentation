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

#include <vtkCommand.h>
#include <vtkObject.h>

#ifndef __srepCreator_SRepProgressHelper_h
#define __srepCreator_SRepProgressHelper_h

template <class PROGRESS_WIDGET>
class SRepProgressHelper {
public:
  using WidgetType = PROGRESS_WIDGET;
  SRepProgressHelper(vtkObject& object, WidgetType* widget)
    : Object(object)
    , ObservationTag(Object.AddObserver(vtkCommand::ProgressEvent, this, &SRepProgressHelper::OnProgress))
    , Widget(widget)
  {
    if (this->Widget) {
      this->Widget->setMinimum(0);
      this->Widget->setMaximum(100);
      this->Widget->setValue(0);
    }
  }
  ~SRepProgressHelper() {
    this->Object.RemoveObserver(this->ObservationTag);
  }
private:
  void OnProgress(vtkObject *caller, unsigned long event, void* callData) {
    if (this->Widget && caller == &(this->Object) && event == vtkCommand::ProgressEvent) {
      double progress = *reinterpret_cast<double*>(callData);
      this->Widget->setValue(static_cast<int>(progress * 100));
    } else {
      std::cerr << "Unexpected event callback in SRepProgressHelper" << std::endl;
    }
  }

  vtkObject& Object;
  unsigned long ObservationTag;
  WidgetType* Widget;
};

#endif
