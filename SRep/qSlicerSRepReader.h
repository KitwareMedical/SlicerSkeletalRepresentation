/*==============================================================================

  Program: 3D Slicer

  Copyright (c) BWH

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

==============================================================================*/

#ifndef __qSlicerSRepReader
#define __qSlicerSRepReader

// Slicer includes
#include "qSlicerFileReader.h"

class qSlicerSRepReaderPrivate;
class vtkSlicerSRepLogic;

//----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_SRep
class qSlicerSRepReader
  : public qSlicerFileReader
{
  Q_OBJECT
public:
  typedef qSlicerFileReader Superclass;
  qSlicerSRepReader(QObject* parent = nullptr);
  /// Creates a qSlicerSRepReader
  ///
  /// An ownership claim on logic is taken.
  qSlicerSRepReader(vtkSlicerSRepLogic* logic, QObject* parent = nullptr);
  ~qSlicerSRepReader() override;

  vtkSlicerSRepLogic* srepLogic()const;
  void setSRepLogic(vtkSlicerSRepLogic* logic);

  QString description()const override;
  IOFileType fileType()const override;
  QStringList extensions()const override;

  bool load(const IOProperties& properties) override;

protected:
  QScopedPointer<qSlicerSRepReaderPrivate> d_ptr;

private:
  Q_DECLARE_PRIVATE(qSlicerSRepReader);
  Q_DISABLE_COPY(qSlicerSRepReader);
};

#endif
