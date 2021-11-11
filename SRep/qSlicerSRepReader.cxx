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

// Qt includes
#include <QFileInfo>

// Slicer includes
#include "qSlicerSRepReader.h"

// Logic includes
#include <vtkSlicerApplicationLogic.h>
#include "vtkSlicerSRepLogic.h"

// MRML includes
#include "vtkMRMLMessageCollection.h"

// VTK includes
#include <vtkNew.h>
#include <vtkSmartPointer.h>

//-----------------------------------------------------------------------------
class qSlicerSRepReaderPrivate
{
  public:
  vtkSmartPointer<vtkSlicerSRepLogic> SRepLogic;
};

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_Annotations
//-----------------------------------------------------------------------------
qSlicerSRepReader::qSlicerSRepReader(QObject* _parent)
  : Superclass(_parent)
  , d_ptr(new qSlicerSRepReaderPrivate)
{
}

//-----------------------------------------------------------------------------
qSlicerSRepReader::qSlicerSRepReader(vtkSlicerSRepLogic* logic, QObject* _parent)
  : Superclass(_parent)
  , d_ptr(new qSlicerSRepReaderPrivate)
{
  this->setSRepLogic(logic);
}

//-----------------------------------------------------------------------------
qSlicerSRepReader::~qSlicerSRepReader() = default;

//-----------------------------------------------------------------------------
void qSlicerSRepReader::setSRepLogic(vtkSlicerSRepLogic* logic)
{
  Q_D(qSlicerSRepReader);
  d->SRepLogic = logic;
}

//-----------------------------------------------------------------------------
vtkSlicerSRepLogic* qSlicerSRepReader::srepLogic()const
{
  Q_D(const qSlicerSRepReader);
  return d->SRepLogic.GetPointer();
}

//-----------------------------------------------------------------------------
QString qSlicerSRepReader::description()const
{
  return "SRep";
}

//-----------------------------------------------------------------------------
qSlicerIO::IOFileType qSlicerSRepReader::fileType()const
{
  return QString("SRepFile");
}

//-----------------------------------------------------------------------------
QStringList qSlicerSRepReader::extensions()const
{
  return QStringList() << "SRep (*.srep.json)";
}

//-----------------------------------------------------------------------------
bool qSlicerSRepReader::load(const IOProperties& properties)
{
  Q_D(qSlicerSRepReader);

  // get the properties
  Q_ASSERT(properties.contains("fileName"));
  QString fileName = properties["fileName"].toString();

  QString name;
  if (properties.contains("name"))
    {
    name = properties["name"].toString();
    }

  if (d->SRepLogic.GetPointer() == nullptr)
    {
    return false;
    }

  // pass to logic to do the loading, LoadSRep will load at most one node
  const char * nodeID = d->SRepLogic->LoadSRep(fileName.toUtf8(), name.toUtf8());
  QStringList list;
  if (nodeID) {
    list << nodeID;
  }
  this->setLoadedNodes(list);
  return nodeID != nullptr;
}
