#ifndef __vtkMeshSRepInterface_h
#define __vtkMeshSRepInterface_h

#include "vtkSRepSpokeMesh.h"

#include <vtkObject.h>
#include <vtkSmartPointer.h>

#include "vtkSlicerSRepModuleMRMLExport.h"

// not exporting because it is an interface that cannot be instantiated directly
class VTK_SLICER_SREP_MODULE_MRML_EXPORT vtkMeshSRepInterface
  : public vtkObject
{
public:
  /// Standard methods for a VTK class.
  vtkAbstractTypeMacro(vtkMeshSRepInterface, vtkObject);

  using IndexType = vtkSRepSpokeMesh::IndexType;

  /// Makes a deep copy clone of object.
  ///
  /// \returns Owning pointer to deep copy clone of the object.
  virtual VTK_NEWINSTANCE vtkMeshSRepInterface* Clone() const = 0;

  /// Returns true if there are no spokes in the SRep, false otherwise.
  virtual bool IsEmpty() const = 0;

  /// Gets all spokes in the up orientation.
  ///
  /// Will never return nullptr. If there are no spokes in this orientation,
  /// an empty vtkSRepSpokeMesh will be returned.
  virtual const vtkSRepSpokeMesh* GetUpSpokes() const = 0;

  /// Gets all spokes in the down orientation.
  ///
  /// Will never return nullptr. If there are no spokes in this orientation,
  /// an empty vtkSRepSpokeMesh will be returned.
  virtual const vtkSRepSpokeMesh* GetDownSpokes() const = 0;

  /// Gets all spokes in the crest orientation.
  ///
  /// Will never return nullptr. If there are no spokes in this orientation,
  /// an empty vtkSRepSpokeMesh will be returned.
  virtual const vtkSRepSpokeMesh* GetCrestSpokes() const = 0;

  /// Gets the connections from the crest to the skeleton.
  ///
  /// \returns Connections from crest to skeleton. The index in the list corresponds to
  ///          the index into the value of GetCrestSpokes(). The value at each index is
  ///          the index into the value of GetUpSpokes().
  virtual const std::vector<IndexType>& GetCrestToUpSpokeConnections() const = 0;

  /// Gets the connections from the crest to the skeleton.
  ///
  /// \returns Connections from crest to skeleton. The index in the list corresponds to
  ///          the index into the value of GetCrestSpokes(). The value at each index is
  ///          the index into the value of GetUpSpokes().
  virtual const std::vector<IndexType>& GetCrestToDownSpokeConnections() const = 0;

  /// Gets the indices of the spine.
  ///
  /// \returns The spine. The value of each index is an index into the value of
  ///          GetUpSpokes(). The order is the spine connection. I.e.
  ///          GetUpSpine()[0] connects to GetUpSpine()[1] which connects to GetUpSpine()[2]
  ///          and so on to create the whole spine.
  virtual const std::vector<IndexType>& GetUpSpine() const = 0;

  /// Gets the indices of the spine.
  ///
  /// \returns The spine. The value of each index is an index into the value of
  ///          GetDownSpokes(). The order is the spine connection. I.e.
  ///          GetDownSpine()[0] connects to GetDownSpine()[1] which connects to GetDownSpine()[2]
  ///          and so on to create the whole spine.
  virtual const std::vector<IndexType>& GetDownSpine() const = 0;
};

#endif
