#ifndef __vtkSRepSpokeMesh_h
#define __vtkSRepSpokeMesh_h

#include "vtkSRepSpoke.h"

#include <vtkObject.h>
#include <vtkSmartPointer.h>

#include "vtkSlicerSRepModuleMRMLExport.h"

class VTK_SLICER_SREP_MODULE_MRML_EXPORT vtkSRepSpokeMesh
  : public vtkObject
{
public:
  using IndexType = long;
  using NeighborList = std::vector<IndexType>;

  static vtkSRepSpokeMesh* New();
  ~vtkSRepSpokeMesh();

  /// Standard methods for a VTK class.
  vtkTypeMacro(vtkSRepSpokeMesh, vtkObject);

  void PrintSelf(std::ostream& os, vtkIndent indent) override;

  /// Gets the number of spokes in the mesh.
  IndexType GetNumberOfSpokes() const;

  /// Returns true if there are no spokes in the mesh, false otherwise.
  bool IsEmpty() const;

  /// Adding spokes may invalidate references
  /// If index is valid, returned pointer will never be nullptr
  /// @throws std::out_of_range exception if index < 0 or index >= GetNumberOfSpokes
  const vtkSRepSpoke* At(IndexType index) const
    VTK_EXPECTS(0 <= index && index < GetNumberOfSpokes());

  /// Adding spokes may invalidate references
  /// If index is valid, returned pointer will never be nullptr
  /// Friendly reminder that because the return type is "pointer-by-value" that something like
  /// `spokeMesh.At(i) = srepSpoke;` will not work.
  /// @throws std::out_of_range exception if index < 0 or index >=
  vtkSRepSpoke* At(IndexType index)
    VTK_EXPECTS(0 <= index && index < GetNumberOfSpokes());

  /// Gets the spoke at the given index.
  /// 
  /// Adding spokes may invalidate references
  /// If index is valid, returned pointer will never be nullptr.
  const vtkSRepSpoke* operator[](IndexType index) const;

  /// Gets the spoke at the given index.
  ///
  /// Adding spokes may invalidate references.
  /// If index is valid, returned pointer will never be nullptr.
  /// Friendly reminder that because the return type is "pointer-by-value" that something like
  /// `spokeMesh[i] = srepSpoke;` will not work.
  vtkSRepSpoke* operator[](IndexType index);

  /// Removes all spoke and neighbors from the mesh.
  void Clear();

  /// Gets the neighbors of the spoke at the given index.
  /// \throws std::out_of_range if index < 0 or index >= GetNumberOfSpokes()
  const NeighborList& GetNeighbors(IndexType index) const
    VTK_EXPECTS(0 <= index && index < GetNumberOfSpokes());

  /// Sets the neighbors of the spoke at the given index.
  /// \throws std::out_of_range if index < 0 or index >= GetNumberOfSpokes()
  void SetNeighbors(IndexType index, NeighborList neighbors)
    VTK_EXPECTS(0 <= index && index < GetNumberOfSpokes());

  /// Adds a spoke to the mesh.
  /// \throws std::invalid_argument if spoke is nullptr.
  IndexType AddSpoke(vtkSRepSpoke* spoke, NeighborList neighbors = NeighborList{})
    VTK_EXPECTS(nullptr != spoke);

  /// Sets the spoke at the given index.
  /// \throws std::invalid_argument if spoke is nullptr.
  /// \throws std::out_of_range if index < 0 or index >= GetNumberOfSpokes()
  void SetSpoke(IndexType index, vtkSRepSpoke* spoke)
    VTK_EXPECTS(0 <= index && index < GetNumberOfSpokes() && nullptr != spoke);

protected:
  vtkSRepSpokeMesh();
  vtkSRepSpokeMesh(const vtkSRepSpokeMesh&) = delete;
  vtkSRepSpokeMesh(vtkSRepSpokeMesh&&) = delete;
  vtkSRepSpokeMesh& operator=(const vtkSRepSpokeMesh&) = delete;
  vtkSRepSpokeMesh& operator=(vtkSRepSpokeMesh&&) = delete;
private:
  std::vector<vtkSmartPointer<vtkSRepSpoke>> Spokes;
  std::vector<unsigned long> SpokeObservationTags;
  std::vector<NeighborList> Neighbors;

  void onSpokeModified(vtkObject *caller, unsigned long event, void* callData);
};

#endif
