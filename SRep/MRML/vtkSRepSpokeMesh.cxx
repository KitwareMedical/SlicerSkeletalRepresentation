#include <vtkSRepSpokeMesh.h>

#include <vtkCommand.h>
#include <vtkObjectFactory.h>

//----------------------------------------------------------------------
vtkStandardNewMacro(vtkSRepSpokeMesh);

//----------------------------------------------------------------------
vtkSRepSpokeMesh::vtkSRepSpokeMesh() = default;

//----------------------------------------------------------------------
vtkSRepSpokeMesh::~vtkSRepSpokeMesh() {
  this->Clear();
}

//----------------------------------------------------------------------
void vtkSRepSpokeMesh::PrintSelf(ostream& os, vtkIndent indent) {
  this->Superclass::PrintSelf(os, indent);
  os << indent << "Number of Spokes: " << this->GetNumberOfSpokes();
}

//----------------------------------------------------------------------
vtkSRepSpokeMesh::IndexType vtkSRepSpokeMesh::GetNumberOfSpokes() const {
  return static_cast<IndexType>(this->Spokes.size());
}

//----------------------------------------------------------------------
bool vtkSRepSpokeMesh::IsEmpty() const {
  return 0 == this->GetNumberOfSpokes();
}

//----------------------------------------------------------------------
const vtkSRepSpoke* vtkSRepSpokeMesh::At(IndexType index) const {
  return this->Spokes.at(index);
}

//----------------------------------------------------------------------
vtkSRepSpoke* vtkSRepSpokeMesh::At(IndexType index) {
  return this->Spokes.at(index);
}

//----------------------------------------------------------------------
const vtkSRepSpoke* vtkSRepSpokeMesh::operator[](const IndexType index) const {
  return this->Spokes[index];
}

//----------------------------------------------------------------------
vtkSRepSpoke* vtkSRepSpokeMesh::operator[](const IndexType index) {
  return this->Spokes[index];
}

//----------------------------------------------------------------------
void vtkSRepSpokeMesh::Clear() {
  for (size_t i = 0; i < this->Spokes.size(); ++i) {
    this->Spokes[i]->RemoveObserver(this->SpokeObservationTags[i]);
  }
  this->Spokes.clear();
  this->SpokeObservationTags.clear();
  this->Neighbors.clear();
  this->Modified();
}

//----------------------------------------------------------------------
const vtkSRepSpokeMesh::NeighborList& vtkSRepSpokeMesh::GetNeighbors(const IndexType index) const {
  return this->Neighbors.at(index);
}

//----------------------------------------------------------------------
void vtkSRepSpokeMesh::SetNeighbors(IndexType index, NeighborList neighbors) {
  this->Neighbors.at(index) = std::move(neighbors);
  this->Modified();
}

//----------------------------------------------------------------------
vtkSRepSpokeMesh::IndexType vtkSRepSpokeMesh::AddSpoke(vtkSRepSpoke* spoke, NeighborList neighbors) {
  if (!spoke) {
    throw std::invalid_argument("Cannot add nullptr spoke to SpokeMesh");
  }

  this->Spokes.push_back(spoke);
  this->SpokeObservationTags.push_back(this->Spokes.back()->AddObserver(vtkCommand::ModifiedEvent, this, &vtkSRepSpokeMesh::onSpokeModified));
  this->Neighbors.push_back(std::move(neighbors));
  this->Modified();
  return this->Spokes.size() - 1;
}

//----------------------------------------------------------------------
void vtkSRepSpokeMesh::SetSpoke(IndexType index, vtkSRepSpoke* spoke) {
  if (!spoke) {
    throw std::invalid_argument("Cannot add nullptr spoke to SpokeMesh");
  }
  this->Spokes.at(index)->RemoveObserver(this->SpokeObservationTags.at(index));
  this->Spokes[index] = spoke;
  this->SpokeObservationTags[index] = this->Spokes[index]->AddObserver(vtkCommand::ModifiedEvent, this, &vtkSRepSpokeMesh::onSpokeModified);
  this->Modified();
}

//----------------------------------------------------------------------
void vtkSRepSpokeMesh::onSpokeModified(vtkObject */*caller*/, unsigned long /*event*/, void* /*callData*/) {
  this->Modified();
}
