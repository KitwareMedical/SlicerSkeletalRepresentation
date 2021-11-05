#include <srep/SpokeMesh.h>

namespace srep {

size_t SpokeMesh::GetNumberOfSpokes() const {
  return this->Spokes.size();
}

const Spoke& SpokeMesh::at(const IndexType index) const {
  return this->Spokes.at(index);
}

Spoke& SpokeMesh::at(const IndexType index) {
  return this->Spokes.at(index);
}

const Spoke& SpokeMesh::operator[](const IndexType index) const {
  return this->Spokes[index];
}

Spoke& SpokeMesh::operator[](const IndexType index) {
  return this->Spokes[index];
}

void SpokeMesh::Clear() {
  this->Spokes.clear();
  this->Neighbors.clear();
}

const SpokeMesh::NeighborList& SpokeMesh::GetNeighbors(const IndexType index) const {
  return this->Neighbors.at(index);
}

SpokeMesh::NeighborList& SpokeMesh::GetNeighbors(const IndexType index) {
  return this->Neighbors.at(index);
}

void SpokeMesh::SetNeighbors(IndexType index, const NeighborList& neighbors) {
  this->Neighbors.at(index) = neighbors;
}

void SpokeMesh::SetNeighbors(IndexType index, NeighborList&& neighbors) {
  this->Neighbors.at(index) = neighbors;
}

SpokeMesh::IndexType SpokeMesh::AddSpoke(const Spoke& spoke) {
  this->Spokes.push_back(spoke);
  this->Neighbors.push_back({});
  return this->Spokes.size() - 1;
}

SpokeMesh::IndexType SpokeMesh::AddSpoke(const Spoke& spoke, const NeighborList& neighbors) {
  this->Spokes.push_back(spoke);
  this->Neighbors.push_back(neighbors);
  return this->Spokes.size() - 1;
}

SpokeMesh::IndexType SpokeMesh::AddSpoke(const Spoke& spoke, NeighborList&& neighbors) {
  this->Spokes.push_back(spoke);
  this->Neighbors.push_back(neighbors);
  return this->Spokes.size() - 1;
}

} //namespace srep
