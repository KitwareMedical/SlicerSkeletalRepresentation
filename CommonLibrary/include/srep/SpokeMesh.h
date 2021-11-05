#include <vector>

#include <srep/Spoke.h>

namespace srep {

class SpokeMesh {
public:
  using IndexType = size_t;
  using NeighborList = std::vector<IndexType>;

  size_t GetNumberOfSpokes() const;
  /// Adding spokes may invalidate references
  const Spoke& at(IndexType index) const;
  /// Adding spokes may invalidate references
  Spoke& at(IndexType index);
  /// Adding spokes may invalidate references
  const Spoke& operator[](IndexType index) const;
  /// Adding spokes may invalidate references
  Spoke& operator[](IndexType index);
  void Clear();

  /// Adding spokes may invalidate references
  const NeighborList& GetNeighbors(IndexType index) const;
  /// Adding spokes may invalidate references
  NeighborList& GetNeighbors(IndexType index);

  void SetNeighbors(IndexType index, const NeighborList& neighbors);
  void SetNeighbors(IndexType index, NeighborList&& neighbors);

  IndexType AddSpoke(const Spoke& spoke);
  IndexType AddSpoke(const Spoke& spoke, const NeighborList& neighbors);
  IndexType AddSpoke(const Spoke& spoke, NeighborList&& neighbors);

private:
  std::vector<Spoke> Spokes;
  std::vector<NeighborList> Neighbors;
};

} //namespace srep
