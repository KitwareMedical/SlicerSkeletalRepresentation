#ifndef __srep_MeshSRepInterface_h
#define __srep_MeshSRepInterface_h

#include <srep/Spoke.h>
#include <srep/SpokeMesh.h>
#include <srep/Util.h>

namespace srep {

/// The most general interface of an SRep
///
/// 
class MeshSRepInterface {
public:
  using IndexType = SpokeMesh::IndexType;

  //disallow copy and move of the interface
  MeshSRepInterface() = default;
  MeshSRepInterface(const MeshSRepInterface&) = delete;
  MeshSRepInterface& operator=(const MeshSRepInterface&) = delete;
  MeshSRepInterface(MeshSRepInterface&&) = delete;
  MeshSRepInterface& operator=(MeshSRepInterface&&) = delete;

  virtual ~MeshSRepInterface() = default;

  /// Clones the SRep.
  ///
  /// Deep copy is made. The reason that owner is used instead of a unique_ptr is to allow
  /// the return type to be covariant.
  virtual util::owner<MeshSRepInterface*> Clone() const = 0;

  virtual bool IsEmpty() const = 0;

  /// Mesh of up spokes.
  ///
  /// This is a parallel mesh to GetDownSpokes.
  /// \sa GetDownSpokes, GetCrestSpokes
  virtual const SpokeMesh& GetUpSpokes() const = 0;

  /// Mesh of down spokes.
  ///
  /// This is a parallel mesh to GetDownSpokes.
  /// \sa GetUpSpokes, GetCrestSpokes
  virtual const SpokeMesh& GetDownSpokes() const = 0;

  /// Mesh of crest spokes.
  ///
  /// \sa GetUpSpokes, GetDownSpokes
  virtual const SpokeMesh& GetCrestSpokes() const = 0;

  /// Gets the connections from the crest to the skeleton.
  ///
  /// This is a parallel list to GetCrestSpokes.
  /// \sa GetCrestSpokes
  virtual const std::vector<IndexType>& GetCrestSkeletalConnections() const = 0;

  /// Indices out of GetUpSpokes and GetDownSpokes that make up the spine.
  virtual const std::vector<IndexType>& GetSpine() const = 0;
};

};

#endif
