#ifndef __srep_EllipticalSRep_h
#define __srep_EllipticalSRep_h

#include <srep/MeshSRepInterface.h>
#include <srep/SkeletalPoint.h>
#include <srep/Exceptions.h>

namespace srep {

/// An SRep with an elliptical grid
///
/// Grid is somewhat elliptical, but not in polar coordinates.
/// Instead of all lines coming from outer boundary of the ellipse
/// meeting at the center, there is a "spine" in the ellipse that
/// the lines go to.
///
/// Instead of a polar system like this (where all lines meet at the center .)
///   \ | /
/// --- . ---
///   / | \                (characters to prevent multi-line comment)
///
/// It looks more like this (. == spine)
///       \  |  /
///  --- . . . . . ---
///       /  |  \          (characters to prevent multi-line comment)
///
/// where each point on the spine has 2 lines going out of it,
/// with the exception of the poles on the longer axis, which have 1
class EllipticalSRep : public MeshSRepInterface {
public:
  /// Unrolled elliptical grid. first "row" is the spine.
  /// There will be duplicate spine points as there is a "column" for
  /// both sides of the spine. These duplicate points will not show
  /// up as two separate points in the mesh representation.
  /// The last "row" in the grid will be the crest (aka fold).
  /// First line is expected to be one of the poles
  using LineOutFromSpine = std::vector<SkeletalPoint>;
  using UnrolledEllipticalGrid = std::vector<LineOutFromSpine>;
  EllipticalSRep();
  EllipticalSRep(UnrolledEllipticalGrid skeleton);
  ~EllipticalSRep() = default;

  //copy and move deleted
  EllipticalSRep(const EllipticalSRep&) = delete;
  EllipticalSRep& operator=(const EllipticalSRep&) = delete;
  EllipticalSRep(EllipticalSRep&&) = delete;
  EllipticalSRep& operator=(EllipticalSRep&&) = delete;

  const UnrolledEllipticalGrid& GetSkeleton() const;

  ///////////////////////////////////////////////////////////
  //
  // MeshSRepInterface overrides
  //
  ///////////////////////////////////////////////////////////
  util::owner<EllipticalSRep*> Clone() const override;
  bool IsEmpty() const override;
  const SpokeMesh& GetUpSpokes() const override;
  const SpokeMesh& GetDownSpokes() const override;
  const SpokeMesh& GetCrestSpokes() const override;
  const std::vector<UpDownIndices>& GetCrestSkeletalConnections() const override;
  const std::vector<UpDownIndices>& GetSpine() const override;
private:
  static void Validate(const UnrolledEllipticalGrid& skeleton);
  void CreateMeshRepresentation();

  UnrolledEllipticalGrid Skeleton;

  struct MeshRepresentation {
      SpokeMesh UpSpokes;
      SpokeMesh DownSpokes;
      SpokeMesh CrestSpokes;
      std::vector<UpDownIndices> CrestSkeletalConnections;
      std::vector<UpDownIndices> Spine;
  };

  MeshRepresentation SkeletonAsMesh;
};

}

#endif
