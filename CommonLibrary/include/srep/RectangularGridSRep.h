#ifndef __srep_RectangularGridSRep_h
#define __srep_RectangularGridSRep_h

#include <functional>
#include <memory>
#include <srep/MeshSRepInterface.h>
#include <srep/SkeletalPoint.h>

namespace srep {

class InvalidSkeletalGridException : public std::invalid_argument {
    using std::invalid_argument::invalid_argument;
};

/// @note direction u is row direction
/// @note direction v is column direction
class RectangularGridSRep : public MeshSRepInterface {
public:
    using SkeletalGrid = std::vector<std::vector<SkeletalPoint>>;

    RectangularGridSRep();
    /// Construct from a 2D grid of skeleton points
    /// \throws InvalidSkeletalGridException if all the interior vectors are not the same size
    RectangularGridSRep(const SkeletalGrid& skeleton);
    /// Construct from a 2D grid of skeleton points
    /// \throws InvalidSkeletalGridException if all the interior vectors are not the same size
    RectangularGridSRep(SkeletalGrid&& skeleton);

    //copy and move defined and valid
    RectangularGridSRep(const RectangularGridSRep&) = delete;
    RectangularGridSRep& operator=(const RectangularGridSRep&) = delete;
    RectangularGridSRep(RectangularGridSRep&&) = delete;
    RectangularGridSRep& operator=(RectangularGridSRep&&) = delete;
    ~RectangularGridSRep() = default;

    util::owner<RectangularGridSRep*> Clone() const override;

    /// Returns true if RectangularGridSRep has no rows or columns.
    bool IsEmpty() const override;

    /// Gets the number of rows.
    size_t GetNumRows() const;

    /// Gets the number of columns.
    size_t GetNumCols() const;

    inline size_t GetNumCrestPoints() const {
        return NumCrestPoints(this->GetNumRows(), this->GetNumCols());
    }

    inline size_t GetNumSkeletalPoints() const {
        return this->GetNumRows() * this->GetNumCols();
    }

    /// Gets the grid of skeleton points.
    const SkeletalGrid& GetSkeletalPoints() const;
    /// Gets skeletal point without bounds checking
    const SkeletalPoint& GetSkeletalPoint(size_t row, size_t col) const;
    /// Gets skeletal point with bounds checking
    const SkeletalPoint& GetSkeletalPointAt(size_t row, size_t col) const;

    // TODO: is this actualy right if there is only 1 col?
    // is it even valid to have only one col?
    /// Gets the number of crest points an RectangularGridSRep will have given the number of rows and columns.
    static size_t NumCrestPoints(size_t rows, size_t cols) {
        return rows * 2 + (cols-2) * 2;
    }

    ///////////////////////////////////////////////////////////
    //
    // MeshSRepInterface overrides
    //
    ///////////////////////////////////////////////////////////
    const SpokeMesh& GetUpSpokes() const override;
    const SpokeMesh& GetDownSpokes() const override;
    const SpokeMesh& GetCrestSpokes() const override;
    const std::vector<UpDownIndices>& GetCrestSkeletalConnections() const override;
    const std::vector<UpDownIndices>& GetSpine() const override;

private:
    static void Validate(const SkeletalGrid& skeleton);
    void CreateMeshRepresentation();

    SkeletalGrid Skeleton;

    struct MeshRepresentation {
        SpokeMesh UpSpokes;
        SpokeMesh DownSpokes;
        SpokeMesh CrestSpokes;
        std::vector<UpDownIndices> CrestSkeletalConnections;
        std::vector<UpDownIndices> Spine;
    };

    MeshRepresentation SkeletonAsMesh;
};

bool operator==(const RectangularGridSRep& a, const RectangularGridSRep& b);
std::ostream& operator<<(std::ostream& os, const RectangularGridSRep& srep);

/// Create a RectangularGridSRep.
///
/// Crest spokes will be attached to the closest skeletal point.
///
/// @param rows Number of rows in the new RectangularGridSRep.
/// @param cols Number of columns in the new RectangularGridSRep.
/// @param upSpokes Spokes in the up direction. Vector must be rows * cols in length, laid out in row major order.
///                 Each spoke must have the same skeletal point as the spoke in downSpokes at the same index.
/// @param downSpokes Spokes in the down direction. Vector must be rows * cols in length, laid out in row major order.
///                   Each spoke must have the same skeletal point as the spoke in upSpokes at the same index.
/// @param crestSpokes Crest spokes. Vector must be RectangularGridSRep::NumCrestPoints(rows, cols) in length.
/// @throws std::invalid_argument if parameters do not match up as described above.
/// @throws std::runtime_error if not all crest points have a crest spoke attached to them. This can happen if
///         a single skeletal point is the closest point to two crest spokes.
std::unique_ptr<RectangularGridSRep> MakeRectangularGridSRep(size_t rows, size_t cols, const std::vector<Spoke>& upSpokes, const std::vector<Spoke>& downSpokes, const std::vector<Spoke>& crestSpokes);

/// Runs a function over every crest point in an RectangularGridSRep skeletal grid.
void foreachCrestPoint(RectangularGridSRep::SkeletalGrid& grid, const std::function<void(SkeletalPoint& crestPoint)>& call);
/// Runs a function over every crest point in an RectangularGridSRep skeletal grid.
void foreachCrestPoint(const RectangularGridSRep::SkeletalGrid& grid, const std::function<void(const SkeletalPoint& crestPoint)>& call);
/// Runs a function over every point in an RectangularGridSRep skeletal grid.
void foreachPoint(RectangularGridSRep::SkeletalGrid& grid, const std::function<void(SkeletalPoint& point)>& call);
/// Runs a function over every point in an RectangularGridSRep skeletal grid.
void foreachPoint(const RectangularGridSRep::SkeletalGrid& grid, const std::function<void(const SkeletalPoint& point)>& call);

/// Runs a function over every crest point in an RectangularGridSRep.
inline void foreachCrestPoint(const RectangularGridSRep& srep, const std::function<void(const SkeletalPoint& crestPoint)>& call) {
    foreachCrestPoint(srep.GetSkeletalPoints(), call);
}
/// Runs a function over every point in an RectangularGridSRep.
inline void foreachPoint(const RectangularGridSRep& srep, const std::function<void(const SkeletalPoint& point)>& call) {
    foreachPoint(srep.GetSkeletalPoints(), call);
}

}

#endif
