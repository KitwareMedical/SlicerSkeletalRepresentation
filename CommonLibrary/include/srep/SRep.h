#ifndef __srep_SRep_h
#define __srep_SRep_h

#include <functional>
#include <srep/SkeletalPoint.h>

namespace srep {

class InvalidSkeletalGridException : public std::invalid_argument {
    using std::invalid_argument::invalid_argument;
};

/// @note direction u is row direction
/// @note direction v is column direction
class SRep {
public:
    using SkeletalGrid = std::vector<std::vector<SkeletalPoint>>;

    SRep();
    /// Construct from a 2D grid of skeleton points
    /// \throws InvalidSkeletalGridException if all the interior vectors are not the same size
    SRep(const SkeletalGrid& skeleton);
    /// Construct from a 2D grid of skeleton points
    /// \throws InvalidSkeletalGridException if all the interior vectors are not the same size
    SRep(SkeletalGrid&& skeleton);

    //copy and move defined and valid
    SRep(const SRep&) = default;
    SRep& operator=(const SRep&) = default;
    SRep(SRep&&) = default;
    SRep& operator=(SRep&&) = default;
    ~SRep() = default;

    /// Returns true if SRep has no rows or columns.
    bool IsEmpty() const;

    /// Gets the number of rows.
    size_t GetNumRows() const;

    /// Gets the number of columns.
    size_t GetNumCols() const;

    /// Gets the grid of skeleton points.
    const SkeletalGrid& GetSkeletalPoints() const;
    /// Gets skeletal point without bounds checking
    const SkeletalPoint& GetSkeletalPoint(size_t row, size_t col) const;
    /// Gets skeletal point with bounds checking
    const SkeletalPoint& GetSkeletalPointAt(size_t row, size_t col) const;

    // TODO: is this actualy right if there is only 1 col?
    // is it even valid to have only one col?
    /// Gets the number of crest points an SRep will have given the number of rows and columns.
    static size_t NumCrestPoints(size_t rows, size_t cols) {
        return rows * 2 + (cols-2) * 2;
    }
private:
    static void Validate(const SkeletalGrid& skeleton);

    SkeletalGrid Skeleton;
};

/// Create a SRep.
///
/// Crest spokes will be attached to the closest skeletal point.
///
/// @param rows Number of rows in the new SRep.
/// @param cols Number of columns in the new SRep.
/// @param upSpokes Spokes in the up direction. Vector must be rows * cols in length, laid out in row major order.
///                 Each spoke must have the same skeletal point as the spoke in downSpokes at the same index.
/// @param downSpokes Spokes in the down direction. Vector must be rows * cols in length, laid out in row major order.
///                   Each spoke must have the same skeletal point as the spoke in upSpokes at the same index.
/// @param crestSpokes Crest spokes. Vector must be SRep::NumCrestPoints(rows, cols) in length.
/// @throws std::invalid_argument if parameters do not match up as described above.
/// @throws std::runtime_error if not all crest points have a crest spoke attached to them. This can happen if
///         a single skeletal point is the closest point to two crest spokes.
SRep MakeSRep(size_t rows, size_t cols, const std::vector<Spoke>& upSpokes, const std::vector<Spoke>& downSpokes, const std::vector<Spoke>& crestSpokes);

/// Runs a function over every crest point in an SRep skeletal grid.
void foreachCrestPoint(SRep::SkeletalGrid& grid, const std::function<void(SkeletalPoint& crestPoint)>& call);
/// Runs a function over every crest point in an SRep skeletal grid.
void foreachCrestPoint(const SRep::SkeletalGrid& grid, const std::function<void(const SkeletalPoint& crestPoint)>& call);
/// Runs a function over every point in an SRep skeletal grid.
void foreachPoint(SRep::SkeletalGrid& grid, const std::function<void(SkeletalPoint& point)>& call);
/// Runs a function over every point in an SRep skeletal grid.
void foreachPoint(const SRep::SkeletalGrid& grid, const std::function<void(const SkeletalPoint& point)>& call);

/// Runs a function over every crest point in an SRep.
inline void foreachCrestPoint(const SRep& srep, const std::function<void(const SkeletalPoint& crestPoint)>& call) {
    foreachCrestPoint(srep.GetSkeletalPoints(), call);
}
/// Runs a function over every point in an SRep.
inline void foreachPoint(const SRep& srep, const std::function<void(const SkeletalPoint& point)>& call) {
    foreachPoint(srep.GetSkeletalPoints(), call);
}

}

#endif