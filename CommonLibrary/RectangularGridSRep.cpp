#include <srep/RectangularGridSRep.h>
#include <iostream>
#include <limits>
#include <memory>

namespace srep {

void foreachCrestPoint(RectangularGridSRep::SkeletalGrid& grid, const std::function<void(SkeletalPoint& crestPoint)>& call) {
    if (grid.size() == 0) {
        return;
    }
    //top row
    for (auto& point : grid.front()) {
        call(point);
    }
    //outside cols on interior rows
    //starting at 1 ensures we skip the first one row
    //ending at size - 1 ensures we skip the last row
    for (size_t i = 1; i < grid.size() - 1; ++i) {
        if (grid[i].size() > 0) {
            //first col
            call(grid[i].front());
            //last col if we have more than 1 cols
            if (grid[i].size() > 1) {
                call(grid[i].back());
            }
        }
    }
    //bottom row if we have more than 1 row
    if (grid.size() > 1) {
        for (auto& point : grid.back()) {
            call(point);
        }
    }
}

void foreachCrestPoint(const RectangularGridSRep::SkeletalGrid& grid, const std::function<void(const SkeletalPoint& crestPoint)>& call) {
    // yes, const_cast is bad, but we know that foreachCrestPoint (non-const) does not modify the grid, it only
    // leaves the option for the callback function to modify the point, which we do not allow by routing through
    // the wrapper callback here.
    foreachCrestPoint(const_cast<RectangularGridSRep::SkeletalGrid&>(grid), [&](SkeletalPoint& crestPoint) {
        call(crestPoint);
    });
}

void foreachPoint(RectangularGridSRep::SkeletalGrid& grid, const std::function<void(SkeletalPoint& point)>& call) {
    for (auto& row : grid) {
        for (auto& skeletalPoint : row) {
            call(skeletalPoint);
        }
    }
}

void foreachPoint(const RectangularGridSRep::SkeletalGrid& grid, const std::function<void(const SkeletalPoint& point)>& call) {
    for (auto& row : grid) {
        for (auto& skeletalPoint : row) {
            call(skeletalPoint);
        }
    }
}


RectangularGridSRep::RectangularGridSRep()
    : Skeleton()
    , SkeletonAsMesh()
{}

RectangularGridSRep::RectangularGridSRep(const SkeletalGrid& skeleton)
    : Skeleton()
    , SkeletonAsMesh()
{
    this->Validate(skeleton);
    this->Skeleton = skeleton;
    this->CreateMeshRepresentation();
}

RectangularGridSRep::RectangularGridSRep(SkeletalGrid&& skeleton)
    : Skeleton()
    , SkeletonAsMesh()
{
    this->Validate(skeleton);
    this->Skeleton = std::move(skeleton);
    this->CreateMeshRepresentation();
}

void RectangularGridSRep::Validate(const SkeletalGrid& skeleton) {
    //validate it is a valid grid
    if (!skeleton.empty()) {
        const auto numCols = skeleton.front().size();
        if (numCols == 0) {
            throw InvalidSkeletalGridException("SkeletalGrid must have at minimum one column per row");
        }
        for (const auto& row : skeleton) {
            if (row.size() != numCols) {
                throw InvalidSkeletalGridException("All rows of SkeletalGrid must have same number of columns");
            }
        }
    }

    //validate that the crest is correct
    for (size_t row = 0; row < skeleton.size(); ++row) {
        for (size_t col = 0; col < skeleton[row].size(); ++col) {
            if (row == 0 || col == 0 || row == skeleton.size() - 1 || col ==  skeleton[row].size() - 1) {
                if (!skeleton[row][col].IsCrest()) {
                    throw InvalidSkeletalGridException("Outer boundary are not all crest points");
                }
            } else {
                if (skeleton[row][col].IsCrest()) {
                    throw InvalidSkeletalGridException("Interior of the grid should not contain crest points");
                }
            }
        }
    }
}

util::owner<RectangularGridSRep*> RectangularGridSRep::Clone() const {
    // use a unique ptr right up until the very end in case there are any exceptions
    std::unique_ptr<RectangularGridSRep> cloned(new RectangularGridSRep);
    cloned->Skeleton = this->Skeleton;
    cloned->SkeletonAsMesh = this->SkeletonAsMesh;
    return cloned.release();
}

bool RectangularGridSRep::IsEmpty() const {
    return this->Skeleton.empty();
}

size_t RectangularGridSRep::GetNumRows() const {
    return this->Skeleton.size();
}

size_t RectangularGridSRep::GetNumCols() const {
    return this->Skeleton.empty() ? 0 : this->Skeleton.front().size();
}

const RectangularGridSRep::SkeletalGrid& RectangularGridSRep::GetSkeletalPoints() const {
    return this->Skeleton;
}

const SkeletalPoint& RectangularGridSRep::GetSkeletalPoint(size_t row, size_t col) const {
    return this->Skeleton[row][col];
}

const SkeletalPoint& RectangularGridSRep::GetSkeletalPointAt(size_t row, size_t col) const {
    return this->Skeleton.at(row).at(col);
}

void RectangularGridSRep::CreateMeshRepresentation() {
    this->SkeletonAsMesh.UpSpokes.Clear();
    this->SkeletonAsMesh.DownSpokes.Clear();
    this->SkeletonAsMesh.CrestSpokes.Clear();
    this->SkeletonAsMesh.CrestSkeletalConnections.clear();
    this->SkeletonAsMesh.Spine.clear();

    /////////////////////////////////////////////////////////////////////////////
    //
    // UpSpokes, DownSpokes, and Spine
    //
    /////////////////////////////////////////////////////////////////////////////
    // UpSpokes/DownSpokes will be row major
    const auto toUpDownMeshIndex = [this](const size_t row, const size_t col) {
        return row * this->GetNumCols() + col;
    };

    for (size_t row = 0; row < this->GetNumRows(); ++row) {
        for (size_t col = 0; col < this->GetNumCols(); ++col) {
            std::vector<IndexType> neighbors;
            if (row > 0) neighbors.push_back(toUpDownMeshIndex(row - 1, col));
            if (row < this->GetNumRows() - 1) neighbors.push_back(toUpDownMeshIndex(row + 1, col));
            if (col > 0) neighbors.push_back(toUpDownMeshIndex(row, col - 1));
            if (col < this->GetNumCols() - 1) neighbors.push_back(toUpDownMeshIndex(row, col + 1));

            const auto& skeletalPoint = this->Skeleton[row][col];

            this->SkeletonAsMesh.UpSpokes.AddSpoke(skeletalPoint.GetUpSpoke(), neighbors);
            this->SkeletonAsMesh.DownSpokes.AddSpoke(skeletalPoint.GetDownSpoke(), neighbors);
        }
    }

    //Spine is the middle row if there is an even number of rows, then no spine?
    if (this->GetNumRows() % 2 == 1) {
        const auto spineRow = this->GetNumRows() / 2;
        for (size_t col = 0; col < this->GetNumCols(); ++col) {
            const auto index = toUpDownMeshIndex(spineRow, col);
            // for this class the up and down spokes are at equivalent indices
            this->SkeletonAsMesh.Spine.push_back(std::make_pair(index, index));
        }
    }

    /////////////////////////////////////////////////////////////////////////////
    //
    // CrestSpokes and CrestSkeletalConnections
    //
    /////////////////////////////////////////////////////////////////////////////
    // CrestSpokes will be clockwise from top left
    const auto numCrestPoints = this->GetNumCrestPoints();
    const auto addCrestPoint = [this, numCrestPoints, toUpDownMeshIndex](const size_t row, const size_t col) {
        // neighbors will be the prev and next index wrapping at the end.
        // this->SkeletonAsMesh.CrestSpokes.GetNumberOfSpokes() is the index to being added
        std::vector<IndexType> neighbors;
        neighbors.push_back((numCrestPoints + this->SkeletonAsMesh.CrestSpokes.GetNumberOfSpokes() - 1) % numCrestPoints);
        neighbors.push_back((numCrestPoints + this->SkeletonAsMesh.CrestSpokes.GetNumberOfSpokes() + 1) % numCrestPoints);

        const auto skeletalPoint = this->Skeleton[row][col];
        this->SkeletonAsMesh.CrestSpokes.AddSpoke(skeletalPoint.GetCrestSpoke(), neighbors);
        const auto index = toUpDownMeshIndex(row, col);
        this->SkeletonAsMesh.CrestSkeletalConnections.push_back(std::make_pair(index, index));
    };

    //top row, traverse right
    for (size_t row = 0; row < this->GetNumRows(); ++row) {
        addCrestPoint(row, 0);
    }
    //rightmost col, except top and bottom, traverse down
    for (size_t col = 1; col < this->GetNumCols() - 1; ++col) {
        addCrestPoint(this->GetNumRows() - 1, col);
    }
    //bottom row, traverse left (long long to allow negative for loop break)
    for (long long row = this->GetNumRows() - 1; row >= 0; --row) {
        addCrestPoint(row, this->GetNumCols() - 1);
    }
    //left col, except top and bottom, traverse up 
    for (size_t col = this->GetNumCols() - 2; col >= 1; --col){
        addCrestPoint(0, col);
    }
}

const SpokeMesh& RectangularGridSRep::GetUpSpokes() const {
    return this->SkeletonAsMesh.UpSpokes;
}
const SpokeMesh& RectangularGridSRep::GetDownSpokes() const {
    return this->SkeletonAsMesh.DownSpokes;
}
const SpokeMesh& RectangularGridSRep::GetCrestSpokes() const {
    return this->SkeletonAsMesh.CrestSpokes;
}
const std::vector<RectangularGridSRep::UpDownIndices>& RectangularGridSRep::GetCrestSkeletalConnections() const {
    return this->SkeletonAsMesh.CrestSkeletalConnections;
}
const std::vector<RectangularGridSRep::UpDownIndices>& RectangularGridSRep::GetSpine() const {
    return this->SkeletonAsMesh.Spine;
}


std::unique_ptr<RectangularGridSRep> MakeRectangularGridSRep(
    const size_t rows,
    const size_t cols,
    const std::vector<Spoke>& upSpokes,
    const std::vector<Spoke>& downSpokes,
    const std::vector<Spoke>& crestSpokes)
{
    //do a bit of up front validation
    if (upSpokes.size() != downSpokes.size()) {
        throw std::invalid_argument("Must have same number of up and down spokes: "
            + std::to_string(upSpokes.size()) + " != " + std::to_string(downSpokes.size()));
    }
    if (upSpokes.size() != rows * cols) {
        throw std::invalid_argument("Number of up/down spokes must equal the grid size: "
            + std::to_string(upSpokes.size()) + " != " + std::to_string(rows * cols));
    }
    for (size_t i = 0; i < upSpokes.size(); ++i) {
        if (upSpokes[i].GetSkeletalPoint() != downSpokes[i].GetSkeletalPoint()) {
            throw std::invalid_argument("Expecting up and down spokes to be parallel lists at the same skeletal points");
        }
    }
    const auto numCrestPoints = RectangularGridSRep::NumCrestPoints(rows, cols);
    if (crestSpokes.size() != numCrestPoints) {
        throw std::invalid_argument("The number of crest points found does not match what there should be for the grid size: "
            + std::to_string(crestSpokes.size()) + " != " + std::to_string(numCrestPoints));
    }

    if (upSpokes.size() == 0) {
        //just shortcut this here to make things nice further down
        return {};
    }

    RectangularGridSRep::SkeletalGrid grid;
    for (size_t row = 0; row < rows; ++row) {
        std::vector<SkeletalPoint> skeletalPointRow;
        for (size_t col = 0; col < cols; ++col) {
            const size_t index = row * cols + col;
            skeletalPointRow.emplace_back(upSpokes[index], downSpokes[index]);
        }
        grid.emplace_back(std::move(skeletalPointRow));
    }

    // attach crest spokes to their nearest skeletal point. Must not attach to a interior point.
    // if up and down spokes are not the same location by a large margin, this may not work well
    for (const auto& spoke : crestSpokes) {
        double minDist = Point3d::Distance(spoke.GetSkeletalPoint(), grid[0][0].GetUpSpoke().GetSkeletalPoint());
        SkeletalPoint* bestSkeletalPoint = &(grid[0][0]);
        foreachCrestPoint(grid, [&minDist, &bestSkeletalPoint, &spoke](SkeletalPoint& crestPoint) {
            const double dist = Point3d::Distance(spoke.GetSkeletalPoint(), crestPoint.GetUpSpoke().GetSkeletalPoint());
            if (dist < minDist) {
                minDist = dist;
                bestSkeletalPoint = &crestPoint;
            }
        });
        bestSkeletalPoint->SetCrestSpoke(&spoke);
    }

    //all crest spokes should be attached. Verify that every crest point got its crest spoke
    foreachCrestPoint(grid, [](const SkeletalPoint& crestPoint) {
        if (!crestPoint.IsCrest()) {
            throw std::runtime_error("Error attaching crest points");
        }
    });

    return std::unique_ptr<RectangularGridSRep>(new RectangularGridSRep(grid));
}

bool operator==(const RectangularGridSRep& a, const RectangularGridSRep& b) {
    return a.GetSkeletalPoints() == b.GetSkeletalPoints();
}

std::ostream& operator<<(std::ostream& os, const RectangularGridSRep& srep) {
    os << "RectangularGridSRep {" << std::endl
       << "  rows: " << srep.GetNumRows() << std::endl
       << "  cols: " << srep.GetNumCols() << std::endl;
    const auto& grid = srep.GetSkeletalPoints();
    for (size_t i = 0; i < grid.size(); ++i) {
        for (size_t j = 0; j < grid[i].size(); ++j) {
            os << "(" << i << ", " << j << ") " << grid[i][j] << "," << std::endl; 
        }
    }
    os << "}";
    return os;
}

}
