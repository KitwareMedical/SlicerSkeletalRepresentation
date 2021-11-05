#include <srep/RectangularGridSRep.h>
#include <iostream>
#include <limits>

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
{}

RectangularGridSRep::RectangularGridSRep(const SkeletalGrid& skeleton)
    : Skeleton()
{
    this->Validate(skeleton);
    this->Skeleton = skeleton;
}

RectangularGridSRep::RectangularGridSRep(SkeletalGrid&& skeleton)
    : Skeleton()
{
    this->Validate(skeleton);
    this->Skeleton = std::move(skeleton);
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

RectangularGridSRep MakeRectangularGridSRep(const size_t rows,
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
        return RectangularGridSRep();
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

    return RectangularGridSRep(grid);
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
