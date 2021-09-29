#include <srep/Spoke.h>
#include <srep/SRepIO.h>
#include <srep/Util.h>
#include "ReadHeaderFile.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkXMLPolyDataReader.h>

using namespace srep::io::detail;

namespace srep {
namespace io {

namespace {

const char * const RadiusArrayName = "spokeLength";
const char * const DirectionArrayName = "spokeDirection";

std::vector<Spoke> readSpokeFile(const std::string& filename) {
    if (filename.empty()) {
        throw std::invalid_argument("readSpokeFile requires a filename");
    }

    std::vector<Spoke> spokes;

    vtkNew<vtkXMLPolyDataReader> reader;
    reader->SetFileName(filename.c_str());
    reader->Update();

    auto* spokeData = reader->GetOutput();
    auto* pointData = spokeData->GetPointData();
    const auto numberOfArrays = pointData->GetNumberOfArrays();

    if (numberOfArrays == 0) {
        std::cerr << "File does not contain data: " << filename << std::endl;
    }

    const auto numberOfPoints = spokeData->GetNumberOfPoints();
    
    auto* radiusArray = pointData->GetArray(RadiusArrayName);
    if (radiusArray->GetNumberOfComponents() != 1) {
        throw std::runtime_error("Malformed srep file, expected 1 component in radius array: " + filename);
    }
    auto* directionArray = pointData->GetArray(DirectionArrayName);
    if (directionArray->GetNumberOfComponents() != 3) {
        throw std::runtime_error("Malformed srep file, expected 3 components in radius array: " + filename);
    }
    for (int i = 0 ; i < numberOfPoints; ++i) {
        double skeletalPoint[3];
        spokeData->GetPoint(i, skeletalPoint);

        const double* direction = directionArray->GetTuple3(i);
        const double radius = radiusArray->GetTuple1(i);

        spokes.emplace_back(Point3d(skeletalPoint),
            Vector3d(direction[0] * radius, direction[1] * radius, direction[2] * radius));
    }

    return spokes;
}

} // namespace {}

SRep ReadSRep(const std::string& filename) {
    const auto headerParams = readHeaderFile(filename);
    const auto upSpokes = readSpokeFile(headerParams.upSpoke);
    const auto downSpokes = readSpokeFile(headerParams.downSpoke);
    const auto crestSpokes = readSpokeFile(headerParams.crestSpoke);

    //TODO: use the rest of the headerParams

    return MakeSRep(headerParams.nRows, headerParams.nCols, upSpokes, downSpokes, crestSpokes);
}

void WriteSRep(const std::string& /*filename*/, const SRep& /*srep*/) {
    //TODO: implement
}

}
}
