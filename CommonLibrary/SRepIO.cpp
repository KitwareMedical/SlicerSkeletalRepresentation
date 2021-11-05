#include <srep/Spoke.h>
#include <srep/SRepIO.h>
#include <srep/Util.h>
#include "ReadHeaderFile.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <vtkDoubleArray.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>

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

void writeHeaderFile(const RectangularGridSRep& srep,
                     const std::string& headerFilename,
                     const std::string& upFilename,
                     const std::string& downFilename,
                     const std::string& crestFilename)
{
    std::ofstream out(headerFilename);
    out << "<s-rep>" << std::endl
        << "  <nRows>" << srep.GetNumRows() << "</nRows>" << std::endl
        << "  <nCols>" << srep.GetNumCols() << "</nCols>" << std::endl
        << "  <meshType>Quad</meshType>" << std::endl
        << "  <upSpoke>" << upFilename << "</upSpoke>" << std::endl
        << "  <downSpoke>" << downFilename << "</downSpoke>" << std::endl
        << "  <crestSpoke>" << crestFilename << "</crestSpoke>" << std::endl
        << "</s-rep>" << std::endl
    ;
}

void writeSpokeFiles(const RectangularGridSRep& srep,
                     const std::string& upFilename,
                     const std::string& downFilename,
                     const std::string& crestFilename)
{
    vtkNew<vtkPolyData> upPolyData;
    vtkNew<vtkPolyData> downPolyData;
    vtkNew<vtkPolyData> crestPolyData;

    const auto initPoly = [](vtkPolyData* poly) {
        vtkSmartPointer<vtkDoubleArray> spokeDirection = vtkSmartPointer<vtkDoubleArray>::New();
        spokeDirection->SetName(DirectionArrayName);
        spokeDirection->SetNumberOfComponents(3);

        vtkSmartPointer<vtkDoubleArray> spokeLengths = vtkSmartPointer<vtkDoubleArray>::New();
        spokeLengths->SetName(RadiusArrayName);
        spokeLengths->SetNumberOfComponents(1);

        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

        poly->SetPoints(points);
        poly->GetPointData()->AddArray(spokeDirection);
        poly->GetPointData()->SetActiveVectors(DirectionArrayName);
        poly->GetPointData()->AddArray(spokeLengths);
        poly->GetPointData()->SetActiveScalars(RadiusArrayName);
    };

    initPoly(upPolyData);
    initPoly(downPolyData);
    initPoly(crestPolyData);

    const auto insertSpoke = [](vtkPolyData* poly, const Spoke& spoke) {
        poly->GetPoints()->InsertNextPoint(spoke.GetSkeletalPoint().AsArray().data());
        poly->GetPointData()->GetArray(DirectionArrayName)->InsertNextTuple(spoke.GetDirection().Unit().AsArray().data());
        poly->GetPointData()->GetArray(RadiusArrayName)->InsertNextTuple1(spoke.GetDirection().GetLength());
    };

    //relying on fact foreachPoint goes in row major order
    foreachPoint(srep, [&](const SkeletalPoint& point) {
        insertSpoke(upPolyData, point.GetUpSpoke());
        insertSpoke(downPolyData, point.GetDownSpoke());
        if (point.IsCrest()) {
            insertSpoke(crestPolyData, point.GetCrestSpoke());
        }
    });

    vtkNew<vtkXMLPolyDataWriter> writer;
    writer->SetDataModeToAscii();
    const auto writePoly = [&writer](const std::string& filename, vtkPolyData* poly) {
        writer->SetFileName(filename.c_str());
        writer->SetInputData(poly);
        writer->Update();
    };
    writePoly(upFilename, upPolyData);
    writePoly(downFilename, downPolyData);
    writePoly(crestFilename, crestPolyData);
}

} // namespace {}

RectangularGridSRep ReadRectangularGridSRep(const std::string& filename) {
    const auto headerParams = readHeaderFile(filename);
    const auto upSpokes = readSpokeFile(headerParams.upSpoke);
    const auto downSpokes = readSpokeFile(headerParams.downSpoke);
    const auto crestSpokes = readSpokeFile(headerParams.crestSpoke);

    //TODO: use the rest of the headerParams

    return MakeRectangularGridSRep(headerParams.nRows, headerParams.nCols, upSpokes, downSpokes, crestSpokes);
}

void WriteSRep(const RectangularGridSRep& srep,
               const std::string& headerFilename,
               const std::string& upFilename,
               const std::string& downFilename,
               const std::string& crestFilename)
{
    writeHeaderFile(srep, headerFilename, upFilename, downFilename, crestFilename);
    writeSpokeFiles(srep, upFilename, downFilename, crestFilename);
}

}
}
