#ifndef __srep_ReadHeaderFile_h
#define __srep_ReadHeaderFile_h

#include <iostream>
#include <string>

namespace srep {
namespace io {
namespace detail {

enum class MeshType {
    Quad
};

struct Color {
    double red;
    double green;
    double blue;
};
std::ostream& operator<<(std::ostream& os, const Color& color);

struct HeaderParams {
    size_t nRows;
    size_t nCols;
    MeshType meshType;
    Color color;
    bool isMean;
    std::string meanStatPath;
    std::string upSpoke;
    std::string downSpoke;
    std::string crestSpoke;
};
std::ostream& operator<<(std::ostream& os, const HeaderParams& params);

HeaderParams readHeaderFile(const std::string& filename);

}
}
}

#endif
