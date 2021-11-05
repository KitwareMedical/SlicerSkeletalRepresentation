#ifndef __srep_SRepIO_h
#define __srep_SRepIO_h

#include <string>
#include <srep/RectangularGridSRep.h>

namespace srep {
namespace io {

/// Reads RectangularGridSRep from XML file.
RectangularGridSRep ReadRectangularGridSRep(const std::string& filename);

/// Writes RectangularGridSRep to XML file.
void WriteSRep(const RectangularGridSRep& srep,
               const std::string& headerFilename,
               const std::string& upFilename,
               const std::string& downFilename,
               const std::string& crestFilename);

}
}

#endif
