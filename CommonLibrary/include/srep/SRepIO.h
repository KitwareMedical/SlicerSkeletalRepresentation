#ifndef __srep_SRepIO_h
#define __srep_SRepIO_h

#include <string>
#include <srep/SRep.h>

namespace srep {
namespace io {

/// Reads SRep from XML file.
SRep ReadSRep(const std::string& filename);

/// Writes SRep to XML file.
void WriteSRep(const std::string& filename, const SRep& srep);

}
}

#endif
