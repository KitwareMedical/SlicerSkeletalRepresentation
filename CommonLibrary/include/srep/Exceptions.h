#ifndef __srep_Exceptions_h
#define __srep_Exceptions_h

#include <stdexcept>

namespace srep {

class InvalidSkeletalGridException : public std::invalid_argument {
    using std::invalid_argument::invalid_argument;
};

}

#endif
