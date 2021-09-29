#include "ReadHeaderFile.h"
#include <srep/Util.h>

#include <algorithm>
#include <fstream>
#include <iterator>

#include <vtk_expat.h>

using namespace srep::util;

namespace srep {
namespace io {
namespace detail {

namespace {


// https://stackoverflow.com/questions/9358718/similar-function-in-c-to-pythons-strip
template <std::ctype_base::mask mask>
class IsNot
{
    std::locale myLocale;       // To ensure lifetime of facet...
    std::ctype<char> const* myCType;
public:
    IsNot( std::locale const& l = std::locale() )
        : myLocale( l )
        , myCType( &std::use_facet<std::ctype<char> >( l ) )
    {
    }
    bool operator()( char ch ) const
    {
        return ! myCType->is( mask, ch );
    }
};

typedef IsNot<std::ctype_base::space> IsNotSpace;

std::string trim( std::string const& original )
{
    std::string::const_iterator right = std::find_if( original.rbegin(), original.rend(), IsNotSpace() ).base();
    std::string::const_iterator left = std::find_if(original.begin(), right, IsNotSpace() );
    return std::string( left, right );
}
// end https://stackoverflow.com/questions/9358718/similar-function-in-c-to-pythons-strip

struct UserData {
    HeaderParams params;
    std::string xmlpath;
    std::string buffer;
};

bool parseBool(const std::string& text) {
    if ("True" == text
        || "true" == text
        || "1" == text)
    {
        return true;
    }
    return false;
}

void startElement(void* vdata, const char* element_cstr, const char** /*attribute*/) {
    if (vdata) {
        UserData* data = reinterpret_cast<UserData*>(vdata);
        data->xmlpath += std::string("/") + element_cstr;
    } else {
        std::cerr << "Bug! empty data in startElement" << std::endl;
    }
}
void endElement(void* vdata, const char* element_cstr) {
    if (vdata) {
        UserData* data = reinterpret_cast<UserData*>(vdata);
        const std::string element(element_cstr);
        if (data->xmlpath.rfind(element) != data->xmlpath.size() - element.size()) {
            std::cerr << "Bug in setting up xmlpath! " << element << " " << data->xmlpath << std::endl;
        }

        //things are much easier without leading/trailing whitespace
        data->buffer = trim(data->buffer);

        if ("/s-rep" == data->xmlpath) {
            //no-op
        }
        else if ("/s-rep/nRows" == data->xmlpath) {
            data->params.nRows = std::stoul(data->buffer);
        }
        else if ("/s-rep/nCols" == data->xmlpath) {
            data->params.nCols = std::stoul(data->buffer);
        } else if ("/s-rep/meshType" == data->xmlpath) {
            if ("Quad" == data->buffer) {
                data->params.meshType = MeshType::Quad;
            } else {
                std::cerr << "Found unexpected meshType: " << data->buffer << std::endl;
            }
        } else if ("/s-rep/color" == data->xmlpath) {
            //no-op
        } else if ("/s-rep/color/red" == data->xmlpath) {
            data->params.color.red = std::stod(data->buffer);
        } else if ("/s-rep/color/green" == data->xmlpath) {
            data->params.color.green = std::stod(data->buffer);
        } else if ("/s-rep/color/blue" == data->xmlpath) {
            data->params.color.blue = std::stod(data->buffer);
        } else if ("/s-rep/isMean" == data->xmlpath) {
            data->params.isMean = parseBool(data->buffer);
        } else if ("/s-rep/meanStatPath" == data->xmlpath) {
            data->params.meanStatPath.swap(data->buffer);
        } else if ("/s-rep/upSpoke" == data->xmlpath) {
            data->params.upSpoke.swap(data->buffer);
        } else if ("/s-rep/downSpoke" == data->xmlpath) {
            data->params.downSpoke.swap(data->buffer);
        } else if ("/s-rep/crestSpoke" == data->xmlpath) {
            data->params.crestSpoke.swap(data->buffer);
        } else {
            std::cerr << "Found unexpected element! " << data->xmlpath << std::endl;
        }

        //peel this element off the path since we are ending it
        data->xmlpath.resize(data->xmlpath.size() - element.size() - 1); //extra -1 for the trailing /
        data->buffer.clear();
    } else {
        std::cerr << "Bug! empty data in endElement" << std::endl;
    }
}

void handleData(void *vdata, const char *content, int length) {
    if (vdata) {
        UserData* data = reinterpret_cast<UserData*>(vdata);
        data->buffer.insert(data->buffer.end(), content, content + length);
    } else {
        std::cerr << "Bug! empty data in handleData" << std::endl;
    }
}

} // namespace {}

HeaderParams readHeaderFile(const std::string& filename) {
    UserData userData;
    userData.buffer.reserve(256);
    XML_Parser parser = XML_ParserCreate(nullptr);
    const auto fin = finally([&parser](){
        XML_ParserFree(parser);
    });
    XML_SetUserData(parser, &userData);
    XML_SetElementHandler(parser, startElement, endElement);
    XML_SetCharacterDataHandler(parser, handleData);

    //these headers aren't very big, so just read and parse it all at once
    std::ifstream in(filename);
    const std::string contents((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());

    if (XML_Parse(parser, contents.c_str(), contents.size(), 1) == 0) {
        printf("Error: %s\n", XML_ErrorString(XML_GetErrorCode(parser)));
    }

    return userData.params;
}

std::ostream& operator<<(std::ostream& os, const Color& color) {
    os << "(" << color.red << ", " << color.green << ", " << color.blue << ")";
    return os;
}

std::ostream& operator<<(std::ostream& os, const HeaderParams& params) {
    os << "nRows: " << params.nRows << std::endl
       << "nCols: " << params.nCols << std::endl
       << "meshType: " << static_cast<int>(params.meshType) << std::endl
       << "color: " << params.color << std::endl
       << "isMean: " << (params.isMean ? std::string("true") : std::string("false")) << std::endl
       << "meanStatPath: " << params.meanStatPath << std::endl
       << "upSpoke: " << params.upSpoke << std::endl
       << "downSpoke: " << params.downSpoke << std::endl
       << "crestSpoke: " << params.crestSpoke << std::endl
    ;
    return os;
}

}
}
}
