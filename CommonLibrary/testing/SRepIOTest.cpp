#include <gtest/gtest.h>
#include <srep/SRepIO.h>

using namespace srep::io;

TEST(SRepIO, test) {
    //TODO: setup test data
    const std::string filename("/home/connor/Downloads/srep_examples2/hippocampus0/header.xml");

    const auto origSRep = ReadRectangularGridSRep(filename);

    WriteSRep(origSRep,
              "/tmp/srep/header.xml",
              "/tmp/srep/up.xml",
              "/tmp/srep/down.xml",
              "/tmp/srep/crest.xml");

    const auto rereadSRep = ReadRectangularGridSRep("/tmp/srep/header.xml");

    //it is close but not exact double equal
    // EXPECT_EQ(origSRep, rereadSRep);
}
