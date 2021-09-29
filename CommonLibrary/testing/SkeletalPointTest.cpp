#include <gtest/gtest.h>
#include <srep/SkeletalPoint.h>

using namespace srep;

TEST(SkeletalPoint, test) {
    {
        //default constructor
        SkeletalPoint s;
        EXPECT_EQ(s.GetPoint(), Point3d{});
        EXPECT_EQ(s.GetUpSpoke(), Spoke{});
        EXPECT_EQ(s.GetDownSpoke(), Spoke{});
        EXPECT_FALSE(s.IsCrest());
        EXPECT_THROW(s.GetCrestSpoke(), NotACrestException);
    }
    {
        //not a crest constructor
        const Point3d point(1,2,3);
        const Spoke up(point, Vector3d(0,1,0));
        const Spoke down(point, Vector3d(0,-1,0));
        const SkeletalPoint s(up, down);
        EXPECT_EQ(s.GetPoint(), point);
        EXPECT_EQ(s.GetUpSpoke(), up);
        EXPECT_EQ(s.GetDownSpoke(), down);
        EXPECT_FALSE(s.IsCrest());
    }
    {
        //crest constructor
        //not a crest constructor
        const Point3d point(1,2,3);
        const Spoke up(point, Vector3d(0,1,0));
        const Spoke down(point, Vector3d(0,-1,0));
        const Spoke crest(point, Vector3d(1,0,0));
        const SkeletalPoint s(up, down, crest);
        EXPECT_EQ(s.GetPoint(), point);
        EXPECT_EQ(s.GetUpSpoke(), up);
        EXPECT_EQ(s.GetDownSpoke(), down);
        EXPECT_TRUE(s.IsCrest());
        EXPECT_EQ(s.GetCrestSpoke(), crest);
    }
    {
        //misaligned spokes non-crest
        const Spoke up(Point3d(0,0,0), Vector3d(0,1,0));
        const Spoke down(Point3d(0,0,1), Vector3d(0,-1,0));
        EXPECT_THROW(SkeletalPoint s(up, down), MisalignedSpokeException);
    }
    {
        //misaligned spokes crest
        const Spoke up(Point3d(0,0,0), Vector3d(0,1,0));
        const Spoke down(Point3d(0,0,0), Vector3d(0,-1,0));
        const Spoke crest(Point3d(0.1,0,0), Vector3d(1,0,0));
        EXPECT_THROW(SkeletalPoint s(up, down, crest), MisalignedSpokeException);
    }
}