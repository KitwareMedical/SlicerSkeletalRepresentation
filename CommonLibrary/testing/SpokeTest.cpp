#include <gtest/gtest.h>
#include <srep/Spoke.h>

using namespace srep;

TEST(Spoke, constructors) {
    {
        const Spoke spoke;
        const Point3d  expectedSkeletalPoint(0,0,0);
        const Vector3d expectedDirection(0,0,0);
        const Point3d  expectedBoundaryPoint(0,0,0);
        EXPECT_EQ(spoke.GetRadius(), 0);
        EXPECT_EQ(spoke.GetSkeletalPoint(), expectedSkeletalPoint);
        EXPECT_EQ(spoke.GetDirection(), expectedDirection);
        EXPECT_EQ(spoke.GetBoundaryPoint(), expectedBoundaryPoint);
    }
    {
        const Point3d  skeletalPoint(-5,10,20);
        const Vector3d direction(3,4,-5);
        const Point3d  expectedBoundaryPoint(-2, 14, 15);
        const Spoke spoke(skeletalPoint, direction);
        EXPECT_EQ(spoke.GetRadius(), direction.GetLength());
        EXPECT_EQ(spoke.GetSkeletalPoint(), skeletalPoint);
        EXPECT_EQ(spoke.GetDirection(), direction);
        EXPECT_EQ(spoke.GetBoundaryPoint(), expectedBoundaryPoint);
    }
}

TEST(Spoke, GetSetSkeletalPoint) {
    {
        //test const correctness
        const Point3d skeletalPoint(5, -9, 2.5);
        const Vector3d direction(0,0,0);
        const Spoke spoke(skeletalPoint, direction);
        EXPECT_EQ(spoke.GetSkeletalPoint(), skeletalPoint);
    }
    {
        const Point3d skeletalPoint(5, -9, 2.5);
        const Vector3d direction(0,0,0);
        Spoke spoke(skeletalPoint, direction);
        EXPECT_EQ(spoke.GetSkeletalPoint(), skeletalPoint);

        const Point3d newPoint1(-0.1, 99, 23);
        spoke.SetSkeletalPoint(newPoint1);
        EXPECT_EQ(spoke.GetSkeletalPoint(), newPoint1);

        const Point3d newPoint2(0.1, -99, -23);
        spoke.SetSkeletalPoint(newPoint2);
        EXPECT_EQ(spoke.GetSkeletalPoint(), newPoint2);
    }
}

TEST(Spoke, GetSetDirection) {
    {
        //test const correctness
        const Point3d skeletalPoint(0,0,0);
        const Vector3d direction(3,2,1);
        const Spoke spoke(skeletalPoint, direction);
        EXPECT_EQ(spoke.GetDirection(), direction);
    }
    {
        //test const correctness
        const Point3d skeletalPoint(0,0,0);
        const Vector3d direction(3,-2,1);
        Spoke spoke(skeletalPoint, direction);
        EXPECT_EQ(spoke.GetDirection(), direction);

        const Vector3d newDirection1(5,10,15);
        spoke.SetDirectionAndMagnitude(newDirection1);
        EXPECT_EQ(spoke.GetDirection(), newDirection1);

        const Vector3d newDirection2(-3,21,0);
        spoke.SetDirectionOnly(newDirection2);
        const Vector3d expectedDirection = Vector3d::Resize(newDirection2, newDirection1.GetLength());
        EXPECT_EQ(spoke.GetDirection(), expectedDirection);
    }
}

TEST(Spoke, Equality) {
    const Point3d skeletalPoint(0.12345, 0.0009, 8970.);
    const Vector3d direction(3.865, 2.872, 1.984686);

    const Spoke spoke1(skeletalPoint, direction);
    Spoke spoke2(skeletalPoint, direction);

    EXPECT_EQ(spoke1, spoke2);
}