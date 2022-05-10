#include <gtest/gtest.h>
#include <srepPoint3d.h>

using namespace srep;

constexpr double tol = 1e-5;

TEST(Point3d, defaultConstructor) {
    const Point3d p1;
    const Point3d p2;
    EXPECT_NEAR(p1.GetX(), 0.0, tol);
    EXPECT_NEAR(p1.GetY(), 0.0, tol);
    EXPECT_NEAR(p1.GetZ(), 0.0, tol);
    EXPECT_EQ(p1, p2);
}

TEST(Point3d, threeDoubleConstructor) {
    const Point3d p1(-1.234, 3.21, 9.99);
    const Point3d p2(-1.234, 3.21, 9.99);
    EXPECT_NEAR(p1.GetX(), -1.234, tol);
    EXPECT_NEAR(p1.GetY(), 3.21, tol);
    EXPECT_NEAR(p1.GetZ(), 9.99, tol);
    EXPECT_EQ(p1, p2);
}

TEST(Point3d, doubleArrayConstructor) {
    const double p[3] = {12.234, -31.21, -9.99};
    const Point3d p1(p);
    const Point3d p2(p);
    EXPECT_NEAR(p1.GetX(), 12.234, tol);
    EXPECT_NEAR(p1.GetY(), -31.21, tol);
    EXPECT_NEAR(p1.GetZ(), -9.99, tol);
    EXPECT_EQ(p1, p2);
}

TEST(Point3d, sets) {
    Point3d p;

    EXPECT_NEAR(p.GetX(), 0.0, tol);
    EXPECT_NEAR(p.GetY(), 0.0, tol);
    EXPECT_NEAR(p.GetZ(), 0.0, tol);

    p.SetX(5.55);

    EXPECT_NEAR(p.GetX(), 5.55, tol);
    EXPECT_NEAR(p.GetY(), 0.0, tol);
    EXPECT_NEAR(p.GetZ(), 0.0, tol);

    p.SetZ(2.22);

    EXPECT_NEAR(p.GetX(), 5.55, tol);
    EXPECT_NEAR(p.GetY(), 0.0, tol);
    EXPECT_NEAR(p.GetZ(), 2.22, tol);

    p.SetY(-0.0001);

    EXPECT_NEAR(p.GetX(), 5.55, tol);
    EXPECT_NEAR(p.GetY(), -0.0001, tol);
    EXPECT_NEAR(p.GetZ(), 2.22, tol);
}

TEST(Point3d, copyAndCopyAssign) {
    const Point3d p1(1,2,3);
    Point3d p2(p1);
    EXPECT_EQ(p1, p2);

    Point3d p3 = p1;
    EXPECT_EQ(p1, p3);
    EXPECT_EQ(p2, p3);

    //deep copy
    p3.SetX(0);

    EXPECT_EQ(p1, p2);
    EXPECT_NE(p1, p3);

    p2.SetY(0);

    EXPECT_NE(p1, p2);
    EXPECT_NE(p1, p3);
    EXPECT_NE(p2, p3);
}

TEST(Point3d, moveAndMoveAssign) {
    Point3d p1(1,2,3);

    EXPECT_NEAR(p1.GetX(), 1, tol);
    EXPECT_NEAR(p1.GetY(), 2, tol);
    EXPECT_NEAR(p1.GetZ(), 3, tol);

    Point3d p2(std::move(p1));

    EXPECT_NEAR(p2.GetX(), 1, tol);
    EXPECT_NEAR(p2.GetY(), 2, tol);
    EXPECT_NEAR(p2.GetZ(), 3, tol);

    //test p1 in unspecified but valid state (i.e. can still be used)
    p1.SetX(4);
    p1.SetY(5);
    p1.SetZ(6);

    EXPECT_NEAR(p1.GetX(), 4, tol);
    EXPECT_NEAR(p1.GetY(), 5, tol);
    EXPECT_NEAR(p1.GetZ(), 6, tol);

    Point3d p3 = std::move(p1);

    EXPECT_NEAR(p3.GetX(), 4, tol);
    EXPECT_NEAR(p3.GetY(), 5, tol);
    EXPECT_NEAR(p3.GetZ(), 6, tol);
}

TEST(Point3d, relationalOperators) {
    const Point3d p1(1,2,3);
    const Point3d p1b(p1);
    const Point3d p2(2,2,3);
    const Point3d p3(1,3,3);
    const Point3d p4(1,2,4);

    // ==
    EXPECT_TRUE (p1 == p1b);
    EXPECT_FALSE(p1 == p2);
    EXPECT_FALSE(p1 == p3);
    EXPECT_FALSE(p1 == p4);

    // !=
    EXPECT_FALSE(p1 != p1b);
    EXPECT_TRUE (p1 != p2);
    EXPECT_TRUE (p1 != p3);
    EXPECT_TRUE (p1 != p4);

    // <
    EXPECT_FALSE(p1 < p1b);
    EXPECT_TRUE (p1 < p2);
    EXPECT_TRUE (p1 < p3);
    EXPECT_TRUE (p1 < p4);
    EXPECT_FALSE(p2 < p1);
    EXPECT_FALSE(p3 < p1);
    EXPECT_FALSE(p4 < p1);

    // >
    EXPECT_FALSE(p1 > p1b);
    EXPECT_FALSE(p1 > p2);
    EXPECT_FALSE(p1 > p3);
    EXPECT_FALSE(p1 > p4);
    EXPECT_TRUE (p2 > p1);
    EXPECT_TRUE (p3 > p1);
    EXPECT_TRUE (p4 > p1);

    // <=
    EXPECT_TRUE (p1 <= p1b);
    EXPECT_TRUE (p1 <= p2);
    EXPECT_TRUE (p1 <= p3);
    EXPECT_TRUE (p1 <= p4);
    EXPECT_FALSE(p2 <= p1);
    EXPECT_FALSE(p3 <= p1);
    EXPECT_FALSE(p4 <= p1);

    // >=
    EXPECT_TRUE(p1 >= p1b);
    EXPECT_FALSE(p1 >= p2);
    EXPECT_FALSE(p1 >= p3);
    EXPECT_FALSE(p1 >= p4);
    EXPECT_TRUE (p2 >= p1);
    EXPECT_TRUE (p3 >= p1);
    EXPECT_TRUE (p4 >= p1);
}
