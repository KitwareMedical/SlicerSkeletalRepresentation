#include <gtest/gtest.h>
#include <srepVector3d.h>

using namespace srep;

constexpr double tol = 1e-5;

#define EXPECT_VEC_NEAR(a, b) \
    EXPECT_NEAR((a).GetX(), (b).GetX(), tol); \
    EXPECT_NEAR((a).GetY(), (b).GetY(), tol); \
    EXPECT_NEAR((a).GetZ(), (b).GetZ(), tol)

TEST(Vector3d, defaultConstructor) {
    const Vector3d p1;
    const Vector3d p2;
    EXPECT_NEAR(p1.GetX(), 0.0, tol);
    EXPECT_NEAR(p1.GetY(), 0.0, tol);
    EXPECT_NEAR(p1.GetZ(), 0.0, tol);
    EXPECT_EQ(p1, p2);
}

TEST(Vector3d, threeDoubleConstructor) {
    const Vector3d p1(-1.234, 3.21, 9.99);
    const Vector3d p2(-1.234, 3.21, 9.99);
    EXPECT_NEAR(p1.GetX(), -1.234, tol);
    EXPECT_NEAR(p1.GetY(), 3.21, tol);
    EXPECT_NEAR(p1.GetZ(), 9.99, tol);
    EXPECT_EQ(p1, p2);
}

TEST(Vector3d, doubleArrayConstructor) {
    const double p[3] = {12.234, -31.21, -9.99};
    const Vector3d p1(p);
    const Vector3d p2(p);
    EXPECT_NEAR(p1.GetX(), 12.234, tol);
    EXPECT_NEAR(p1.GetY(), -31.21, tol);
    EXPECT_NEAR(p1.GetZ(), -9.99, tol);
    EXPECT_EQ(p1, p2);
}

TEST(Vector3d, sets) {
    Vector3d p;

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

TEST(Vector3d, copyAndCopyAssign) {
    const Vector3d p1(1,2,3);
    Vector3d p2(p1);
    EXPECT_EQ(p1, p2);

    Vector3d p3 = p1;
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

TEST(Vector3d, moveAndMoveAssign) {
    Vector3d p1(1,2,3);

    EXPECT_NEAR(p1.GetX(), 1, tol);
    EXPECT_NEAR(p1.GetY(), 2, tol);
    EXPECT_NEAR(p1.GetZ(), 3, tol);

    Vector3d p2(std::move(p1));

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

    Vector3d p3 = std::move(p1);

    EXPECT_NEAR(p3.GetX(), 4, tol);
    EXPECT_NEAR(p3.GetY(), 5, tol);
    EXPECT_NEAR(p3.GetZ(), 6, tol);
}

TEST(Vector3d, relationalOperators) {
    const Vector3d p1(1,2,3);
    const Vector3d p1b(p1);
    const Vector3d p2(2,2,3);
    const Vector3d p3(1,3,3);
    const Vector3d p4(1,2,4);

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

TEST(Vector3d, PointConstructor) {
    const Point3d p1(1,2,3);
    const Point3d p2(1,2,4);

    const Vector3d v(p1, p2);
    const Vector3d truth(0,0,1);

    EXPECT_VEC_NEAR(truth, v);
}

TEST(Vector3d, GetLength) {
    {
        const Vector3d v(5,5,15);
        EXPECT_NEAR(v.GetLength(), 16.583123951777, tol);
    }
    {
        const Vector3d v(-5,-5,15);
        EXPECT_NEAR(v.GetLength(), 16.583123951777, tol);
    }
    {
        const Vector3d v(0,0,0);
        EXPECT_NEAR(v.GetLength(), 0, tol);
    }
    {
        const Vector3d v(1,2,3);
        EXPECT_NEAR(v.GetLength(), 3.7416573867739413, tol);
    }
    {
        constexpr double length = 16.583123951777;
        const Vector3d v(5/length, 5/length, 15/length);
        EXPECT_NEAR(v.GetLength(), 1, tol);
    }
}

TEST(Vector3d, Resize) {
    {
        Vector3d v(0,0,0);
        EXPECT_THROW(v.Resize(1), std::runtime_error);
        EXPECT_THROW(Vector3d::Resize(v, 1), std::runtime_error);
    }
    {
        Vector3d v(1,-1,1);
        constexpr double newLength = 3;
        const Vector3d truth(1.73205080757, -1.73205080757, 1.73205080757); //3 / sqrt(3)
        EXPECT_VEC_NEAR(truth, Vector3d::Resize(v, newLength));
        v.Resize(newLength);
        EXPECT_VEC_NEAR(truth, v);
    }
    {
        Vector3d v(1,2,3);
        constexpr double newLength = 0;
        const Vector3d truth(0, 0, 0);
        EXPECT_VEC_NEAR(truth, Vector3d::Resize(v, newLength));
        v.Resize(newLength);
        EXPECT_VEC_NEAR(truth, v);
    }
    {
        Vector3d v(3,4,0);
        constexpr double newLength = 10;
        const Vector3d truth(6,8,0);
        EXPECT_VEC_NEAR(truth, Vector3d::Resize(v, newLength));
        v.Resize(newLength);
        EXPECT_VEC_NEAR(truth, v);
    }
}

TEST(Vector3d, Unit) {
    {
        constexpr double length = 16.583123951777;
        const Vector3d v(5,5,15);

        const Vector3d truth(5/length, 5/length, 15/length);
        EXPECT_VEC_NEAR(truth, v.Unit());
    }
    {
        // test when we are a unit vector
        constexpr double length = 16.583123951777;
        const Vector3d v(5/length, 5/length, 15/length);

        EXPECT_VEC_NEAR(v, v.Unit());
    }
    {
        constexpr double length = 948.6854062332782;
        const Vector3d v(900, -2, 300);

        const Vector3d truth(900/length, -2/length, 300/length);
        EXPECT_VEC_NEAR(truth, v.Unit());
    }
    {
        //test when length is 0
        const Vector3d v(0,0,0);
        EXPECT_THROW(v.Unit(), std::runtime_error);
    }
}

TEST(Vector3d, add) {
    {
        const Vector3d v1(1,2,3);
        const Vector3d v2(4,5,6);
        const Vector3d truth(5,7,9);
        EXPECT_VEC_NEAR(truth, v1 + v2);
    }
    {
        const Vector3d v1(-1,2,3);
        const Vector3d v2(4,-5,6);
        const Vector3d truth(3,-3,9);
        EXPECT_VEC_NEAR(truth, v1 + v2);
    }
    {
        const Vector3d v1(1,2,3);
        const Vector3d v2;
        const Vector3d truth(1,2,3);
        EXPECT_VEC_NEAR(truth, v1 + v2);
    }
}

TEST(Vector3d, subtract) {
    {
        const Vector3d v1(1,2,3);
        const Vector3d v2(4,5,6);
        const Vector3d truth(-3,-3,-3);
        EXPECT_VEC_NEAR(truth, v1 - v2);
    }
    {
        const Vector3d v1(-1,2,3);
        const Vector3d v2(4,-5,6);
        const Vector3d truth(-5,7,-3);
        EXPECT_VEC_NEAR(truth, v1 - v2);
    }
    {
        const Vector3d v1(1,2,3);
        const Vector3d v2;
        const Vector3d truth(1,2,3);
        EXPECT_VEC_NEAR(truth, v1 - v2);
    }
}

TEST(Vector3d, multiply) {
    {
        const Vector3d v(1,2,3);
        const double multiplier = 7;
        const Vector3d truth(7,14,21);
        EXPECT_VEC_NEAR(truth, v * multiplier);
    }
    {
        const Vector3d v(-1,2,3);
        const double multiplier = -7;
        const Vector3d truth(7,-14,-21);
        EXPECT_VEC_NEAR(truth, v * multiplier);
    }
    {
        const Vector3d v(1,2,3);
        const double multiplier = 0;
        const Vector3d truth(0,0,0);
        EXPECT_VEC_NEAR(truth, v * multiplier);
    }
    {
        const Vector3d v(1,2,3);
        const double multiplier = 0.5;
        const Vector3d truth(0.5, 1.0, 1.5);
        EXPECT_VEC_NEAR(truth, v * multiplier);
    }
    {
        const Vector3d v(1,2,3);
        const double multiplier = 1;
        const Vector3d truth(1,2,3);
        EXPECT_VEC_NEAR(truth, v * multiplier);
    }
}

TEST(Vector3d, divide) {
    {
        const Vector3d v(1,2,3);
        const double divisor = 2;
        const Vector3d truth(0.5, 1.0, 1.5);
        EXPECT_VEC_NEAR(truth, v / divisor);
    }
    {
        const Vector3d v(-1,2,3);
        const double divisor = -5;
        const Vector3d truth(0.2,-0.4,-0.6);
        EXPECT_VEC_NEAR(truth, v / divisor);
    }
    {
        const Vector3d v(1,2,3);
        const double divisor = 0.5;
        const Vector3d truth(2,4,6);
        EXPECT_VEC_NEAR(truth, v / divisor);
    }
    {
        const Vector3d v(1,2,3);
        const double divisor = 1;
        const Vector3d truth(1,2,3);
        EXPECT_VEC_NEAR(truth, v / divisor);
    }
}

TEST(Vector3d, addPoint) {
    {
        const Point3d p(1,2,3);
        const Vector3d v(4,5,6);
        const Point3d truth(5,7,9);
        EXPECT_VEC_NEAR(truth, p + v);
    }
    {
        const Point3d p(-1,2,3);
        const Vector3d v(4,-5,6);
        const Point3d truth(3,-3,9);
        EXPECT_VEC_NEAR(truth, p + v);
    }
    {
        const Vector3d p(1,2,3);
        const Vector3d v;
        const Vector3d truth(1,2,3);
        EXPECT_VEC_NEAR(truth, p + v);
    }
}

TEST(Vector3d, subtractPoint) {
    {
        const Point3d p(1,2,3);
        const Vector3d v(4,5,6);
        const Point3d truth(-3,-3,-3);
        EXPECT_VEC_NEAR(truth, p - v);
    }
    {
        const Point3d p(-1,2,3);
        const Vector3d v(4,-5,6);
        const Point3d truth(-5,7,-3);
        EXPECT_VEC_NEAR(truth, p - v);
    }
    {
        const Point3d p(1,2,3);
        const Vector3d v;
        const Point3d truth(1,2,3);
        EXPECT_VEC_NEAR(truth, p - v);
    }
}
