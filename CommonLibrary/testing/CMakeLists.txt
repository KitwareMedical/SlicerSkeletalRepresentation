find_package(GTest REQUIRED CONFIG)

add_executable(SRepTest
    Point3dTest.cpp
    SkeletalPointTest.cpp
    SpokeTest.cpp
    SRepIOTest.cpp
    Vector3dTest.cpp
)
target_link_libraries(SRepTest
PRIVATE
    srep::srep
    GTest::gtest GTest::gmock GTest::gmock_main
)
