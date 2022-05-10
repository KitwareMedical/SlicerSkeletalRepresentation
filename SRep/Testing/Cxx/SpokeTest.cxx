#include <gtest/gtest.h>
#include <vtkSRepSpoke.h>
#include <srepUtil.h>
#include <vtkCommand.h>

#include "SRepUnitTestHelpers.h"

using namespace srepUnitTestHelpers;

TEST(SpokeTest, Factories) {
  {
    const auto spoke = vtkSRepSpoke::New();
    auto fin = srep::util::finally([spoke](){ spoke->Delete(); });
    EXPECT_EQ(srep::Point3d(0,0,0), spoke->GetSkeletalPoint());
    EXPECT_EQ(srep::Point3d(0,0,0), spoke->GetBoundaryPoint());
    EXPECT_EQ(srep::Vector3d(0,0,0), spoke->GetDirection());
    EXPECT_EQ(0, spoke->GetRadius());
  }
  {
    const srep::Point3d skeletalPoint(1, 4, 77);
    const srep::Vector3d direction(-2, 55, -0.1);
    const auto boundaryPoint = skeletalPoint + direction;
    const auto radius = direction.GetLength();

    const vtkVector3d vSkeletalPoint(1, 4, 77);
    const vtkVector3d vDirection(-2, 55, -0.1);

    const auto spoke = vtkSRepSpoke::Create(skeletalPoint, direction);
    auto fin = srep::util::finally([spoke](){ spoke->Delete(); });
    const auto smartSRep = vtkSRepSpoke::SmartCreate(skeletalPoint, direction);
    const auto vSRep = vtkSRepSpoke::CreatePointAndDirection(vSkeletalPoint, vDirection);
    auto vfin = srep::util::finally([vSRep](){ vSRep->Delete(); });
    EXPECT_SPOKE_EQ(spoke, smartSRep);
    EXPECT_SPOKE_EQ(spoke, vSRep);

    EXPECT_EQ(skeletalPoint, spoke->GetSkeletalPoint());
    EXPECT_EQ(boundaryPoint, spoke->GetBoundaryPoint());
    EXPECT_EQ(direction, spoke->GetDirection());
    EXPECT_EQ(radius, spoke->GetRadius());
  }
  {
    const srep::Point3d skeletalPoint(0, 0, 0);
    const srep::Point3d boundaryPoint(2, 3, 60);
    const srep::Vector3d direction(2, 3, 60);
    const auto radius = direction.GetLength();

    const vtkVector3d vSkeletalPoint(0, 0, 0);
    const vtkVector3d vBoundaryPoint(2, 3, 60);

    const auto spoke = vtkSRepSpoke::Create(skeletalPoint, boundaryPoint);
    auto fin = srep::util::finally([spoke](){ spoke->Delete(); });
    const auto smartSRep = vtkSRepSpoke::SmartCreate(skeletalPoint, boundaryPoint);
    const auto vSRep = vtkSRepSpoke::CreatePointToPoint(vSkeletalPoint, vBoundaryPoint);
    auto vfin = srep::util::finally([vSRep](){ vSRep->Delete(); });
    EXPECT_SPOKE_EQ(spoke, smartSRep);
    EXPECT_SPOKE_EQ(spoke, vSRep);

    EXPECT_EQ(skeletalPoint, spoke->GetSkeletalPoint());
    EXPECT_EQ(boundaryPoint, spoke->GetBoundaryPoint());
    EXPECT_EQ(direction, spoke->GetDirection());
    EXPECT_EQ(radius, spoke->GetRadius());
  }
}

TEST(SpokeTest, Clone) {
  const auto spoke = vtkSRepSpoke::SmartCreate(srep::Point3d(1, 4, 77), srep::Vector3d(-2, 55, -0.1));

  const auto clone = spoke->Clone();
  auto fin = srep::util::finally([clone](){ clone->Delete(); });
  const auto smartClone = spoke->SmartClone();
  EXPECT_SPOKE_EQ(spoke, clone);
  EXPECT_SPOKE_EQ(spoke, smartClone);
  EXPECT_NE(spoke->GetMTime(), clone->GetMTime());
  EXPECT_NE(spoke->GetMTime(), smartClone->GetMTime());

  // make sure clones are not related to the original
  clone->SetRadius(1234);
  EXPECT_NE(clone->GetRadius(), spoke->GetRadius());

  spoke->SetDirectionOnly(srep::Vector3d(1, 1, 1));
  EXPECT_NE(spoke->GetDirection(), smartClone->GetDirection());
}

TEST(SpokeTest, CloneWithObservers) {
  // observers don't copy with clone
  const auto spoke = vtkSRepSpoke::SmartCreate(srep::Point3d(1, 4, 77), srep::Vector3d(-2, 55, -0.1));

  TestObserver obs;
  const auto tag = spoke->AddObserver(vtkCommand::ModifiedEvent, &obs, &TestObserver::callback);
  auto tagFin = srep::util::finally([tag, spoke](){ spoke->RemoveObserver(tag); });

  EXPECT_EQ(0, obs.numCalls());

  spoke->SetRadius(1);

  EXPECT_EQ(1, obs.numCalls());

  const auto clone = spoke->Clone();
  auto fin = srep::util::finally([clone](){ clone->Delete(); });
  const auto smartClone = spoke->SmartClone();

  clone->SetRadius(9);
  smartClone->SetDirectionOnly(srep::Vector3d(4,5,6));

  EXPECT_EQ(1, obs.numCalls());
}

TEST(SpokeTest, SetGetSRepClasses) {
  auto spoke = vtkSmartPointer<vtkSRepSpoke>::New();

  spoke->SetSkeletalPoint(srep::Point3d(1,2,3));
  EXPECT_EQ(srep::Point3d(1,2,3), spoke->GetSkeletalPoint());

  spoke->SetDirectionAndMagnitude(srep::Vector3d(-1, -2, -3));
  EXPECT_EQ(srep::Vector3d(-1, -2, -3), spoke->GetDirection());
  EXPECT_EQ(srep::Point3d(0,0,0), spoke->GetBoundaryPoint());
  
  spoke->SetDirectionOnly(srep::Vector3d(7, 8, 9));
  const auto expectedDirection = srep::Vector3d(7, 8, 9).Unit() * srep::Vector3d(-1, -2, -3).GetLength();
  EXPECT_EQ(expectedDirection, spoke->GetDirection());
  EXPECT_EQ(srep::Point3d(1,2,3) + expectedDirection, spoke->GetBoundaryPoint());
}

TEST(SpokeTest, SetGetVTKVector3d) {
  auto spoke = vtkSmartPointer<vtkSRepSpoke>::New();
  vtkVector3d output;

  spoke->SetSkeletalPoint(vtkVector3d(1,2,3));
  EXPECT_EQ(srep::Point3d(1,2,3), spoke->GetSkeletalPoint());
  spoke->GetSkeletalPoint(output);
  EXPECT_EQ(vtkVector3d(1,2,3), output);

  spoke->SetDirectionAndMagnitude(vtkVector3d(-1, -2, -3));
  EXPECT_EQ(srep::Vector3d(-1, -2, -3), spoke->GetDirection());
  spoke->GetDirection(output);
  EXPECT_EQ(vtkVector3d(-1, -2, -3), output);

  EXPECT_EQ(srep::Point3d(0,0,0), spoke->GetBoundaryPoint());
  spoke->GetBoundaryPoint(output);
  EXPECT_EQ(vtkVector3d(0,0,0), output);
  
  spoke->SetDirectionOnly(vtkVector3d(7, 8, 9));
  const auto expectedDirection = srep::Vector3d(7, 8, 9).Unit() * srep::Vector3d(-1, -2, -3).GetLength();
  EXPECT_EQ(expectedDirection, spoke->GetDirection());
  spoke->GetDirection(output);
  EXPECT_EQ(vtkVector3d(expectedDirection[0], expectedDirection[1], expectedDirection[2]), output);

  EXPECT_EQ(srep::Point3d(1,2,3) + expectedDirection, spoke->GetBoundaryPoint());
  spoke->GetBoundaryPoint(output);
  EXPECT_EQ(vtkVector3d(1 + expectedDirection[0], 2 + expectedDirection[1], 3 + expectedDirection[2]), output);
}

TEST(SpokeTest, Observation) {
  auto spoke = vtkSRepSpoke::SmartCreate(srep::Point3d(1, 4, 77), srep::Vector3d(-2, 55, -0.1));

  TestObserver obs;
  const auto tag = spoke->AddObserver(vtkCommand::ModifiedEvent, &obs, &TestObserver::callback);
  auto tagFin = srep::util::finally([tag, spoke](){ spoke->RemoveObserver(tag); });

  EXPECT_EQ(0, obs.numCalls());
  spoke->SetRadius(1);
  EXPECT_EQ(1, obs.numCalls());
  spoke->SetSkeletalPoint(srep::Point3d(1,2,3));
  EXPECT_EQ(2, obs.numCalls());
  spoke->SetSkeletalPoint(vtkVector3d(4,5,6));
  EXPECT_EQ(3, obs.numCalls());
  spoke->SetDirectionOnly(srep::Vector3d(5,5,5));
  EXPECT_EQ(4, obs.numCalls());
  spoke->SetDirectionOnly(vtkVector3d(6,6,7));
  EXPECT_EQ(5, obs.numCalls());
  spoke->SetDirectionAndMagnitude(srep::Vector3d(5,5,5));
  EXPECT_EQ(6, obs.numCalls());
  spoke->SetDirectionAndMagnitude(vtkVector3d(6,6,7));
  EXPECT_EQ(7, obs.numCalls());

  // these shouldn't cause the modified
  spoke->GetRadius();

  vtkVector3d v;
  spoke->GetSkeletalPoint();
  spoke->GetSkeletalPoint(v);
  spoke->GetDirection();
  spoke->GetDirection(v);
  spoke->GetBoundaryPoint();
  spoke->GetBoundaryPoint(v);

  std::stringstream ss;
  spoke->PrintSelf(ss, vtkIndent());

  auto clone = vtkSmartPointer<vtkSRepSpoke>::Take(spoke->Clone());
  EXPECT_NE(spoke, clone);
  auto smartClone = spoke->SmartClone();
  EXPECT_NE(spoke, smartClone);

  EXPECT_EQ(7, obs.numCalls());
}

TEST(SpokeTest, PrintSelf) {
  const auto spoke = vtkSRepSpoke::SmartCreate(srep::Point3d(1, 4, 77), srep::Vector3d(-2, 55, -0.1));
  PrintSelfTest(*spoke);
}
