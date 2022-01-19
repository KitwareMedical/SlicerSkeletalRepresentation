#include <gtest/gtest.h>
#include <vtkSRepSkeletalPoint.h>
#include <srepUtil.h>
#include <vtkCommand.h>

#include "SRepUnitTestHelpers.h"

using namespace srepUnitTestHelpers;

TEST(SkeletalPointTest, Factories) {
  {
    const auto emptySpoke = vtkSmartPointer<vtkSRepSpoke>::New();
    const auto pt = vtkSRepSkeletalPoint::New();
    auto fin = srep::util::finally([pt](){ pt->Delete(); });
    EXPECT_SPOKE_EQ(emptySpoke, pt->GetUpSpoke());
    EXPECT_SPOKE_EQ(emptySpoke, pt->GetDownSpoke());
    EXPECT_EQ(nullptr, pt->GetCrestSpoke());
    EXPECT_FALSE(pt->IsCrest());
  }
  {
    auto upSpoke = vtkSRepSpoke::SmartCreate(srep::Point3d(0,0,0), srep::Point3d(1,1,1));
    auto downSpoke = vtkSRepSpoke::SmartCreate(srep::Point3d(0,0,0), srep::Point3d(-1,-1,-1));

    auto pt = vtkSRepSkeletalPoint::Create(upSpoke, downSpoke);
    auto fin = srep::util::finally([pt](){ pt->Delete(); });
    auto smartPt = vtkSRepSkeletalPoint::SmartCreate(upSpoke, downSpoke);

    // shallow copies of the spokes, note EXPECT_EQ not EXPECT_SPOKE_EQ
    EXPECT_EQ(pt->GetUpSpoke(), upSpoke);
    EXPECT_EQ(pt->GetDownSpoke(), downSpoke);
    EXPECT_EQ(pt->GetCrestSpoke(), nullptr);
    EXPECT_FALSE(pt->IsCrest());
    EXPECT_EQ(smartPt->GetUpSpoke(), upSpoke);
    EXPECT_EQ(smartPt->GetUpSpoke(), upSpoke);
    EXPECT_EQ(smartPt->GetCrestSpoke(), nullptr);
    EXPECT_FALSE(smartPt->IsCrest());
  }
  {
    auto upSpoke = vtkSRepSpoke::SmartCreate(srep::Point3d(0,0,0), srep::Point3d(1,1,1));
    auto downSpoke = vtkSRepSpoke::SmartCreate(srep::Point3d(0,0,0), srep::Point3d(-1,-1,-1));
    auto crestSpoke = vtkSRepSpoke::SmartCreate(srep::Point3d(0,1,0), srep::Point3d(0,5,0));

    auto pt = vtkSRepSkeletalPoint::Create(upSpoke, downSpoke, crestSpoke);
    auto fin = srep::util::finally([pt](){ pt->Delete(); });
    auto smartPt = vtkSRepSkeletalPoint::SmartCreate(upSpoke, downSpoke, crestSpoke);

    // shallow copies of the spokes, note EXPECT_EQ not EXPECT_SPOKE_EQ
    EXPECT_EQ(pt->GetUpSpoke(), upSpoke);
    EXPECT_EQ(pt->GetDownSpoke(), downSpoke);
    EXPECT_EQ(pt->GetCrestSpoke(), crestSpoke);
    EXPECT_TRUE(pt->IsCrest());
    EXPECT_EQ(smartPt->GetUpSpoke(), upSpoke);
    EXPECT_EQ(smartPt->GetUpSpoke(), upSpoke);
    EXPECT_EQ(smartPt->GetCrestSpoke(), crestSpoke);
    EXPECT_TRUE(smartPt->IsCrest());
  }
}

TEST(SkeletalPointTest, SetGet) {
  using SO = vtkSRepSkeletalPoint::SpokeOrientation;
  auto upSpoke = vtkSRepSpoke::SmartCreate(srep::Point3d(0,0,0), srep::Point3d(1,1,1));
  auto downSpoke = vtkSRepSpoke::SmartCreate(srep::Point3d(0,0,0), srep::Point3d(-1,-1,-1));
  auto crestSpoke = vtkSRepSpoke::SmartCreate(srep::Point3d(0,1,0), srep::Point3d(0,5,0));
  auto pt = vtkSmartPointer<vtkSRepSkeletalPoint>::New();

  EXPECT_NE(upSpoke, pt->GetUpSpoke());
  pt->SetUpSpoke(upSpoke);
  EXPECT_EQ(upSpoke, pt->GetUpSpoke());
  EXPECT_EQ(upSpoke, pt->GetSpoke(SO::UpOrientation));

  EXPECT_NE(downSpoke, pt->GetDownSpoke());
  pt->SetDownSpoke(downSpoke);
  EXPECT_EQ(downSpoke, pt->GetDownSpoke());
  EXPECT_EQ(downSpoke, pt->GetSpoke(SO::DownOrientation));

  EXPECT_NE(crestSpoke, pt->GetCrestSpoke());
  EXPECT_FALSE(pt->IsCrest());
  pt->SetCrestSpoke(crestSpoke);
  EXPECT_EQ(crestSpoke, pt->GetCrestSpoke());
  EXPECT_EQ(crestSpoke, pt->GetSpoke(SO::CrestOrientation));
  EXPECT_TRUE(pt->IsCrest());

  // throw if nullptr set up or down
  EXPECT_THROW(pt->SetUpSpoke(nullptr), std::invalid_argument);
  EXPECT_EQ(upSpoke, pt->GetUpSpoke());
  EXPECT_EQ(upSpoke, pt->GetSpoke(SO::UpOrientation));
  EXPECT_THROW(pt->SetDownSpoke(nullptr), std::invalid_argument);
  EXPECT_EQ(downSpoke, pt->GetDownSpoke());
  EXPECT_EQ(downSpoke, pt->GetSpoke(SO::DownOrientation));
  EXPECT_NO_THROW(pt->SetCrestSpoke(nullptr));
  EXPECT_EQ(nullptr, pt->GetCrestSpoke());
  EXPECT_EQ(nullptr, pt->GetSpoke(SO::CrestOrientation));
  EXPECT_FALSE(pt->IsCrest());
}

TEST(SkeletalPointTest, SetGetBySpokeOrientation) {
  using SO = vtkSRepSkeletalPoint::SpokeOrientation;
  auto upSpoke = vtkSRepSpoke::SmartCreate(srep::Point3d(0,0,0), srep::Point3d(1,1,1));
  auto downSpoke = vtkSRepSpoke::SmartCreate(srep::Point3d(0,0,0), srep::Point3d(-1,-1,-1));
  auto crestSpoke = vtkSRepSpoke::SmartCreate(srep::Point3d(0,1,0), srep::Point3d(0,5,0));
  auto pt = vtkSmartPointer<vtkSRepSkeletalPoint>::New();

  EXPECT_NE(upSpoke, pt->GetSpoke(SO::UpOrientation));
  pt->SetSpoke(SO::UpOrientation, upSpoke);
  EXPECT_EQ(upSpoke, pt->GetUpSpoke());
  EXPECT_EQ(upSpoke, pt->GetSpoke(SO::UpOrientation));

  EXPECT_NE(downSpoke, pt->GetSpoke(SO::DownOrientation));
  pt->SetSpoke(SO::DownOrientation, downSpoke);
  EXPECT_EQ(downSpoke, pt->GetDownSpoke());
  EXPECT_EQ(downSpoke, pt->GetSpoke(SO::DownOrientation));

  EXPECT_NE(crestSpoke, pt->GetSpoke(SO::CrestOrientation));
  EXPECT_FALSE(pt->IsCrest());
  pt->SetSpoke(SO::CrestOrientation, crestSpoke);
  EXPECT_EQ(crestSpoke, pt->GetCrestSpoke());
  EXPECT_EQ(crestSpoke, pt->GetSpoke(SO::CrestOrientation));
  EXPECT_TRUE(pt->IsCrest());

  // throw if nullptr set up or down
  EXPECT_THROW(pt->SetSpoke(SO::UpOrientation, nullptr), std::invalid_argument);
  EXPECT_EQ(upSpoke, pt->GetUpSpoke());
  EXPECT_EQ(upSpoke, pt->GetSpoke(SO::UpOrientation));
  EXPECT_THROW(pt->SetSpoke(SO::DownOrientation, nullptr), std::invalid_argument);
  EXPECT_EQ(downSpoke, pt->GetDownSpoke());
  EXPECT_EQ(downSpoke, pt->GetSpoke(SO::DownOrientation));
  EXPECT_NO_THROW(pt->SetSpoke(SO::CrestOrientation, nullptr));
  EXPECT_EQ(nullptr, pt->GetCrestSpoke());
  EXPECT_EQ(nullptr, pt->GetSpoke(SO::CrestOrientation));
  EXPECT_FALSE(pt->IsCrest());
}

TEST(SkeletalPointTest, Observation) {
  // observation of default
  auto pt = vtkSmartPointer<vtkSRepSkeletalPoint>::New();
  TestObserver obs;
  const auto tag = pt->AddObserver(vtkCommand::ModifiedEvent, &obs, &TestObserver::callback);
  auto tagFin = srep::util::finally([tag, pt](){ pt->RemoveObserver(tag); });
  EXPECT_EQ(0, obs.numCalls());

  // changing spokes calls modified
  pt->GetUpSpoke()->SetDirectionAndMagnitude(srep::Vector3d(1,1,1));
  EXPECT_EQ(1, obs.numCalls());
  pt->GetDownSpoke()->SetDirectionAndMagnitude(srep::Vector3d(1,1,1));
  EXPECT_EQ(2, obs.numCalls());

  // gettings spokes and not changing them doesn't
  vtkSmartPointer<vtkSRepSpoke> originalUpSpoke = pt->GetUpSpoke();
  EXPECT_EQ(2, obs.numCalls());
  vtkSmartPointer<vtkSRepSpoke> originalDownSpoke = pt->GetDownSpoke();
  EXPECT_EQ(2, obs.numCalls());

  // setting new spokes does
  auto newUpSpoke = vtkSmartPointer<vtkSRepSpoke>::New();
  pt->SetUpSpoke(newUpSpoke);
  EXPECT_EQ(3, obs.numCalls());

  auto newDownSpoke = vtkSmartPointer<vtkSRepSpoke>::New();
  pt->SetDownSpoke(newDownSpoke);
  EXPECT_EQ(4, obs.numCalls());

  auto newCrestSpoke = vtkSmartPointer<vtkSRepSpoke>::New();
  pt->SetCrestSpoke(newCrestSpoke);
  EXPECT_EQ(5, obs.numCalls());

  // updating the old spokes doesn't
  originalUpSpoke->SetDirectionAndMagnitude(srep::Vector3d(1,1,2));
  EXPECT_EQ(5, obs.numCalls());
  originalDownSpoke->SetDirectionAndMagnitude(srep::Vector3d(1,1,2));
  EXPECT_EQ(5, obs.numCalls());

  // but updating the new ones does
  newUpSpoke->SetDirectionAndMagnitude(srep::Vector3d(2,1,2));
  EXPECT_EQ(6, obs.numCalls());
  newDownSpoke->SetDirectionAndMagnitude(srep::Vector3d(2,1,2));
  EXPECT_EQ(7, obs.numCalls());

  // all methods that shouldn't update the modified
  EXPECT_TRUE(pt->IsCrest());
  EXPECT_EQ(7, obs.numCalls());
  std::stringstream ss;
  pt->PrintSelf(ss, vtkIndent());

  auto clone = vtkSmartPointer<vtkSRepSkeletalPoint>::Take(pt->Clone());
  EXPECT_NE(pt, clone);
  auto smartClone = pt->SmartClone();
  EXPECT_NE(pt, smartClone);

  EXPECT_EQ(7, obs.numCalls());
}

TEST(SkeletalPointTest, ObservationBySpokeOrientation) {
  using SO = vtkSRepSkeletalPoint::SpokeOrientation;

  // observation of default
  auto pt = vtkSmartPointer<vtkSRepSkeletalPoint>::New();
  TestObserver obs;
  const auto tag = pt->AddObserver(vtkCommand::ModifiedEvent, &obs, &TestObserver::callback);
  auto tagFin = srep::util::finally([tag, pt](){ pt->RemoveObserver(tag); });
  EXPECT_EQ(0, obs.numCalls());

  // changing spokes calls modified
  pt->GetSpoke(SO::UpOrientation)->SetDirectionAndMagnitude(srep::Vector3d(1,1,1));
  EXPECT_EQ(1, obs.numCalls());
  pt->GetSpoke(SO::DownOrientation)->SetDirectionAndMagnitude(srep::Vector3d(1,1,1));
  EXPECT_EQ(2, obs.numCalls());

  // gettings spokes and not changing them doesn't
  vtkSmartPointer<vtkSRepSpoke> originalUpSpoke = pt->GetSpoke(SO::UpOrientation);
  EXPECT_EQ(2, obs.numCalls());
  vtkSmartPointer<vtkSRepSpoke> originalDownSpoke = pt->GetSpoke(SO::DownOrientation);
  EXPECT_EQ(2, obs.numCalls());

  // setting new spokes does
  auto newUpSpoke = vtkSmartPointer<vtkSRepSpoke>::New();
  pt->SetSpoke(SO::UpOrientation, newUpSpoke);
  EXPECT_EQ(3, obs.numCalls());

  auto newDownSpoke = vtkSmartPointer<vtkSRepSpoke>::New();
  pt->SetSpoke(SO::DownOrientation, newDownSpoke);
  EXPECT_EQ(4, obs.numCalls());

  auto newCrestSpoke = vtkSmartPointer<vtkSRepSpoke>::New();
  pt->SetSpoke(SO::CrestOrientation, newCrestSpoke);
  EXPECT_EQ(5, obs.numCalls());

  // updating the old spokes doesn't
  originalUpSpoke->SetDirectionAndMagnitude(srep::Vector3d(1,1,2));
  EXPECT_EQ(5, obs.numCalls());
  originalDownSpoke->SetDirectionAndMagnitude(srep::Vector3d(1,1,2));
  EXPECT_EQ(5, obs.numCalls());

  // but updating the new ones does
  newUpSpoke->SetDirectionAndMagnitude(srep::Vector3d(2,1,2));
  EXPECT_EQ(6, obs.numCalls());
  newDownSpoke->SetDirectionAndMagnitude(srep::Vector3d(2,1,2));
  EXPECT_EQ(7, obs.numCalls());
}

TEST(SkeletalPointTest, Clone) {
  using SO = vtkSRepSkeletalPoint::SpokeOrientation;

  auto upSpoke = vtkSRepSpoke::SmartCreate(srep::Point3d(0,0,0), srep::Point3d(1,1,1));
  auto downSpoke = vtkSRepSpoke::SmartCreate(srep::Point3d(0,0,0), srep::Point3d(-1,-1,-1));
  auto crestSpoke = vtkSRepSpoke::SmartCreate(srep::Point3d(0,1,0), srep::Point3d(0,5,0));

  auto pt = vtkSRepSkeletalPoint::SmartCreate(upSpoke, downSpoke, crestSpoke);
  auto clone = pt->Clone();
  auto fin = srep::util::finally([clone](){ clone->Delete(); });
  auto smartClone = pt->SmartClone();

  // make sure they are equal
  EXPECT_SKELETAL_POINT_EQ(pt, clone);
  EXPECT_SKELETAL_POINT_EQ(pt, smartClone);

  // but actually deep copies
  EXPECT_NE(pt, clone);
  EXPECT_NE(pt, smartClone);
  EXPECT_NE(pt->GetUpSpoke(), clone->GetUpSpoke());
  EXPECT_NE(pt->GetDownSpoke(), clone->GetDownSpoke());
  EXPECT_NE(pt->GetCrestSpoke(), clone->GetCrestSpoke());
  EXPECT_NE(pt->GetSpoke(SO::UpOrientation), clone->GetSpoke(SO::UpOrientation));
  EXPECT_NE(pt->GetSpoke(SO::DownOrientation), clone->GetSpoke(SO::DownOrientation));
  EXPECT_NE(pt->GetSpoke(SO::CrestOrientation), clone->GetSpoke(SO::CrestOrientation));
}

TEST(SkeletalPointTest, CloneWithObservers) {
  auto upSpoke = vtkSRepSpoke::SmartCreate(srep::Point3d(0,0,0), srep::Point3d(1,1,1));
  auto downSpoke = vtkSRepSpoke::SmartCreate(srep::Point3d(0,0,0), srep::Point3d(-1,-1,-1));
  auto crestSpoke = vtkSRepSpoke::SmartCreate(srep::Point3d(0,1,0), srep::Point3d(0,5,0));

  auto pt = vtkSRepSkeletalPoint::SmartCreate(upSpoke, downSpoke, crestSpoke);
  TestObserver obs;
  const auto tag = pt->AddObserver(vtkCommand::ModifiedEvent, &obs, &TestObserver::callback);
  auto tagFin = srep::util::finally([tag, pt](){ pt->RemoveObserver(tag); });
  EXPECT_EQ(0, obs.numCalls());

  auto clone = pt->Clone();
  auto fin = srep::util::finally([clone](){ clone->Delete(); });
  auto smartClone = pt->SmartClone();

  EXPECT_EQ(0, obs.numCalls());

  // clones don't affect original observer
  clone->SetCrestSpoke(nullptr);
  EXPECT_NE(nullptr, pt->GetCrestSpoke());
  EXPECT_NE(nullptr, smartClone->GetCrestSpoke());
  smartClone->GetUpSpoke()->SetRadius(1);

  EXPECT_EQ(0, obs.numCalls());

  // original doesn't affect clones
  pt->GetUpSpoke()->SetRadius(2);
  EXPECT_EQ(1, obs.numCalls());
  EXPECT_EQ(1, smartClone->GetUpSpoke()->GetRadius());
}

TEST(SkeletalPointTest, PrintSelf) {
  auto upSpoke = vtkSRepSpoke::SmartCreate(srep::Point3d(0,0,0), srep::Point3d(1,1,1));
  auto downSpoke = vtkSRepSpoke::SmartCreate(srep::Point3d(0,0,0), srep::Point3d(-1,-1,-1));
  auto crestSpoke = vtkSRepSpoke::SmartCreate(srep::Point3d(0,1,0), srep::Point3d(0,5,0));

  auto pt = vtkSRepSkeletalPoint::SmartCreate(upSpoke, downSpoke, crestSpoke);
  PrintSelfTest(*pt);
}
