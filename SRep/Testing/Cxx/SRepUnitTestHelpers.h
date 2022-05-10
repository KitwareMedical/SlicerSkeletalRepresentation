#ifndef srepModuleUnitTestHelpers_h
#define srepModuleUnitTestHelpers_h

#include <vtkObject.h>
#include <vector>

#define EXPECT_SPOKE_EQ(S1, S2) \
  EXPECT_EQ((S1)->GetSkeletalPoint(), (S2)->GetSkeletalPoint()); \
  EXPECT_EQ((S1)->GetBoundaryPoint(), (S2)->GetBoundaryPoint()); \
  EXPECT_EQ((S1)->GetDirection(), (S2)->GetDirection()); \
  EXPECT_EQ((S1)->GetRadius(), (S2)->GetRadius())

#define EXPECT_SKELETAL_POINT_EQ(S1, S2) \
  do { \
    EXPECT_SPOKE_EQ((S1)->GetUpSpoke(), (S2)->GetUpSpoke()); \
    EXPECT_SPOKE_EQ((S1)->GetDownSpoke(), (S2)->GetDownSpoke()); \
    EXPECT_EQ((S1)->IsCrest(), (S2)->IsCrest()); \
    if ((S1)->IsCrest() && (S2)->IsCrest()) { \
      EXPECT_SPOKE_EQ((S1)->GetCrestSpoke(), (S2)->GetCrestSpoke()); \
    } \
  } while(false)

namespace srepUnitTestHelpers {

struct TestObserver {
  void callback(vtkObject *caller, unsigned long event, void* callData) {
    callers.push_back(caller);
    events.push_back(event);
    callDatas.push_back(callData);
  }

  std::vector<vtkObject*> callers;
  std::vector<unsigned long> events;
  std::vector<void*> callDatas;

  size_t numCalls() const {
    return callers.size();
  }
};

template <class T>
void PrintSelfTest(T& t) {
  // pretty much just make sure it prints something more than its superclass
  std::stringstream ss;
  t.PrintSelf(ss, vtkIndent());

  std::stringstream superSS;
  t.T::Superclass::PrintSelf(superSS, vtkIndent());

  EXPECT_GT(ss.str().length(), superSS.str().length());
}

}

#endif
