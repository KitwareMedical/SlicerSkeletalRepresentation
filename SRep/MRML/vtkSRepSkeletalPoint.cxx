#include <vtkSRepSkeletalPoint.h>

#include <vtkCommand.h>
#include <vtkObjectFactory.h>

namespace {
  template <class T>
  T* vtkSmartPointerRelease(vtkSmartPointer<T>& t) {
    // Register so object exists when smart pointer doesn't (object will be a raw owning pointer)
    auto ptr = t.Get();
    ptr->Register(nullptr);
    t = nullptr;
    return ptr;
  }
  template <class T>
  T* vtkSmartPointerRelease(vtkSmartPointer<T>&& t) {
    // Register so object exists when smart pointer doesn't (object will be a raw owning pointer)
    auto ptr = t.Get();
    ptr->Register(nullptr);
    t = nullptr;
    return ptr;
  }
}

//----------------------------------------------------------------------
vtkStandardNewMacro(vtkSRepSkeletalPoint);

//----------------------------------------------------------------------
vtkSRepSkeletalPoint::vtkSRepSkeletalPoint()
  : UpSpoke()
  , DownSpoke()
  , CrestSpoke()
{
  this->SetUpSpoke(vtkSmartPointer<vtkSRepSpoke>::New());
  this->SetDownSpoke(vtkSmartPointer<vtkSRepSpoke>::New());
}

//----------------------------------------------------------------------
vtkSRepSkeletalPoint::~vtkSRepSkeletalPoint() {
  if (this->UpSpoke) {
    this->UpSpoke->RemoveObserver(this->UpObservationTag);
    this->UpSpoke = nullptr;
  }
  if (this->DownSpoke) {
    this->DownSpoke->RemoveObserver(this->DownObservationTag);
    this->DownSpoke = nullptr;
  }
  if (this->CrestSpoke) {
    this->CrestSpoke->RemoveObserver(this->CrestObservationTag);
    this->CrestSpoke = nullptr;
  }
};

//----------------------------------------------------------------------
void vtkSRepSkeletalPoint::PrintSelf(std::ostream& os, vtkIndent indent) {
  this->Superclass::PrintSelf(os, indent);
  os << indent << *this;
}

//----------------------------------------------------------------------
std::ostream& operator<<(std::ostream& os, const vtkSRepSkeletalPoint& point) {
  os << "SkeletalPoint {" << std::endl
    << "  Up:    " << *point.GetUpSpoke() << std::endl
    << "  Down:  " << *point.GetDownSpoke() << std::endl;
  if (point.IsCrest()) {
    os << "  Crest: " << *point.GetCrestSpoke() << std::endl;
  } else {
    os << "  Crest: None" << std::endl;
  }
  os << "}";
  return os;
}

//----------------------------------------------------------------------
bool vtkSRepSkeletalPoint::IsCrest() const {
  return nullptr != this->GetCrestSpoke();
}

//----------------------------------------------------------------------
#define SPOKE_OPERATIONS(spokeType, isOptional) \
  const vtkSRepSpoke* vtkSRepSkeletalPoint::Get##spokeType##Spoke() const { \
    return this->spokeType##Spoke; \
  } \
  vtkSRepSpoke* vtkSRepSkeletalPoint::Get##spokeType##Spoke() { \
    return this->spokeType##Spoke; \
  } \
  void vtkSRepSkeletalPoint::Set##spokeType##Spoke(vtkSRepSpoke* spoke) { \
    if (!(isOptional) && !spoke) { \
      throw std::invalid_argument(#spokeType " spoke must not be nullptr"); \
    } \
    if (spoke != this->spokeType##Spoke.Get()) { \
      if (this->spokeType##Spoke) { \
        this->spokeType##Spoke->RemoveObserver(this->spokeType##ObservationTag); \
      } \
      this->spokeType##Spoke = spoke; \
      if (this->spokeType##Spoke) { \
        this->spokeType##ObservationTag = \
          this->spokeType##Spoke->AddObserver(vtkCommand::ModifiedEvent, this, &vtkSRepSkeletalPoint::onSpokeModified); \
      } \
      this->Modified(); \
    } \
  }

SPOKE_OPERATIONS(Up, false)
SPOKE_OPERATIONS(Down, false)
SPOKE_OPERATIONS(Crest, true)

//----------------------------------------------------------------------
const vtkSRepSpoke* vtkSRepSkeletalPoint::GetSpoke(SpokeOrientation spokeType) const {
  if (spokeType == UpOrientation) {
    return this->GetUpSpoke();
  } else if (spokeType == DownOrientation) {
    return this->GetDownSpoke();
  } else if (spokeType == CrestOrientation) {
    return this->GetCrestSpoke();
  } else {
    throw std::invalid_argument("Unknown spoke type: " + std::to_string(static_cast<int>(spokeType)));
  }
}

//----------------------------------------------------------------------
vtkSRepSpoke* vtkSRepSkeletalPoint::GetSpoke(SpokeOrientation spokeType) {
  return const_cast<vtkSRepSpoke*>(const_cast<const vtkSRepSkeletalPoint*>(this)->GetSpoke(spokeType));
}

//----------------------------------------------------------------------
void vtkSRepSkeletalPoint::SetSpoke(SpokeOrientation spokeType, vtkSRepSpoke* spoke) {
  if (spokeType == UpOrientation) {
    this->SetUpSpoke(spoke);
  } else if (spokeType == DownOrientation) {
    this->SetDownSpoke(spoke);
  } else if (spokeType == CrestOrientation) {
    this->SetCrestSpoke(spoke);
  } else {
    throw std::invalid_argument("Unknown spoke type: " + std::to_string(static_cast<int>(spokeType)));
  }
}

//----------------------------------------------------------------------
void vtkSRepSkeletalPoint::onSpokeModified(vtkObject */*caller*/, unsigned long /*event*/, void* /*callData*/) {
  this->Modified();
}

//----------------------------------------------------------------------
vtkSRepSkeletalPoint* vtkSRepSkeletalPoint::Create(vtkSRepSpoke* upSpoke, vtkSRepSpoke* downSpoke, vtkSRepSpoke* crestSpoke) {
  return vtkSmartPointerRelease(SmartCreate(upSpoke, downSpoke, crestSpoke));
}

//----------------------------------------------------------------------
vtkSmartPointer<vtkSRepSkeletalPoint> vtkSRepSkeletalPoint::SmartCreate(vtkSRepSpoke* upSpoke, vtkSRepSpoke* downSpoke, vtkSRepSpoke* crestSpoke) {
  // start off with smart point so if there are any exceptions we don't get a memory leak
  auto smartSkeletalPoint = vtkSmartPointer<vtkSRepSkeletalPoint>::New();
  smartSkeletalPoint->SetUpSpoke(upSpoke);
  smartSkeletalPoint->SetDownSpoke(downSpoke);
  smartSkeletalPoint->SetCrestSpoke(crestSpoke);

  return smartSkeletalPoint;
}

//----------------------------------------------------------------------
vtkSRepSkeletalPoint* vtkSRepSkeletalPoint::Clone() const {
  return vtkSmartPointerRelease(this->SmartClone());
}

//----------------------------------------------------------------------
vtkSmartPointer<vtkSRepSkeletalPoint> vtkSRepSkeletalPoint::SmartClone() const {
  auto clone = vtkSmartPointer<vtkSRepSkeletalPoint>::New();
  clone->SetUpSpoke(vtkSmartPointer<vtkSRepSpoke>::Take(this->GetUpSpoke()->Clone()));
  clone->SetDownSpoke(vtkSmartPointer<vtkSRepSpoke>::Take(this->GetDownSpoke()->Clone()));
  if (this->IsCrest()) {
    clone->SetCrestSpoke(vtkSmartPointer<vtkSRepSpoke>::Take(this->GetCrestSpoke()->Clone()));
  }
  return clone;
}
