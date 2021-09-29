#include <srep/SkeletalPoint.h>

namespace srep {

SkeletalPoint::SkeletalPoint()
    : UpSpoke()
    , DownSpoke()
    , CrestSpoke()
    , HasCrestSpoke(false)
{}

SkeletalPoint::SkeletalPoint(const Spoke& upSpoke, const Spoke& downSpoke)
    : UpSpoke(upSpoke)
    , DownSpoke(downSpoke)
    , CrestSpoke()
    , HasCrestSpoke(false)
{
    this->Validate();
}

SkeletalPoint::SkeletalPoint(const Spoke& upSpoke, const Spoke& downSpoke, const Spoke& crestSpoke)
    : UpSpoke(upSpoke)
    , DownSpoke(downSpoke)
    , CrestSpoke(crestSpoke)
    , HasCrestSpoke(true)
{
    this->Validate();
}

void SkeletalPoint::Validate() const {
    if (this->GetUpSpoke().GetSkeletalPoint() != this->GetDownSpoke().GetSkeletalPoint())
    {
        throw MisalignedSpokeException("Up and down spokes on must have the same skeletal point");
    }
}

const Spoke& SkeletalPoint::GetUpSpoke() const {
    return this->UpSpoke;
}
const Spoke& SkeletalPoint::GetDownSpoke() const {
    return this->DownSpoke;
}
const Spoke& SkeletalPoint::GetCrestSpoke() const {
    if (!this->IsCrest()) {
        throw NotACrestException("This is not a crest point");
    }
    return this->CrestSpoke;
}


void SkeletalPoint::SetUpSpoke(const Spoke& spoke) {
    Spoke oldSpoke = this->UpSpoke;
    try {
        this->UpSpoke = spoke;
        this->Validate();
    } catch (...) {
        this->UpSpoke = oldSpoke;
        throw;
    }
}
void SkeletalPoint::SetDownSpoke(const Spoke& spoke) {
    Spoke oldSpoke = this->DownSpoke;
    try {
        this->DownSpoke = spoke;
        this->Validate();
    } catch (...) {
        this->DownSpoke = oldSpoke;
        throw;
    }
}
void SkeletalPoint::SetCrestSpoke(const Spoke* spoke) {
    if (spoke) {
        this->CrestSpoke = *spoke;
        HasCrestSpoke = true;
    } else {
        this->CrestSpoke = Spoke();
        HasCrestSpoke = false;
    }
}

bool SkeletalPoint::IsCrest() const {
    return this->HasCrestSpoke;
}

Point3d SkeletalPoint::GetPoint() const {
    // leveraging fact that all are at the same point
    return this->GetUpSpoke().GetSkeletalPoint();
}

}
