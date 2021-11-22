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
{}

SkeletalPoint::SkeletalPoint(const Spoke& upSpoke, const Spoke& downSpoke, const Spoke& crestSpoke)
    : UpSpoke(upSpoke)
    , DownSpoke(downSpoke)
    , CrestSpoke(crestSpoke)
    , HasCrestSpoke(true)
{}

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
    this->UpSpoke = spoke;
}
void SkeletalPoint::SetDownSpoke(const Spoke& spoke) {
    this->DownSpoke = spoke;
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
void SkeletalPoint::SetCrestSpoke(const Spoke& spoke) {
    this->SetCrestSpoke(&spoke);
}

bool SkeletalPoint::IsCrest() const {
    return this->HasCrestSpoke;
}

bool operator==(const SkeletalPoint& a, const SkeletalPoint& b) {
    bool ret = a.GetUpSpoke() == b.GetUpSpoke()
        && a.GetDownSpoke() == b.GetDownSpoke()
        && a.IsCrest() == b.IsCrest();

    if (ret && a.IsCrest()) {
        ret = ret && a.GetCrestSpoke() == b.GetCrestSpoke();
    }
    return ret;
}
bool operator!=(const SkeletalPoint& a, const SkeletalPoint& b) {
    return !(a == b);
}
bool operator< (const SkeletalPoint& a, const SkeletalPoint& b) {
    if (a.GetUpSpoke() != b.GetUpSpoke()) {
        return a.GetUpSpoke() < b.GetUpSpoke();
    }
    if (a.GetDownSpoke() != b.GetDownSpoke()) {
        return a.GetDownSpoke() < b.GetDownSpoke();
    }
    if (a.IsCrest() != b.IsCrest()) {
        return a.IsCrest() ? true : false;
    }
    //if both are crest
    if (a.IsCrest()) {
        return a.GetCrestSpoke() < b.GetCrestSpoke();
    }
    //both are not crest and up and down are same, therefore a and b are equal
    return false;
}
bool operator> (const SkeletalPoint& a, const SkeletalPoint& b) {
    return b < a;
}
bool operator<=(const SkeletalPoint& a, const SkeletalPoint& b) {
    return !(b < a);
}
bool operator>=(const SkeletalPoint& a, const SkeletalPoint& b) {
    return !(a < b);
}

std::ostream& operator<<(std::ostream& os, const SkeletalPoint& point) {
    os << "SkeletalPoint {" << std::endl
       << "  Up:    " << point.GetUpSpoke() << std::endl
       << "  Down:  " << point.GetDownSpoke() << std::endl;
    if (point.IsCrest()) {
        os << "  Crest: " << point.GetCrestSpoke() << std::endl;
    } else {
        os << "  Crest: None" << std::endl;
    }
    os << "}";
    return os;
}

}
