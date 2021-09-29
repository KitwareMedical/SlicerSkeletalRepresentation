#ifndef __srep_SkeletalPoint_h
#define __srep_SkeletalPoint_h

#include <stdexcept>

#include <srep/Spoke.h>

namespace srep {

class MisalignedSpokeException : public std::invalid_argument {
public:
    using std::invalid_argument::invalid_argument;
};

class NotACrestException : public std::invalid_argument {
public:
    using std::invalid_argument::invalid_argument;
};

class SkeletalPoint {
public:
    SkeletalPoint();
    /// Constructs a non-crest skeletal point.
    /// @throws MisalignedSpokeException if upSpoke and downSpoke are not at the same skeletalPoint.
    SkeletalPoint(const Spoke& upSpoke, const Spoke& downSpoke);
    /// Constructs a crest skeletal point.
    /// @throws MisalignedSpokeException if upSpoke and downSpokeare are not at the same skeletalPoint.
    /// @note Crest spoke are also known as fold spokes.
    /// @note The crestSpoke does not need to be at the same skeletal point as the upSpoke and downSpoke.
    SkeletalPoint(const Spoke& upSpoke, const Spoke& downSpoke, const Spoke& crestSpoke);

    //copy and move defined and valid
    SkeletalPoint(const SkeletalPoint&) = default;
    SkeletalPoint& operator=(const SkeletalPoint&) = default;
    SkeletalPoint(SkeletalPoint&&) = default;
    SkeletalPoint& operator=(SkeletalPoint&&) = default;
    ~SkeletalPoint() = default;

    /// Gets the spoke in the up direction.
    const Spoke& GetUpSpoke() const;

    /// Getst the spoke in the down direction.
    const Spoke& GetDownSpoke() const;

    /// Returns whether or not this is a crest skeletal point.
    bool IsCrest() const;

    /// Gets the crest spoke.
    /// @throws NotACrestException if IsCrest returns false
    /// @note Crest spoke are also known as fold spokes.
    const Spoke& GetCrestSpoke() const;

    /// Sets the up spoke. If the incoming spoke is not aligned to the existing down spoke, this 
    /// object will not change.
    /// @throws MisalignedSpokeException if upSpoke and downSpoke are not at the same skeletalPoint
    void SetUpSpoke(const Spoke& spoke);
    /// Sets the down spoke. If the incoming spoke is not aligned to the existing up spoke, this 
    /// object will not change.
    /// @throws MisalignedSpokeException if upSpoke and downSpoke are not at the same skeletalPoint.
    void SetDownSpoke(const Spoke& spoke);
    /// Sets crest spoke. If nullptr is passed in, will remove crest spoke.
    /// @note Crest spoke are also known as fold spokes.
    /// @note This function does not take ownership of the pointer.
    void SetCrestSpoke(const Spoke* spoke);

    /// Gets the point on the skeleton for this skeletal point.
    Point3d GetPoint() const;
private:
    void Validate() const;

    Spoke UpSpoke;
    Spoke DownSpoke;
    Spoke CrestSpoke;
    bool HasCrestSpoke;
};

}

#endif