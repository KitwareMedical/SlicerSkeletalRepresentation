#ifndef __srep_Util_h
#define __srep_Util_h

namespace srep {
namespace util {

/// FinalAction class based off of guidelines support library.
/// http://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines#e19-use-a-final_action-object-to-express-cleanup-if-no-suitable-resource-handle-is-available
template <class T>
class FinalAction {
public:
    FinalAction(T action)
        : m_action(action)
    {}

    ~FinalAction() {
        m_action();
    }
private:
    T m_action;
};

/// finally function based off of guidelines support library.
///
/// Allows an arbitary action to be called when returned FinalAction goes out of scope.
/// http://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines#e19-use-a-final_action-object-to-express-cleanup-if-no-suitable-resource-handle-is-available
template <class T>
FinalAction<T> finally(T action) {
    return FinalAction<T>(action);
}

}
}

#endif