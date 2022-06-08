#include <fantom/math.hpp>

using namespace fantom;

class Utils {

    Utils(){
        return;
    }
    
    Vector<3> scaleVector(Vector<3> v, double scale) {
        return Vector<3>(v[0] * scale, v[1] * scale, v[2] * scale);
    }
};