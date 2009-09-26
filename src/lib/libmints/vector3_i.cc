#include <cmath>

inline Vector3 operator*(double d, const Vector3& x)
{
    Vector3 result;
    result[0] = d * x[0];
    result[1] = d * x[1];
    result[2] = d * x[2];
    return result;
}

inline double Vector3::distance(const Vector3& s) const
{
    double x = v_[0] - s.v_[0];
    double y = v_[1] - s.v_[1];
    double z = v_[2] - s.v_[2];
    return sqrt(x*x + y*y + z*z);
}

inline void Vector3::normalize()
{
    double temp=0.0;
    int i;
    
    for (i=0; i<3; ++i)
        temp += v_[i] * v_[i];
    temp = 1.0 / sqrt(temp);
    for (i=0; i<3; ++i)
        v_[i] *= temp;
}

inline Vector3 Vector3::cross(const Vector3& x) const
{
    Vector3 result(v_[1] * x.v_[2] - v_[2] * x.v_[1],
                   v_[2] * x.v_[0] - v_[0] * x.v_[2],
                   v_[0] * x.v_[1] - v_[1] * x.v_[0]);
    return result;
}

