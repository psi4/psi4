#ifndef _psi_src_lib_libmints__vector3_h
#define _psi_src_lib_libmints__vector3_h

#include <cmath>

namespace psi {
    
class Vector3
{
private:
    double v_[3];
    
public:
    Vector3() { v_[0] = v_[1] = v_[2] = 0.0; }
    Vector3(const double p[3]) {
        v_[0] = p[0]; v_[1] = p[1]; v_[2] = p[2];
    }
    Vector3(double d) {
        v_[0] = v_[1] = v_[2] = d;
    }
    Vector3(double x, double y, double z) {
        v_[0] = x; v_[1] = y; v_[2] = z;
    }
    Vector3(const Vector3& p) {
        v_[0] = p.v_[0]; v_[1] = p.v_[1]; v_[2] = p.v_[2];
    }
    
    void operator=(const double *x) {
        v_[0] = x[0];
        v_[1] = x[1];
        v_[2] = x[2];
    }
    void operator=(const Vector3& x) {
        v_[0] = x.v_[0];
        v_[1] = x.v_[1];
        v_[2] = x.v_[2];
    }
    void operator=(double d) {
        v_[0] = d;
        v_[1] = d;
        v_[2] = d;
    }
    void operator+=(const Vector3& x) {
        v_[0] += x.v_[0];
        v_[1] += x.v_[1];
        v_[2] += x.v_[2];
    }
    void operator-=(const Vector3& x) {
        v_[0] -= x.v_[0];
        v_[1] -= x.v_[1];
        v_[2] -= x.v_[2];
    }
    void operator*=(double m) {
        v_[0] *= m;
        v_[1] *= m;
        v_[2] *= m;
    }
    Vector3 operator+(const Vector3& x) {
        Vector3 result;
        result.v_[0] = v_[0] + x.v_[0];
        result.v_[1] = v_[1] + x.v_[1];
        result.v_[2] = v_[2] + x.v_[2];
        return result;
    }
    Vector3 operator-(const Vector3& x) {
        Vector3 result;
        result.v_[0] = v_[0] - x.v_[0];
        result.v_[1] = v_[1] - x.v_[1];
        result.v_[2] = v_[2] - x.v_[2];
        return result;
    }
    Vector3 operator-() { return Vector3(-v_[0], -v_[1], -v_[2]); }
    
    double& operator[](int i) { return v_[i]; }
    const double& operator[](int i) const { return v_[i]; }
    
    double dot(const Vector3& x) const {
        return v_[0]*x.v_[0] + v_[1]*x.v_[1] + v_[2]*x.v_[2];
    }
    
    double distance(const Vector3&) const;
    void normalize();
    double norm() const { return sqrt(this->dot(*this)); }
    
    Vector3 cross(const Vector3&) const;
};

Vector3 operator*(double, const Vector3&);

#include "vector3_i.cc"

}

#endif
