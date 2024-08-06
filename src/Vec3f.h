//
// Created by ajy on 24-7-14.
//

#ifndef MESH_SIMPLIFICATION_VEC3F_H
#define MESH_SIMPLIFICATION_VEC3F_H

#include <cmath>

class Vec3f {
public:
    double x, y, z;

    inline Vec3f() : x(0), y(0), z(0) {}

    inline Vec3f(const double a, const double b, const double c)
            : x(a), y(b), z(c) {}

    inline Vec3f operator + (const Vec3f &a ) const {
        return { x + a.x, y + a.y, z + a.z };
    }

    inline Vec3f& operator += (const Vec3f &a) {
        x += a.x;  y += a.y;  z += a.z;
        return *this;
    }

    inline Vec3f operator * (const double a ) const {
        return { x * a, y * a, z * a };
    }

    inline Vec3f operator * (const Vec3f &a ) const {
        return { x * a.x, y * a.y, z * a.z };
    }

    inline Vec3f& operator = (const Vec3f &a ) {
        if (this != &a) {
            x = a.x;
            y = a.y;
            z = a.z;
        }
        return *this;
    }

    inline Vec3f operator / (const Vec3f &a ) const {
        return { x / a.x, y / a.y, z / a.z };
    }

    inline Vec3f operator - (const Vec3f &a ) const {
        return { x - a.x, y - a.y, z - a.z };
    }

    inline Vec3f operator / (const double a ) const {
        return { x / a, y / a, z / a };
    }

    [[nodiscard]] inline double dot( const Vec3f& a ) const {
        return a.x * x + a.y * y + a.z * z;
    }

    inline Vec3f cross(const Vec3f& a , const Vec3f& b )
    {
        x = a.y * b.z - a.z * b.y;
        y = a.z * b.x - a.x * b.z;
        z = a.x * b.y - a.y * b.x;
        return *this;
    }

    inline Vec3f cross(const Vec3f& other) const {
        return {
                y * other.z - z * other.y,
                z * other.x - x * other.z,
                x * other.y - y * other.x
        };
    }

    inline double angle( const Vec3f& v )
    {
        Vec3f a = v , b = *this;
        double dot = v.x*x + v.y*y + v.z*z;
        double len = a.length() * b.length();
        if (len == 0) len = 0.00001f;
        double input = std::max(-1.0, std::min(1.0, dot / len));
        return acos(input);
    }//计算本向量与v的夹角，弧度制

    inline Vec3f rot_x(double a )
    {
        double yy = cos ( a ) * y + sin ( a ) * z;
        double zz = cos ( a ) * z - sin ( a ) * y;
        y = yy; z = zz;
        return *this;
    }//绕X轴旋转a，弧度制

    inline Vec3f rot_y(double a )
    {
        double xx = cos ( -a ) * x + sin ( -a ) * z;
        double zz = cos ( -a ) * z - sin ( -a ) * x;
        x = xx; z = zz;
        return *this;
    }//绕Y轴旋转a，弧度制

    inline Vec3f rot_z(double a )
    {
        double yy = cos ( a ) * y + sin ( a ) * x;
        double xx = cos ( a ) * x - sin ( a ) * y;
        y = yy; x = xx;
        return *this;
    }//绕Z轴旋转a，弧度制

    [[nodiscard]] inline double length() const
    {
        return sqrt(x*x + y*y + z*z);
    }

    inline void normalize() {
        double length = std::sqrt(x * x + y * y + z * z);
        if (length != 0) {
            x /= length;
            y /= length;
            z /= length;
        }
    }
};

#endif //MESH_SIMPLIFICATION_VEC3F_H
