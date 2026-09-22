#include <Eigen/Geometry>
#include <numbers>
#include "GeometryHelper.h"

Vector3F GeometryHelper::rotateVectorAroundAxisAndCenter(const Vector3F &vec, const Vector3F &axis, double angle_rad,
                                                         const Vector3F &center)
{
    // Normalize the rotation axis
    Vector3F unit_axis = axis.normalized();

    // Translate the vector such that the center of rotation becomes the origin
    Vector3F translated_vec = vec - center;

    // Create a quaternion representing the rotation
    Eigen::AngleAxisd rotation(angle_rad, unit_axis);

    // Rotate the vector
    Vector3F rotatedVec = rotation * translated_vec;

    // Translate back
    Vector3F result = rotatedVec + center;

    return result;
}

Vector2F GeometryHelper::rotateVector2D(const Vector2F &vec, double angle_rad)
{
    double x = vec.x();
    double y = vec.y();
    double cost = cos(angle_rad);
    double sint = sin(angle_rad);
    return {x * cost - y * sint, x * sint + y * cost};
}

bool GeometryHelper::intersectLinesPointDirection2D(const Vector2F &p0, const Vector2F &d0, const Vector2F &p1,
                                                    const Vector2F &d1, Vector2F &t)
{
    // Line 1: P = p0 + t0 * d0
    // Line 2: Q = p1 + t1 * d1
    // For intersection: p0 + t0 * d0 = p1 + t1 * d1
    // => t0 * d0 - t1 * d1 = p1 - p0

    Matrix<double, 2, 2> A;
    A << d0, -d1;
    Vector2F b = p1 - p0;

    // Check if lines are parallel (determinant is zero)
    double det = A.determinant();
    if (std::abs(det) < 1e-8)
    {
        // Lines are parallel or coincident, no unique intersection
        return false;
    }

    // Solve for t0 and t1
    t = A.inverse() * b;
    return true;
}

bool GeometryHelper::intersectLineAndPlane3D(const Vector3F &point, const Vector3F &line, const Vector4F &plane, double &t)
{
    Vector3F n = plane.head<3>();
    if (fabs(line.dot(n)) < 1e-8)
    {
        // Line is parallel to plane
        return false;
    }

    t = -(point.dot(n) + plane(3)) / line.dot(n);
    return true;
}

bool GeometryHelper::pointInTriangle3D(const Vector3F &point, const Vector3F &v0, const Vector3F &v1,
                                       const Vector3F &v2)
{
    // Compute vectors
    Vector3F u = v1 - v0;
    Vector3F v = v2 - v0;
    Vector3F w = point - v0;

    // Compute dot products
    double uu = u.dot(u);
    double uv = u.dot(v);
    double vv = v.dot(v);
    double wu = w.dot(u);
    double wv = w.dot(v);

    // Compute barycentric coordinates
    double denominator = uv * uv - uu * vv;
    double s = (uv * wv - vv * wu) / denominator;
    double t = (uv * wu - uu * wv) / denominator;

    // Check if point is in triangle
    return (s >= 0.0) && (t >= 0.0) && (s + t <= 1.0);
}

double GeometryHelper::windingNumber2D(const Vector2F &point, const Vector2F &edge0, const Vector2F &edge1)
{
    double x0 = edge0(0);
    double y0 = edge0(1);
    double x1 = edge1(0);
    double y1 = edge1(1);

    double winding_number = atan2(y1 - point.y(), x1 - point.x()) - atan2(y0 - point.y(), x0 - point.x());
    if (winding_number > std::numbers::pi)
        winding_number -= 2 * std::numbers::pi;
    if (winding_number < -std::numbers::pi)
        winding_number += 2 * std::numbers::pi;
    return winding_number;
}

double GeometryHelper::windingNumber3D(const Vector3F &point, const Vector3F &face0, const Vector3F &face1,
                                       const Vector3F &face2)
{
    Vector3F a = face0 - point;
    Vector3F b = face1 - point;
    Vector3F c = face2 - point;

    double al = a.norm();
    double bl = b.norm();
    double cl = c.norm();

    Vector3F t = b.cross(c);
    double numerator = a.dot(t);
    double denominator = al * bl * cl + a.dot(b) * cl + b.dot(c) * al + c.dot(a) * bl;

    double omega = 2.0 * atan2(numerator, denominator);
    if (omega > 2 * std::numbers::pi)
        omega -= 4 * std::numbers::pi;
    if (omega < -2 * std::numbers::pi)
        omega += 4 * std::numbers::pi;
    return omega;
}

double GeometryHelper::triangleAreaWithOrigin2D(const Vector2F &p0, const Vector2F &p1)
{
    double x0 = p0.x();
    double y0 = p0.y();
    double x1 = p1.x();
    double y1 = p1.y();

    return 0.5 * x0 * y1 - 0.5 * x1 * y0;
}

bool GeometryHelper::triangleOrientation2D(const Vector2F &p0, const Vector2F &p1, const Vector2F &p2)
{
    return (p1(0) - p0(0)) * (p2(1) - p0(1)) - (p1(1) - p0(1)) * (p2(0) - p0(0)) > 0;
}

double GeometryHelper::tetVolumeWithOrigin3D(const Vector3F &p0, const Vector3F &p1, const Vector3F &p2)
{
    double x0 = p0.x();
    double y0 = p0.y();
    double z0 = p0.z();
    double x1 = p1.x();
    double y1 = p1.y();
    double z1 = p1.z();
    double x2 = p2.x();
    double y2 = p2.y();
    double z2 = p2.z();

    return (y1 * z2 - y2 * z1) * x0 / 6.0 + (-y0 * z2 + y2 * z0) * x1 / 6.0 + x2 * (y0 * z1 - y1 * z0) / 6.0;
}
