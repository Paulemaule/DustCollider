#include "typedefs.cuh"

#ifndef CVECTOR
#define CVECTOR

//  -> CPU implementation
inline double cpu_quat_length(quat q)
{
    return sqrt(q.e0 * q.e0 + q.e1 * q.e1 + q.e2 * q.e2 + q.e3 * q.e3);
}

inline void cpu_quat_normalize(quat& q)
{
    double length = cpu_quat_length(q);

    if (length > 0)
    {
        q.e0 /= length;
        q.e1 /= length;
        q.e2 /= length;
        q.e3 /= length;
    }
}

inline double cpu_vec3D_length(vec3D v)
{
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

inline double cpu_vec3D_length_sq(vec3D v)
{
    return v.x * v.x + v.y * v.y + v.z * v.z;
}

inline bool cpu_vec3D_is_zero(vec3D v)
{
    return (v.x + v.y + v.z) == 0;
}

inline double cpu_vec3D_distance(vec3D u, vec3D v)
{
    double dx = u.x - v.x;
    double dy = u.y - v.y;
    double dz = u.z - v.z;

    return sqrt(dx * dx + dy * dy + dz * dz);
}

inline double cpu_vec3D_dot(vec3D u, vec3D v)
{
    return u.x * v.x + u.y * v.y + u.z * v.z;
}

inline double cpu_vec3D_distance_sq(vec3D u, vec3D v)
{
    double dx = u.x - v.x;
    double dy = u.y - v.y;
    double dz = u.z - v.z;

    return dx * dx + dy * dy + dz * dz;
}

inline void cpu_vec3D_normalize(vec3D& v)
{
    double length = cpu_vec3D_length(v);

    if (length > 0)
    {
        v.x /= length;
        v.y /= length;
        v.z /= length;
    }
}

//normal that points from u to v
inline vec3D cpu_vec3D_get_normal(vec3D u, vec3D v)
{
    vec3D res;
    res.x = u.x - v.x;
    res.y = u.y - v.y;
    res.z = u.z - v.z;

    cpu_vec3D_normalize(res);

    return res;
}

//get normal of v
inline vec3D cpu_vec3D_get_normal(vec3D v)
{
    cpu_vec3D_normalize(v);

    return v;
}

//points from u to v
inline vec3D cpu_vec3D_diff(vec3D u, vec3D v)
{
    vec3D res;
    res.x = u.x - v.x;
    res.y = u.y - v.y;
    res.z = u.z - v.z;

    return res;
}

inline vec3D cpu_vec3D_cross(const vec3D u, const vec3D v)
{
    vec3D res;

    res.x = u.y * v.z - u.z * v.y;
    res.y = u.z * v.x - u.x * v.z;
    res.z = u.x * v.y - u.y * v.x;

    return res;
}

inline void cpu_vec3D_set(vec3D& v, double val)
{
    v.x = val;
    v.y = val;
    v.z = val;
}


//  -> GPU implementation

__device__ double gpu_quat_length(quat q)
{
    return sqrt(q.e0 * q.e0 + q.e1 * q.e1 + q.e2 * q.e2 + q.e3 * q.e3);
}

__device__ void gpu_quat_normalize(quat& q)
{
    double length = gpu_quat_length(q);

    if (length > 0)
    {
        q.e0 /= length;
        q.e1 /= length;
        q.e2 /= length;
        q.e3 /= length;
    }
}

__device__ double gpu_vec3D_length(vec3D v)
{
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

__device__ double gpu_vec3D_length_sq(vec3D v)
{
    return v.x * v.x + v.y * v.y + v.z * v.z;
}

__device__ bool gpu_vec3D_is_zero(vec3D v)
{
    return (v.x + v.y + v.z) == 0;
}

__device__ double gpu_vec3D_distance(vec3D u, vec3D v)
{
    double dx = u.x - v.x;
    double dy = u.y - v.y;
    double dz = u.z - v.z;

    return sqrt(dx * dx + dy * dy + dz * dz);
}

__device__ double gpu_vec3D_dot(vec3D u, vec3D v)
{
    return u.x * v.x + u.y * v.y + u.z * v.z;
}

__device__ double gpu_vec3D_distance_sq(vec3D u, vec3D v)
{
    double dx = u.x - v.x;
    double dy = u.y - v.y;
    double dz = u.z - v.z;

    return dx * dx + dy * dy + dz * dz;
}

__device__ void gpu_vec3D_normalize(vec3D& v)
{
    double length = gpu_vec3D_length(v);

    if (length > 0)
    {
        v.x /= length;
        v.y /= length;
        v.z /= length;
    }
}

//normal that points from u to v
__device__ vec3D gpu_vec3D_get_normal(vec3D u, vec3D v)
{
    vec3D res;
    res.x = u.x - v.x;
    res.y = u.y - v.y;
    res.z = u.z - v.z;

    gpu_vec3D_normalize(res);

    return res;
}

//get normal of v
__device__ vec3D gpu_vec3D_get_normal(vec3D v)
{
    gpu_vec3D_normalize(v);

    return v;
}

//points from u to v
__device__ vec3D gpu_vec3D_diff(vec3D u, vec3D v)
{
    vec3D res;
    res.x = u.x - v.x;
    res.y = u.y - v.y;
    res.z = u.z - v.z;

    return res;
}

__device__ vec3D gpu_vec3D_cross(const vec3D u, const vec3D v)
{
    vec3D res;

    res.x = u.y * v.z - u.z * v.y;
    res.y = u.z * v.x - u.x * v.z;
    res.z = u.x * v.y - u.y * v.x;

    return res;
}

__device__ void gpu_vec3D_set(vec3D& v, double val)
{
    v.x = val;
    v.y = val;
    v.z = val;
}

#endif