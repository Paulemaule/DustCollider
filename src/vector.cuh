#pragma once

///////////////////////// HOST CODE /////////////////////////

/** 
 * @brief Calculates the lenght of a quaternion.
 * 
 * @param q: The quaternion whose lenght is to be calculated.
 * @return The lenght of the quaternion.
 */
inline double cpu_quat_length(quat q)
{
    return sqrt(q.e0 * q.e0 + q.e1 * q.e1 + q.e2 * q.e2 + q.e3 * q.e3);
}

/**
 * @brief Normalizes a quaternion.
 * 
 * @param q&: Reference to the qauternion that is to be normalized.
 */
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

/**
 * @brief Calculates the lenght of a vector.
 * 
 * @param v: The vector whos lenght is to be calculated.
 * @return The lenght of the vector. 
 */
inline double cpu_vec3D_length(vec3D v)
{
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

/**
 * @brief Calculates the square of the lenght of a vector.
 * 
 * @param v: The vector whos lenght is to be calculated.
 * @return The squared lenght of the vector.
 */
inline double cpu_vec3D_length_sq(vec3D v)
{
    return v.x * v.x + v.y * v.y + v.z * v.z;
}

/**
 * @brief Checks if the sum of vector components is zero.
 * 
 * @param v: The vector that is to be checked.
 * @return Is the vector zero?
 */
inline bool cpu_vec3D_is_zero(vec3D v)
{
    return (v.x + v.y + v.z) == 0;
}

/**
 * @brief Calculates the lenght of the difference of two vectors.
 * 
 * @param u: One of the two vectors.
 * @param v: One of the two vectors.
 * @return The lenght of the difference.
 */
inline double cpu_vec3D_distance(vec3D u, vec3D v)
{
    double dx = u.x - v.x;
    double dy = u.y - v.y;
    double dz = u.z - v.z;

    return sqrt(dx * dx + dy * dy + dz * dz);
}

/**
 * @brief Calculates the dot product of two vectors.
 * 
 * @param u: One of the two vectors.
 * @param v: One of the two vectors.
 * @return The dot product of the two vectors.
 */
inline double cpu_vec3D_dot(vec3D u, vec3D v)
{
    return u.x * v.x + u.y * v.y + u.z * v.z;
}

/**
 * @brief Calculates the square of the lenght of the difference of two vectors.
 * 
 * @param u: One of the two vectors.
 * @param v: One of the two vectors.
 * @return The square of the lenght of the difference.
 */
inline double cpu_vec3D_distance_sq(vec3D u, vec3D v)
{
    double dx = u.x - v.x;
    double dy = u.y - v.y;
    double dz = u.z - v.z;

    return dx * dx + dy * dy + dz * dz;
}

/**
 * @brief Normalizes a vector.
 * 
 * @param &v: A reference to the vector that is to be normalized.
 * 
 * @warning When the vector is (0,0,0) the function will not do anything!
 */
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

/**
 * @brief Calculates the normal vector poiting from u to v.
 * 
 * @param u: The vector from which the normal points.
 * @param v: The vector to which the normal points.
 * @return The normal vector pointing form u to v.
 * 
 * @warning Undefined behavior when u = v! (Function will return (0,0,0) then.)
 */
inline vec3D cpu_vec3D_get_normal(vec3D u, vec3D v)
{
    vec3D res;
    res.x = u.x - v.x;
    res.y = u.y - v.y;
    res.z = u.z - v.z;

    cpu_vec3D_normalize(res);

    return res;
}

/**
 * @brief Normalizes a vector and returns the result.
 * 
 * @param v: The vector that is to be normalized.
 * @return The normalized vector.
 * 
 * @warning Undefined behavior when v = (0,0,0)!
 */
inline vec3D cpu_vec3D_get_normal(vec3D v)
{
    cpu_vec3D_normalize(v);

    return v;
}

/**
 * @brief Calculates the difference of two vectors,
 * 
 * @param u: The vector from which the difference points.
 * @param v: The vector to which the difference points.
 * @return The difference vector.
 */
//points from u to v
inline vec3D cpu_vec3D_diff(vec3D u, vec3D v)
{
    vec3D res;
    res.x = u.x - v.x;
    res.y = u.y - v.y;
    res.z = u.z - v.z;

    return res;
}

/**
 * @brief Calculates the cross product of two vectors.
 * 
 * @param u: The first vector (the index finger in the right hand rule).
 * @param v: The second vector (the middle finger in the right hand rule).
 * @return The cross product (u x v) of the two vectors (the thumb in the right hand rule).
 */
inline vec3D cpu_vec3D_cross(const vec3D u, const vec3D v)
{
    vec3D res;

    res.x = u.y * v.z - u.z * v.y;
    res.y = u.z * v.x - u.x * v.z;
    res.z = u.x * v.y - u.y * v.x;

    return res;
}

/**
 * @brief Sets the components of a vector to a specified value.
 * 
 * @param &v: Reference to the vector that is to be set.
 * @param val: The value that is to be used.
 */
inline void cpu_vec3D_set(vec3D& v, double val)
{
    v.x = val;
    v.y = val;
    v.z = val;
}

///////////////////////// DEVICE CODE /////////////////////////

/**
 * @brief Calculates the lenght of a quaternion
 * 
 * @param q: The quaternion whos lenght is to be calculated.
 * @return The lenght of the quaternion.
 */
__device__ double gpu_quat_length(quat q)
{
    return sqrt(q.e0 * q.e0 + q.e1 * q.e1 + q.e2 * q.e2 + q.e3 * q.e3);
}

/**
 * @brief Normalizes a quaternion.
 * 
 * @param q: The quaternion that is to be normalized.
 */
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

/**
 * @brief Calculates the lenght of a vector.
 * 
 * @param q: The vector whos lenght is to be calculated.
 * @return The lenght of the vector.
 */
__device__ double gpu_vec3D_length(vec3D v)
{
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

/**
 * @brief Calculates the squared lenght of a vector.
 * 
 * @param q: The vector whos lenght is to be calculated.
 * @return The squared lenght of the quaternion.
 */
__device__ double gpu_vec3D_length_sq(vec3D v)
{
    return v.x * v.x + v.y * v.y + v.z * v.z;
}

/**
 * @brief Checks if the sum of a vectors components is zero.
 * 
 * @param q: The vector that is to be checked.
 * @return Is the sum of the vectors components zero?
 */
__device__ bool gpu_vec3D_is_zero(vec3D v)
{
    return (v.x + v.y + v.z) == 0;
}

/**
 * @brief Calculates the difference between two vectors.
 * 
 * @param u: The vector from which the difference points.
 * @param v: The vector to which the difference points.
 * @return The difference vector (u - v).
 */
__device__ double gpu_vec3D_distance(vec3D u, vec3D v)
{
    double dx = u.x - v.x;
    double dy = u.y - v.y;
    double dz = u.z - v.z;

    return sqrt(dx * dx + dy * dy + dz * dz);
}

/**
 * @brief Calculated the dot product between two vectors.
 * 
 * @param u: The first vector.
 * @param v: The second vector.
 * @return The dot product (u * v).
 */
__device__ double gpu_vec3D_dot(vec3D u, vec3D v)
{
    return u.x * v.x + u.y * v.y + u.z * v.z;
}

/**
 * @brief Calculates the squared lenght of the difference of two vectors.
 * 
 * @param u: The first vector.
 * @param v: The secodn vector.
 * @return The squared lenght of the distance |u - v|².
 */
__device__ double gpu_vec3D_distance_sq(vec3D u, vec3D v)
{
    double dx = u.x - v.x;
    double dy = u.y - v.y;
    double dz = u.z - v.z;

    return dx * dx + dy * dy + dz * dz;
}

/**
 * @brief Normalizes a vector.
 * 
 * @param v: The vector that is to be normalized.
 * 
 * @warning This has undefined behavior at v = (0,0,0)! The function will leave the vector unchanged there.
 */
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

/**
 * @brief Calculates the normal vector pointing between two vectors.
 * 
 * @param u: The vector from which the normal will point.
 * @param v: The vector to which the normal will point.
 * @return The normal vector pointing from u to v.
 */
__device__ vec3D gpu_vec3D_get_normal(vec3D u, vec3D v)
{
    vec3D res;
    res.x = u.x - v.x;
    res.y = u.y - v.y;
    res.z = u.z - v.z;

    gpu_vec3D_normalize(res);

    return res;
}

/**
 * @brief Noramlizes a vector and returns it.
 * 
 * @param v: The vector that is to be normalized.
 * @return The normalized vector.
 */
__device__ vec3D gpu_vec3D_get_normal(vec3D v)
{
    gpu_vec3D_normalize(v);

    return v;
}

/**
 * @brief Calculates the difference between two vectors.
 * 
 * @param u: The vector from which the difference will point.
 * @param v: The vector to which the difference will point.
 * @return The difference vector (u - v).
 */
__device__ vec3D gpu_vec3D_diff(vec3D u, vec3D v)
{
    vec3D res;
    res.x = u.x - v.x;
    res.y = u.y - v.y;
    res.z = u.z - v.z;

    return res;
}

/**
 * @brief Calculates the cross product of two vectors.
 * 
 * @param u: The first vector (the index finger in the right hand rule).
 * @param v: The second vector (the middle finger in the right hand rule).
 * @return The cross product (u x v) of the two vectors (the thumb in the right hand rule).
 */
__device__ vec3D gpu_vec3D_cross(const vec3D u, const vec3D v)
{
    vec3D res;

    res.x = u.y * v.z - u.z * v.y;
    res.y = u.z * v.x - u.x * v.z;
    res.z = u.x * v.y - u.y * v.x;

    return res;
}

/**
 * @brief Sets the components of a vector to a specified value.
 * 
 * @param &v: A reference to the vector whos components are to be set.
 * @param val: The value the components are to be set to.
 */
__device__ void gpu_vec3D_set(vec3D& v, double val)
{
    v.x = val;
    v.y = val;
    v.z = val;
}