#pragma once

/**
 * @brief Calculates the lenght of a quaternion
 * 
 * @param q: The quaternion whos lenght is to be calculated.
 * @return The lenght of the quaternion.
 */
__device__ double quat_lenght(const double4 q)
{
    return sqrt(q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z);
}

/**
 * @brief Normalizes a quaternion.
 * 
 * @param q: The quaternion that is to be normalized.
 */
__device__ void quat_normalize(double4 q)
{
    double length = quat_lenght(q);

    if (length > 0)
    {
        q.w /= length;
        q.x /= length;
        q.y /= length;
        q.z /= length;
    }
}

/**
 * @brief Calculates the lenght of a vector.
 * 
 * @param q: The vector whos lenght is to be calculated.
 * @return The lenght of the vector.
 */
__device__ double vec_lenght(const double3 v)
{
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

/**
 * @brief Calculates the squared lenght of a vector.
 * 
 * @param q: The vector whos lenght is to be calculated.
 * @return The squared lenght of the quaternion.
 */
__device__ double vec_lenght_sq(const double3 v)
{
    return v.x * v.x + v.y * v.y + v.z * v.z;
}

/**
 * @brief Checks if a vector is the zero vector.
 * 
 * @param q: The vector that is to be checked.
 * @return Is the sum of the vectors components zero?
 * 
 * @warning This function will also return true if the sum of the vector components is zero.
 */
__device__ bool vec_is_zero(const double3 v)
{
    return (v.x + v.y + v.z) == 0;
}

/**
 * @brief Calculates the difference between two vectors.
 * 
 * @param u: The vector from which the difference points.
 * @param v: The vector to which the difference points.
 * @return The lenght of the difference vector |u - v|.
 */
__device__ double vec_dist_len(const double3 u, const double3 v)
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
__device__ double vec_dot(const double3 u, const double3 v)
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
__device__ double vec_dist_len_sq(const double3 u, const double3 v)
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
__device__ void vec_normalize(double3 v)
{
    double length = vec_lenght(v);

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
__device__ double3 vec_get_normal(const double3 u, const double3 v)
{
    double3 res;
    res.x = u.x - v.x;
    res.y = u.y - v.y;
    res.z = u.z - v.z;

    vec_normalize(res);

    return res;
}

/**
 * @brief Noramlizes a vector and returns it.
 * 
 * @param v: The vector that is to be normalized.
 * @return The normalized vector.
 */
__device__ double3 vec_get_normalized(double3 v)
{
    vec_normalize(v);

    return v;
}

/**
 * @brief Calculates the difference between two vectors.
 * 
 * @param u: The vector from which the difference will point.
 * @param v: The vector to which the difference will point.
 * @return The difference vector (u - v).
 */
__device__ double3 vec_diff(const double3 u, const double3 v)
{
    double3 res;
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
__device__ double3 vec_cross(const double3 u, const double3 v)
{
    double3 res;
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
__device__ void vec_set(double3 v, const double val)
{
    v.x = val;
    v.y = val;
    v.z = val;
}