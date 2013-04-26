/*
 * tetrahedron.cpp
 *
 *  Created on: Dec 20, 2012
 *      Author: lyang
 */
#include <cmath>
#include <algorithm>
#include <cstdio>

//#include "cuda_runtime.h"
//#include "device_launch_parameters.h"

using namespace std;

#include "point.h"
#include "types.h"
#define EPSILON 1e-6
#define EPSILON1 1e-11

CUDA_CALLABLE_MEMBER Point &  Point::operator=(const Point &rhs){
	x = rhs.x;
	y = rhs.y;
	z = rhs.z;
	return *this;
}

CUDA_CALLABLE_MEMBER Point::Point(const Point &point){
	x = point.x;
	y = point.y;
	z = point.z;
}

CUDA_CALLABLE_MEMBER Point::Point(){
	x = 0;
	y = 0;
	z = 0;
}


CUDA_CALLABLE_MEMBER const Point Point::operator+(const Point &other) const{
    Point result = *this;
    result.x += other.x;
    result.y += other.y;
    result.z += other.z;
    return result;
}

CUDA_CALLABLE_MEMBER const Point Point::operator-(const Point &other) const{
    Point result = *this;
    result.x -= other.x;
    result.y -= other.y;
    result.z -= other.z;
    return result;
}

CUDA_CALLABLE_MEMBER const Point Point::operator*(const REAL &other) const{
    Point result = *this;
    result.x *= other;
    result.y *= other;
    result.z *= other;
    return result;
}

CUDA_CALLABLE_MEMBER const Point Point::operator/(const REAL &other) const{
    Point result = *this;
    result.x /= other;
    result.y /= other;
    result.z /= other;
    return result;
}


CUDA_CALLABLE_MEMBER REAL Point::dot(const Point &other) const{
    Point result = *this;
    result.x *= other.x;
    result.y *= other.y;
    result.z *= other.z;
    return result.x + result.y + result.z;
}

CUDA_CALLABLE_MEMBER const Point Point::cross(const Point &b) const{
    Point a = *this;
    Point res;
    res.x = a.y * b.z - a.z * b.y;
    res.y = a.z * b.x - a.x * b.z;
    res.z = a.x * b.y - a.y * b.x;

    return res;
}

