/*
 * tetrahedron.h
 * This file defines the Tetrahedron structure.
 *
 *  Created on: Dec 17, 2012
 *      Author: lyang
 */

#ifndef TETRAHEDRON_H_
#define TETRAHEDRON_H_
#include <cmath>
#include "types.h"


#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif 
using namespace std;

class Point{
public:
	REAL x;
	REAL y;
	REAL z;
	CUDA_CALLABLE_MEMBER Point & operator=(const Point &rhs);
	CUDA_CALLABLE_MEMBER Point(const Point &point);
	CUDA_CALLABLE_MEMBER Point();
    CUDA_CALLABLE_MEMBER const Point operator+(const Point &other) const;
    CUDA_CALLABLE_MEMBER const Point operator-(const Point &other) const;
    CUDA_CALLABLE_MEMBER const Point operator*(const REAL &other) const;
    CUDA_CALLABLE_MEMBER const Point operator/(const REAL &other) const;
    CUDA_CALLABLE_MEMBER REAL dot(const Point &other) const;
    CUDA_CALLABLE_MEMBER const Point cross(const Point &other) const;
};

#endif /* TETRAHEDRON_H_ */
