/*
 * Projective_Point.hpp
 *
 *  Created on: 8 αιεπ 2014
 *      Author: Eyal
 */

#ifndef PROJECTIVEPOINT_HPP_
#define PROJECTIVEPOINT_HPP_

#include <NTL/GF2E.h>
#include <NTL/ZZ_p.h>
#include <mem.h>
#include "Constants.hpp"

using NTL::GF2E;

class ProjectivePoint;

class ProjectivePoint {
	GF2E x, y, z;
	static ProjectivePoint *INF;

public:
	static const int c = 1, d = 2;

	static void initInf();

	static ProjectivePoint& getINF();

	ProjectivePoint(const GF2E& _x, const GF2E& _y, const GF2E& _z = Constants::one) :
			x(_x), y(_y), z(_z) {

	}

	ProjectivePoint() :
			ProjectivePoint(GF2E(), GF2E(), GF2E()) {
	}

	bool isAffine() const;
	ProjectivePoint toAffine() const;

	bool isINF() const;

	ProjectivePoint operator-() const;
	ProjectivePoint neg() const;

	bool operator==(const ProjectivePoint& other) const;

	bool operator!=(const ProjectivePoint& other) const;

	const GF2E& getX() const;
	const GF2E& getY() const;
	const GF2E& getZ() const;

	void setX(const GF2E& x);
	void setY(const GF2E& y);
	void setZ(const GF2E& z);

};

std::ostream& operator<<(std::ostream& s, const ProjectivePoint& a);

#endif /* PROJECTIVE_POINT_HPP_ */
