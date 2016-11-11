/*
 * Projective_Point.cpp
 *
 *  Created on: 8 αιεπ 2014
 *      Author: Eyal
 */

#include "ProjectivePoint.hpp"
#include "EC_Exceptions.hpp"
#include <NTL/GF2E.h>
#include <exception>
#include <assert.h>

#include <stdio.h>

using namespace NTL;

ProjectivePoint *ProjectivePoint::INF;

void ProjectivePoint::initInf() {
	if (INF) {
		delete INF;
	}
	INF = new ProjectivePoint(Constants::one, Constants::zero, Constants::zero);
}

ProjectivePoint& ProjectivePoint::getINF() {
	return *INF;
}

bool ProjectivePoint::isINF() const { //will get some false positives, but ok for what we need and MUCH faster then ==
	return IsZero(z);
}

bool ProjectivePoint::isAffine() const {
	return (bool) IsOne(z);
}

ProjectivePoint ProjectivePoint::toAffine() const {
	if (this->isINF()) {
		throw AffineException();
	}

	GF2E z_inv;
	inv(z_inv, z);

	return ProjectivePoint(x * z_inv, y * sqr(z_inv), GF2E::zero() + 1);
}

ProjectivePoint ProjectivePoint::operator-() const {
	return neg();
}

ProjectivePoint ProjectivePoint::neg() const {
	assert(c == 1 && d == 2);
	if (this->isINF()) {
		return getINF();
	}
	return ProjectivePoint(x, x + y, z);
}

bool ProjectivePoint::operator==(const ProjectivePoint& other) const {
	assert(c == 1 && d == 2);

	const GF2E &x1 = x, &y1 = y, &z1 = z;
	const GF2E &x2 = other.x, &y2 = other.y, &z2 = other.z;

	if (IsZero(z1) && !IsZero(z2)) {
		assert(!IsZero(x1) && IsZero(y1));
		return false;
	}
	if ((!IsZero(z1) && IsZero(z2))) {
		assert(!IsZero(x2) && IsZero(y2));
		return false;
	}

	if (IsZero(z1) && IsZero(z2)) {
		assert(!IsZero(x1) && IsZero(y1));
		assert(!IsZero(x2) && IsZero(y2));

		return true;
	}

	return (x1 * z2 == x2 * z1) && (y1 * sqr(z2) == y2 * sqr(z1));

	return true;
}

bool ProjectivePoint::operator!=(const ProjectivePoint& other) const {
	return !((*this) == other);
}

const GF2E& ProjectivePoint::getX() const {
	return x;
}

const GF2E& ProjectivePoint::getY() const {
	return y;
}

const GF2E& ProjectivePoint::getZ() const {
	return z;
}

void ProjectivePoint::setX(const GF2E& x) {
	this->x = x;
}

void ProjectivePoint::setY(const GF2E& y) {
	this->y = y;
}

void ProjectivePoint::setZ(const GF2E& z) {
	this->z = z;
}
std::ostream& operator<<(std::ostream& s, const ProjectivePoint& a) {
	return s << a.getX() << std::endl << a.getY() << std::endl << a.getZ();
}
