/*
 * Elliptic_Curve.cpp
 *
 *  Created on: 8 αιεπ 2014
 *      Author: Eyal
 */

#include "EllipticCurve.hpp"
#include "EC_Exceptions.hpp"
#include <NTL/GF2E.h>
#include <exception>
#include <assert.h>

#include <stdio.h>

using namespace NTL;

void EllipticCurve::mAdd(GF2E& x1, GF2E& z1, const GF2E& x2, const GF2E& z2,
		const ProjectivePoint& p, GF2E& temp) const {
	const GF2E &p_x = p.getX();

	mul(temp, x1, z2);
	mul(z1, x2, z1);

	mul(x1, temp, z1);

	add(z1, temp, z1);
	sqr(z1, z1);

	mul(temp, p_x, z1);
	add(x1, temp, x1);
}

void EllipticCurve::mDouble(GF2E& x1, GF2E& z1, const ProjectivePoint& p,
		GF2E& temp) const {
	sqr(temp, z1);
	sqr(x1, x1);
	mul(z1, x1, temp);

	sqr(temp, temp);
	mul(temp, b, temp);
	sqr(x1, x1);
	add(x1, x1, temp);
}

ProjectivePoint EllipticCurve::pointMulMontogmery(const ProjectivePoint& p,
		const ZZ& k) const {
	if (p.isINF()) {
		return ProjectivePoint::getINF();
	}

	if (IsZero(k)) {
		return ProjectivePoint::getINF();
	}

//	ProjectivePoint res;
	GF2E res_x, res_y, res_z = Constants::one;
	const GF2E &p_x = p.getX(), &p_y = p.getY();

	GF2E x1(p_x), z1(GF2E::zero() + 1);
	GF2E z2(sqr(p_x));
	GF2E x2(sqr(z2));
	add(x2, x2, b);

	GF2E temp;

	for (int i = NTL::NumBits(k) - 2; i >= 0; i--) {
		if (NTL::bit(k, i) == 1) {
			//compute x1,z1
			mAdd(x1, z1, x2, z2, p, temp);

			//compute x2,z2
			mDouble(x2, z2, p, temp);

		} else {
			//compute x2,z2
			mAdd(x2, z2, x1, z1, p, temp);

			//compute x1,z1
			mDouble(x1, z1, p, temp);
		}
	}
	if (IsZero(z1)) {
		return ProjectivePoint(Constants::one, Constants::zero, Constants::zero);
	}

	if (IsZero(z2)) {
		res_x = p_x;
		res_y = p_x + p_y;
		return ProjectivePoint(res_x, res_y, res_z);
	}

	sqr(res_y, p_x);
	add(res_y, res_y, p_y);

	mul(temp, z1, z2);
	mul(res_y, res_y, temp);

	mul(z1, p_x, z1);
	add(z1, x1, z1);

	mul(z2, p_x, z2);
	mul(res_x, x1, z2);
	add(z2, x2, z2);

	mul(z1, z1, z2);
	add(res_y, z1, res_y);

	mul(temp, p_x, temp);
	inv(temp, temp);

	mul(res_y, temp, res_y);
	mul(res_x, res_x, temp);

	add(temp, p_x, res_x);
	mul(res_y, temp, res_y);
	add(res_y, res_y, p_y);

	return ProjectivePoint(res_x, res_y, res_z);
}

bool EllipticCurve::isOn(const ProjectivePoint& p) const {
	assert(ProjectivePoint::c == 1 && ProjectivePoint::d == 2);

	if (p.isINF()) {
		return true;
	}

	GF2E t1, t2, t3, t4;
	const GF2E &x = p.getX(), &y = p.getY(), &z = p.getZ();
	mul(t1, x, z);
	add(t3, y, t1);
	mul(t3, y, t3);

	sqr(t2, z);
	sqr(t4, x);
	mul(t1, t1, t4);

	sqr(t2, t2);
	mul(t2, t2, b);
	add(t1, t1, t2);

	return t1 == t3;
}

ProjectivePoint EllipticCurve::pointDoubling(const ProjectivePoint& p) const {
	assert(ProjectivePoint::c == 1 && ProjectivePoint::d == 2);
	assert(isOn(p));

	if (p.isINF()) {
		return p;
	}

	if (p == p.neg()) {
		return ProjectivePoint::getINF();
	}

	ProjectivePoint res;

	const GF2E &x1 = p.getX(), &y1 = p.getY(), &z1 = p.getZ();

	GF2E x3, y3, z3;
	GF2E t1, t2;

	sqr(t1, z1);
	sqr(t2, x1);
	mul(z3, t1, t2);
	sqr(x3, t2);
	sqr(t1, t1);
	mul(t2, t1, b);
	add(x3, x3, t2);
	sqr(t1, y1);

	add(t1, t1, t2);
	mul(y3, x3, t1);
	mul(t1, t2, z3);
	add(y3, y3, t1);

	return ProjectivePoint(x3, y3, z3);
}

ProjectivePoint EllipticCurve::pointAdd(const ProjectivePoint& p1,
		const ProjectivePoint& p2) const {
	assert(ProjectivePoint::c == 1 && ProjectivePoint::d == 2);
	assert(isOn(p1) && isOn(p2));

	if (!p2.isAffine()) {
		if (p1.isAffine()) {
			return pointAdd(p2, p1);
		} else {
			throw AffineException();
		}
	}

	if (p1 == p2.neg()) { //Probably non needed cause we check this in the next part
		return ProjectivePoint::getINF();
	}
	if (p1.isINF()) {
		return p2;
	}
	if (p2.isINF()) {
		return p1;
	}

	GF2E t1, t2, t3;

	const GF2E &x1 = p1.getX(), &y1 = p1.getY(), &z1 = p1.getZ();
	const GF2E &x2 = p2.getX(), &y2 = p2.getY();
	GF2E x3, y3, z3;

	mul(t1, z1, x2);
	sqr(t2, z1);
	add(x3, x1, t1);
	mul(t1, z1, x3);
	mul(t3, t2, y2);
	add(y3, y1, t3);

	if (IsZero(x3)) {
		if (IsZero(y3)) {
			return pointDoubling(p2);
		}
		return ProjectivePoint::getINF();
	}

	sqr(z3, t1);
	mul(t3, t1, y3);

	sqr(t2, x3);
	mul(x3, t2, t1);
	sqr(t2, y3);
	add(x3, x3, t2);
	add(x3, x3, t3);
	mul(t2, x2, z3);
	add(t2, t2, x3);
	sqr(t1, z3);
	add(t3, t3, z3);
	mul(y3, t3, t2);
	add(t2, x2, y2);
	mul(t3, t1, t2);
	add(y3, y3, t3);

	return ProjectivePoint(x3, y3, z3);

}

const GF2E& EllipticCurve::getB() const {
	return b;
}

GF2E EllipticCurve::jInvariant() const {
	return inv(b);
}

GF2E sqrt(const GF2E& x) {
	int n = GF2E::degree();
	GF2E temp = x;
	for (int i = 0; i < n - 1; i++) {
		sqr(temp, temp);
	}
	return temp;
}

bool solveQuadEqOverGF2E(GF2E& z, const GF2E& betta) {
	if (GF2E::degree() % 2 == 1) {
		z = betta;
		for (int i = 1; i <= (GF2E::degree() - 1) / 2; i++) {
			NTL::sqr(z, z);
			NTL::sqr(z, z);
			NTL::add(z, z, betta);
		}
	} else {
		z = GF2E::zero();
		GF2E w, p;
		do {
			NTL::random(p);
			w = p;
			for (int i = 1; i <= GF2E::degree() - 1; i++) {
				sqr(z, z);
				sqr(w, w);
				add(z, z, w * betta);
				add(w, w, p);
			}
		} while (NTL::IsZero(w));
	}

	if (power(z, 2) + z == betta) {
		return true;
	}
	return false;
}

ProjectivePoint EllipticCurve::genPoint() const {
	GF2E x;
	bool found = false;
	GF2E z, alpha, betta, t1;

	while (!found) {
		random(x);
		if (IsZero(x)) {
			return ProjectivePoint(Constants::zero, sqrt(b), Constants::one);
		}

		power(alpha, x, 3);
		add(alpha, alpha, b);
		if (IsZero(alpha)) {
			return ProjectivePoint(x, Constants::zero, Constants::one);
		}

		mul(betta, power(x, -2), alpha);

		found = solveQuadEqOverGF2E(z, betta);
	}
	add(z, z, rand() % 2);
	mul(z, z, x);
	return ProjectivePoint(x, z, Constants::one);

}

bool operator==(const EllipticCurve& e1, const EllipticCurve& e2) {
	return e1.getB() == e2.getB();
}
bool operator!=(const EllipticCurve& e1, const EllipticCurve& e2) {
	return !(e1 == e2);
}
