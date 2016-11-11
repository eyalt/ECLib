/*
 * Elliptic_Curve.hpp
 *
 *  Created on: 8 αιεπ 2014
 *      Author: Eyal
 */

#ifndef ELLIPTICCURVE_HPP_
#define ELLIPTICCURVE_HPP_

#include <NTL/GF2E.h>
#include <NTL/ZZ.h>
#include "ProjectivePoint.hpp"

using NTL::GF2E;

class EllipticCurve {
	GF2E b; //OPT: delete a

	/**
	 * used only in Montgomery multiplication.
	 * x1,z1 coordinates of p1. x2,z2 coordinates of p2.
	 * computes x,z coordinate of p1 + p2 in x1,z1.
	 * temp is used for temp var.
	 */
	void mAdd(GF2E& x1, GF2E& z1, const GF2E& x2, const GF2E& z2,
			const ProjectivePoint& p, GF2E& temp) const;

	void mDouble(GF2E& x1, GF2E& z1, const ProjectivePoint& p,
			GF2E& temp) const;

public:
	EllipticCurve(const GF2E& _b) :
			b(_b) {
	}

	const GF2E& getB() const;

	ProjectivePoint pointAdd(const ProjectivePoint& p1,
			const ProjectivePoint& p2) const;

	ProjectivePoint pointDoubling(const ProjectivePoint& p) const;

	ProjectivePoint pointMulMontogmery(const ProjectivePoint& p,
			const NTL::ZZ& k) const;

	bool isOn(const ProjectivePoint& p) const;

	ProjectivePoint genPoint() const;
	GF2E jInvariant() const;
};
bool operator==(const EllipticCurve& e1, const EllipticCurve& e2);
bool operator!=(const EllipticCurve& e1, const EllipticCurve& e2);

#endif /* ELLIPTIC_CURVE_HPP_ */
