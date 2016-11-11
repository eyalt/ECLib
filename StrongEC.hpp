/*
 * CryptoEC.hpp
 *
 *  Created on: Jul 28, 2014
 *      Author: Itay
 */

#ifndef STRONGEC_HPP_
#define STRONGEC_HPP_

#include <NTL/GF2E.h>
#include <NTL/ZZ.h>
#include "EllipticCurve.hpp"

using NTL::GF2E;
using NTL::ZZ;

#ifndef MAX_ALLOWED_PRIMES
#define MAX_ALLOWED_PRIMES 10000
#endif

#define MILLER_RABIN_TRIALS 10

class StrongEC {
private:
	EllipticCurve curve;
	ZZ n;
	ZZ p;
	ZZ order; //order = n*p
	bool strong;

	ProjectivePoint generator; //a point of order p

	/***
	 * if the order is "big" (n <= MAX_ALLOWED) then stores order = n*p into n and p and returns true.
	 * else, it returns false.
	 */
	void calcStrength();

public:

	StrongEC(); //will build a new RANDOM strong EC using satoh

	StrongEC(const EllipticCurve& _curve) :
			curve(_curve) {
		countPointsSatoh();
		calcStrength();
		setGenerator();
	}

	ProjectivePoint pointAdd(const ProjectivePoint& p1,
			const ProjectivePoint& p2) const;
	ProjectivePoint pointDoubling(const ProjectivePoint& p) const;
	ProjectivePoint pointMul(const ProjectivePoint& p, const NTL::ZZ& k) const;
	bool isOn(const ProjectivePoint& p) const;

	const EllipticCurve& getCurve() const;
	const ZZ& getN() const;
	const ZZ& getP() const;
	const ZZ& getOrder() const;
	bool isStrong() const;

	bool hasGenerator() const;
	ProjectivePoint getGenerator() const; //returns a generator for the strong sub group
	void setGenerator(); //finds and sets a generator for the strong sub group
	void setGenerator(const ProjectivePoint& g); //finds and sets a generator for the strong sub group

	bool isStrongPoint(const ProjectivePoint& g) const;

	ProjectivePoint genStrongPoint() const; //different then the generator
	ProjectivePoint genPoint() const;

	GF2E jInvariant() const;

	void countPointsSatoh();
};

bool operator==(const StrongEC& e1, const StrongEC& e2);
bool operator!=(const StrongEC& e1, const StrongEC& e2);

#endif /* STRONGEC_HPP_ */
