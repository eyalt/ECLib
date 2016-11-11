#include <NTL/tools.h>
#include <NTL/GF2E.h>
#include <NTL/ZZ.h>
#include <NTL/GF2XFactoring.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pE.h>

#include <string>
#include <sstream>
#include <iostream>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <time.h>

#include <EC_Exceptions.hpp>
#include <ProjectivePoint.hpp>
#include <EllipticCurve.hpp>
#include <StrongEC.hpp>
#include <CryptoLib.hpp>

using NTL::GF2E;
using NTL::GF2X;
using NTL::ZZ;

using std::cout;
using std::endl;

//HELPER FUNCTION DECLERATIONS:
bool checkIsOn(GF2E x, GF2E y, EllipticCurve e);

ProjectivePoint sillyAdd(const ProjectivePoint& p1, const ProjectivePoint& p2,
		const EllipticCurve& e);

ProjectivePoint sillyDouble(const ProjectivePoint& p1, const EllipticCurve& e);

ProjectivePoint sillyMul(const ProjectivePoint& p, const ZZ& k,
		const EllipticCurve& e);

GF2E sillySqrt(const GF2E& x);

ProjectivePoint sillyGenPoint(const EllipticCurve& e);

bool sillySolveQuadEqOverGF2E(GF2E& z, const GF2E& betta);

//TEST FUNCTION DECLERATIONS:
bool testEqualityINF();

bool testEquality();

bool testIsOn();

bool testToAffine();

bool testGenPoint();

bool testDouble();

bool testAdd();

bool testMul();

bool testMulRand();

bool testMulTime();

bool testSatoh();

bool testFindStrongEC();

bool testECDH();

void statSatoh();

//HELPER FUNCTIONS:
bool checkIsOn(GF2E x, GF2E y, EllipticCurve e) {

	ProjectivePoint p = ProjectivePoint(x, y, GF2E::zero() + 1);
	if (!e.isOn(p)) {
		return false;
	}

	for (int i = 0; i < 10; i++) {
		NTL::GF2E z;
		do {
			NTL::random(z);
		} while (NTL::IsZero(z));

		p = ProjectivePoint(x * z, y * z * z, z);
		if (!e.isOn(p)) {
			return false;
		}
	}

	return true;
}

ProjectivePoint sillyAdd(const ProjectivePoint& p1, const ProjectivePoint& p2,
		const EllipticCurve& e) {
	if (p1 == ProjectivePoint::getINF()) {
		return p2;
	}
	if (p2 == ProjectivePoint::getINF()) {
		return p1;
	}
	if (p1 == p2) {
		return sillyDouble(p1, e);
	}
	if (p1 == -p2) {
		return ProjectivePoint::getINF();
	}
	GF2E lambda = (p1.getY() + p2.getY()) / (p1.getX() + p2.getX());
	GF2E x3 = lambda * lambda + lambda + p1.getX() + p2.getX();
	GF2E y3 = lambda * (p1.getX() + x3) + x3 + p1.getY();
	return ProjectivePoint(x3, y3, GF2E::zero() + 1);
}

ProjectivePoint sillyDouble(const ProjectivePoint& p1, const EllipticCurve& e) {
	if (p1 == -p1) {
		return ProjectivePoint::getINF();
	}
	GF2E lambda = p1.getX() + p1.getY() / p1.getX();
	GF2E x3 = lambda * lambda + lambda;
	GF2E y3 = p1.getX() * p1.getX() + lambda * x3 + x3;

	return ProjectivePoint(x3, y3, GF2E::zero() + 1);
}

ProjectivePoint sillyMul(const ProjectivePoint& p, const ZZ& k,
		const EllipticCurve& e) {
	if (NTL::IsZero(k)) {
		return ProjectivePoint::getINF();
	}
	ProjectivePoint res = p;
	for (ZZ i = NTL::conv<ZZ>(1); i < k; i++) {
		res = sillyAdd(res, p, e);
	}
	return res;
}

GF2E sillySqrt(const GF2E& x) {
	int n = GF2E::degree();
	GF2E temp = x;
	for (int i = 0; i < n - 1; i++) {
		sqr(temp, temp);
	}
	return temp;
}

bool sillySolveQuadEqOverGF2E(GF2E& z, const GF2E& betta) {
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

ProjectivePoint sillyGenPoint(const EllipticCurve& e) {
	GF2E x;
	bool found = false;
	GF2E z, alpha, betta, t1;
	GF2E y = sillySqrt(e.getB());
	int i = 1;
	while (!found) {
		random(x);
		if (IsZero(x)) {
			return ProjectivePoint(Constants::zero, y, Constants::one);
		}

		power(alpha, x, 3);
		add(alpha, alpha, e.getB());
		if (IsZero(alpha)) {
			return ProjectivePoint(x, Constants::zero, Constants::one);
		}

		mul(betta, power(x, -2), alpha);

		found = sillySolveQuadEqOverGF2E(z, betta);
		i++;
	}
	add(z, z, rand() % 2);
	mul(z, z, x);
	assert(e.isOn(ProjectivePoint(x, z, Constants::one)));
	return ProjectivePoint(x, z, Constants::one);

}

//TEST FUNCTIONS:
bool testEqualityINF() {
	printf("TESTING testEqualityINF:\n");
	NTL::GF2X P;
	NTL::BuildIrred(P, 11);
	GF2E::init(P);
	Constants::initConstants();

	printf("Checking inf==inf CASE1...");
	if (ProjectivePoint::getINF() != ProjectivePoint::getINF()) {
		printf("FAILED\n");
		return false;
	}
	printf("OK\n");

	printf("Checking inf==inf CASE2 ...");
	for (int i = 0; i < 10; i++) {
		NTL::GF2E x;
		NTL::random(x);
		ProjectivePoint p = ProjectivePoint(x, GF2E::zero(), GF2E::zero());
		if (ProjectivePoint::getINF() != p) {
			printf("FAILED\n");
			return false;
		}
		if (p != ProjectivePoint::getINF()) {
			printf("FAILED\n");
			return false;
		}
	}
	printf("OK\n");

	printf("Checking inf==inf CASE3 ...");
	for (int i = 0; i < 10; i++) {
		NTL::GF2E x, y;
		do {
			NTL::random(x);
		} while (NTL::IsZero(x));
		do {
			NTL::random(y);
		} while (NTL::IsZero(y));

		ProjectivePoint p1 = ProjectivePoint(x, GF2E::zero(), GF2E::zero());
		ProjectivePoint p2 = ProjectivePoint(y, GF2E::zero(), GF2E::zero());
		if (p1 != p2) {
			printf("FAILED\n");
			return false;
		}
		if (p2 != p1) {
			printf("FAILED\n");
			return false;
		}
	}
	printf("OK\n");

	printf("Checking inf!=point for any affine point ...");
	for (int i = 0; i < 10; i++) {
		NTL::GF2E x, y;
		NTL::random(x);
		NTL::random(y);

		ProjectivePoint p = ProjectivePoint(x, y, GF2E::zero() + 1);
		if (ProjectivePoint::getINF() == p) {
			printf("FAILED\n");
			return false;
		}
	}
	printf("OK\n");

	printf("Checking inf!=point for any projective point CASE1 ...");
	for (int i = 0; i < 10; i++) {
		NTL::GF2E x, y, z;
		NTL::random(x);
		NTL::random(y);

		do {
			NTL::random(z);
		} while (NTL::IsZero(z));

		ProjectivePoint p = ProjectivePoint(x, y, z);
		if (ProjectivePoint::getINF() == p) {
			printf("FAILED\n");
			return false;
		}
		if (p == ProjectivePoint::getINF()) {
			printf("FAILED\n");
			return false;
		}
	}
	printf("OK\n");

	printf("Checking inf!=point for any projective point CASE2 ...");
	for (int i = 0; i < 10; i++) {
		NTL::GF2E x, y, z;
		NTL::random(x);
		NTL::random(y);
		do {
			NTL::random(z);
		} while (NTL::IsZero(z));

		NTL::GF2E xINF;
		do {
			NTL::random(xINF);
		} while (NTL::IsZero(xINF));

		ProjectivePoint p1 = ProjectivePoint(xINF, GF2E::zero(), GF2E::zero());
		ProjectivePoint p2 = ProjectivePoint(x, y, z);
		if (p1 == p2) {
			printf("FAILED\n");
			return false;
		}
		if (p2 == p1) {
			printf("FAILED\n");
			return false;
		}
	}
	printf("OK\n");

	printf("\n");
	return true;

}

bool testEquality() {
	printf("TESTING testEquality:\n");
	NTL::GF2X P;
	NTL::BuildIrred(P, 11);
	GF2E::init(P);
	Constants::initConstants();

	printf("Checking point==point for any affine point ...");
	for (int i = 0; i < 10; i++) {
		NTL::GF2E x, y;
		NTL::random(x);
		NTL::random(y);
		ProjectivePoint affine = ProjectivePoint(x, y, GF2E::zero() + 1);
		if (affine != affine) {
			printf("FAILED\n");
			return false;
		}
	}
	printf("OK\n");

	printf("Checking point1==point2 for affine+projective ...");
	for (int i = 0; i < 10; i++) {
		NTL::GF2E x, y, z;
		NTL::random(x);
		NTL::random(y);
		do {
			NTL::random(z);
		} while (NTL::IsZero(z));

		ProjectivePoint affine = ProjectivePoint(x, y, GF2E::zero() + 1);
		ProjectivePoint proj = ProjectivePoint(x * z, y * z * z, z);
		if (affine != proj) {
			printf("FAILED\n");
			return false;
		}
		if (proj != affine) {
			printf("FAILED\n");
			return false;
		}
	}
	printf("OK\n");

	printf("Checking point1==point2 for projective+projective ...");
	for (int i = 0; i < 10; i++) {
		NTL::GF2E x, y, z1, z2;
		NTL::random(x);
		NTL::random(y);
		do {
			NTL::random(z1);
		} while (NTL::IsZero(z1));
		do {
			NTL::random(z2);
		} while (NTL::IsZero(z2));

		ProjectivePoint proj1 = ProjectivePoint(x * z1, y * z1 * z1, z1);
		ProjectivePoint proj2 = ProjectivePoint(x * z2, y * z2 * z2, z2);
		if (proj1 != proj2) {
			printf("FAILED\n");
			return false;
		}
		if (proj2 != proj1) {
			printf("FAILED\n");
			return false;
		}
	}
	printf("OK\n");

	printf("\n");
	return true;
}

bool testIsOn() {
	printf("TESTING testIsOn:\n");

	std::stringstream ss;

	NTL::GF2X P;
	ss.clear();
	ss
			<< "0x800000000000000000000000000000000000000000000000000000000000000000010a1";
	ss >> P;
	P.normalize();
	NTL::reverse(P, P);
	GF2E::init(P);
	Constants::initConstants();

	GF2E x, y;

	cout
			<< "Randomizing x,y to calculate b, and checking they are on the curve...";
	for (int i = 0; i < 10; i++) {
		do {
			NTL::random(x);
		} while (NTL::IsZero(x));
		do {
			NTL::random(y);
		} while (NTL::IsZero(y));

		ProjectivePoint p = ProjectivePoint(x, y, Constants::one);

		GF2E t1, t2, t3;
		sqr(t1, y);
		mul(t2, x, y);
		sqr(t3, x);
		mul(t3, t3, x);
		add(t1, t1, t2);
		add(t1, t1, t3);

		EllipticCurve e = EllipticCurve(t1);

		if (!e.isOn(p)) {
			cout << "Error for x=" << x << " y=" << y << endl;
			return false;
		}
	}
	cout << "OK" << endl;
	cout << endl;

	return true;
}

bool testGenPoint() {
	std::cout << "TESTING testGenPoint:" << std::endl;

	NTL::GF2X P;
	NTL::BuildIrred(P, 11);
	GF2E::init(P);
	Constants::initConstants();

	GF2E b;
	std::cout
			<< "Checking order 11: Checking point generation for random EC...";
	for (int i = 0; i < 10; i++) {
		NTL::random(b);
		EllipticCurve e(b);

		ProjectivePoint g;
		for (int j = 0; j < 10; j++) {
			g = e.genPoint();
			if (!e.isOn(g)) {
				std::cout << "Failed" << std::endl;
			}
		}
	}
	std::cout << "OK" << std::endl;

	NTL::BuildIrred(P, 120);
	GF2E::init(P);
	Constants::initConstants();

	std::cout
			<< "Checking order 120: Checking point generation for random EC...";
	for (int i = 0; i < 10; i++) {
		NTL::random(b);
		EllipticCurve e(b);

		ProjectivePoint g;
		for (int j = 0; j < 10; j++) {
			g = e.genPoint();
			if (!e.isOn(g)) {
				std::cout << "Failed" << std::endl;
			}
		}
	}
	std::cout << "OK" << std::endl;
	cout << endl;
	return true;
}

bool testDouble() {
	printf("TESTING testDouble:\n");

	NTL::GF2X P;
	NTL::BuildIrred(P, 11);
	GF2E::init(P);
	Constants::initConstants();

	GF2E b;
	NTL::random(b);
	while (NTL::IsZero(b)) {
		NTL::random(b);
	}
	EllipticCurve e = EllipticCurve(b);

	printf("Checking that double(inf)=inf CASE1...");
	ProjectivePoint p = ProjectivePoint::getINF();
	if (e.pointDoubling(p) != ProjectivePoint::getINF()) {
		printf("FAILED\n");
		return false;
	}
	printf("OK\n");

	printf("Checking that double(inf)=inf CASE2...");
	for (int i = 0; i < 10; i++) {
		NTL::GF2E x;
		NTL::random(x);
		if (NTL::IsZero(x)) {
			continue;
		}
		ProjectivePoint p = ProjectivePoint(x, GF2E::zero(), GF2E::zero());
		if (e.pointDoubling(p) != ProjectivePoint::getINF()) {
			printf("FAILED\n");
			return false;
		}
	}
	printf("OK\n");

	printf(
			"Checking on random points that our double is the same as the regular one ...");
	for (int i = 0; i < 50; i++) {
		ProjectivePoint pAffine = sillyGenPoint(e);
		//Sanity check:
		ProjectivePoint expected = sillyDouble(pAffine, e);
		ProjectivePoint result = e.pointDoubling(pAffine);
		if (expected != result) {
			printf("FAILED\n");
			std::cout << "Expected result: (" << expected.getX() << " , "
					<< expected.getY() << ")" << std::endl;
			std::cout << "Expected result: (" << result.toAffine().getX()
					<< " , " << result.toAffine().getY() << ")" << std::endl;
			std::cout << "Expected result: (" << result.getX() << " , "
					<< result.getY() << " , " << result.getZ() << ")"
					<< std::endl;
			return false;
		}

		//Now for projective
		for (int i = 0; i < 5; i++) {
			GF2E z;
			NTL::random(z);
			if (NTL::IsZero(z)) {
				continue;
			}
			ProjectivePoint pProj = ProjectivePoint(z * pAffine.getX(),
					z * z * pAffine.getY(), z);
			expected = sillyDouble(pAffine, e);
			result = e.pointDoubling(pProj);
			if (expected != result) {
				printf("FAILED\n");
				std::cout << "Expected result: (" << expected.getX() << " , "
						<< expected.getY() << ")" << std::endl;
				std::cout << "Expected result: (" << result.toAffine().getX()
						<< " , " << result.toAffine().getY() << ")"
						<< std::endl;
				return false;
			}
		}
	}
	printf("OK\n");

	printf("\n");
	return true;
}

bool testToAffine() {
	printf("TESTING testToAffine:\n");
	NTL::GF2X P;
	NTL::BuildIrred(P, 11);
	GF2E::init(P);
	Constants::initConstants();

	printf("Checking toAffine(inf) throws CASE1...");
	try {
		ProjectivePoint::getINF().toAffine();
		printf("FAILED\n");
		return false;
	} catch (AffineException& e) {
		printf("OK\n");
	}

	printf("Checking toAffine(inf) throws CASE2 ...");
	for (int i = 0; i < 10; i++) {
		NTL::GF2E x;
		NTL::random(x);
		ProjectivePoint p = ProjectivePoint(x, GF2E::zero(), GF2E::zero());
		try {
			p.toAffine();
			printf("FAILED\n");
			return false;
		} catch (AffineException& e) {
		}
	}
	printf("OK\n");

	printf("Checking toAffine(point) works on affine ...");
	for (int i = 0; i < 10; i++) {
		NTL::GF2E x, y;
		NTL::random(x);
		NTL::random(y);
		ProjectivePoint p = ProjectivePoint(x, y, GF2E::zero() + 1);
		if (p.toAffine().getX() != p.getX()
				|| p.toAffine().getY() != p.getY()) {
			printf("FAILED\n");
			return false;
		}
	}
	printf("OK\n");

	printf("Checking toAffine(point) works on projective ...");
	for (int i = 0; i < 10; i++) {
		NTL::GF2E x, y;
		NTL::random(x);
		NTL::random(y);
		ProjectivePoint affine = ProjectivePoint(x, y, GF2E::zero() + 1);
		for (int i = 0; i < 5; i++) {
			GF2E z;
			NTL::random(z);
			if (NTL::IsZero(z)) {
				continue;
			}
			ProjectivePoint proj = ProjectivePoint(x * z, y * z * z, z);
			if (proj.toAffine().getX() != affine.getX()
					|| proj.toAffine().getY() != affine.getY()) {
				printf("FAILED\n");
				return false;
			}
		}
	}
	printf("OK\n");

	printf("\n");
	return true;

}

bool testAdd() {
	printf("TESTING testAdd:\n");

	NTL::GF2X P;
	NTL::BuildIrred(P, 11);
	GF2E::init(P);
	Constants::initConstants();

	GF2E b;
	NTL::random(b);
	while (NTL::IsZero(b)) {
		NTL::random(b);
	}
	EllipticCurve e = EllipticCurve(b);

	printf("Checking that add(p1,p2) throws when p2 is proj CASE1...");
	for (int i = 0; i < 10; i++) {
		ProjectivePoint p1start = sillyGenPoint(e);
		GF2E x, z;
		NTL::random(x);
		NTL::random(z);
		if (NTL::IsOne(z) || NTL::IsZero(z)) {
			continue;
		}
		ProjectivePoint p1 = ProjectivePoint(p1start.getX() * z,
				p1start.getY() * z * z, z);
		ProjectivePoint p2 = ProjectivePoint(x, GF2E::zero(), GF2E::zero());
		try {
			e.pointAdd(p1, p2);
			printf("FAILED\n");
			return false;
		} catch (AffineException& e) {
		}
	}
	printf("OK\n");

	printf("Checking that add(p1,p2) throws when p2 is proj CASE2...");
	for (int i = 0; i < 10; i++) {
		ProjectivePoint p1start = sillyGenPoint(e);
		GF2E z1, z2;
		NTL::random(z2);
		NTL::random(z1);
		if (NTL::IsOne(z1) || NTL::IsZero(z1) || NTL::IsOne(z2)
				|| NTL::IsZero(z2)) {
			continue;
		}
		ProjectivePoint p1 = ProjectivePoint(p1start.getX() * z1,
				p1start.getY() * z1 * z1, z1);
		ProjectivePoint p2temp = sillyGenPoint(e);

		ProjectivePoint p2 = ProjectivePoint(p2temp.getX() * z2,
				p2temp.getY() * z2 * z2, z2);
		try {
			e.pointAdd(p1, p2);
			printf("FAILED\n");
			return false;
		} catch (AffineException& e) {
		}
	}
	printf("OK\n");

	printf("Checking that add(p1,p1) = double for random points ...");
	for (int i = 0; i < 50; i++) {
		ProjectivePoint pAffine = sillyGenPoint(e);
		//Sanity check:
		ProjectivePoint expected = sillyDouble(pAffine, e);
		ProjectivePoint result = e.pointAdd(pAffine, pAffine);
		if (expected != result) {
			printf("FAILED\n");
			std::cout << "Expected result: (" << expected.getX() << " , "
					<< expected.getY() << ")" << std::endl;
			std::cout << "Expected result: (" << result.toAffine().getX()
					<< " , " << result.toAffine().getY() << ")" << std::endl;
			return false;
		}

		//Now for projective
		for (int i = 0; i < 5; i++) {
			GF2E z;
			NTL::random(z);
			if (NTL::IsZero(z)) {
				continue;
			}
			ProjectivePoint pProj = ProjectivePoint(z * pAffine.getX(),
					z * z * pAffine.getY(), z);
			expected = sillyDouble(pAffine, e);
			result = e.pointAdd(pProj, pAffine);
			if (expected != result) {
				printf("FAILED\n");
				std::cout << "Expected result: (" << expected.getX() << " , "
						<< expected.getY() << ")" << std::endl;
				std::cout << "Actual result: (" << result.toAffine().getX()
						<< " , " << result.toAffine().getY() << ")"
						<< std::endl;
				return false;
			}
		}
	}
	printf("OK\n");

	printf("Checking that add(p1,p2) = sillyAdd for random points ...");
	for (int i = 0; i < 50; i++) {
		ProjectivePoint pAffine1 = sillyGenPoint(e);
		ProjectivePoint pAffine2 = sillyGenPoint(e);
		//Sanity check:
		ProjectivePoint expected = sillyAdd(pAffine1, pAffine2, e);
		ProjectivePoint result = e.pointAdd(pAffine1, pAffine2);
		if (expected != result) {
			printf("FAILED\n");
			std::cout << "P1: (" << pAffine1.getX() << " , " << pAffine1.getY()
					<< ")" << std::endl;
			std::cout << "P2: (" << pAffine2.getX() << " , " << pAffine2.getY()
					<< ")" << std::endl;
			std::cout << "Expected result: (" << expected.getX() << " , "
					<< expected.getY() << ")" << std::endl;
			std::cout << "Actual result: (" << result.toAffine().getX() << " , "
					<< result.toAffine().getY() << ")" << std::endl;
			return false;
		}

		//Now for projective
		for (int i = 0; i < 5; i++) {
			GF2E z;
			NTL::random(z);
			if (NTL::IsZero(z)) {
				continue;
			}
			ProjectivePoint pProj = ProjectivePoint(z * pAffine1.getX(),
					z * z * pAffine1.getY(), z);
			expected = sillyAdd(pAffine1, pAffine2, e);
			result = e.pointAdd(pProj, pAffine2);
			if (expected != result) {
				printf("FAILED\n");
				std::cout << "Expected result: (" << expected.getX() << " , "
						<< expected.getY() << ")" << std::endl;
				std::cout << "Actual result: (" << result.toAffine().getX()
						<< " , " << result.toAffine().getY() << ")"
						<< std::endl;
				return false;
			}
		}
	}
	printf("OK\n");

	printf("\n");
	return true;
}

bool testMul() {
	cout << "Testing Mul:" << endl;
	std::stringstream ss;

	NTL::GF2X P;

	ss << "[1 0 0 1 0 0 0 0 0 0 1]";
	ss >> P;
	GF2E::init(P);
	Constants::initConstants();

	ss << "[0 0 0 1 1 1 1 1 1]";
	ss >> P;
	EllipticCurve e(NTL::conv<GF2E>(P));

	ss << "[1 1 0 1 1 1 1 1]";
	ss >> P;
	ProjectivePoint gen(GF2E::zero(), NTL::conv<GF2E>(P), GF2E::zero() + 1);

	cout << "Checking that mul(([] , [1 1 0 1 1 1 1 1]), i) is correct:...";
	for (int i = 1; i <= 10; i++) {
		ProjectivePoint expected, result;
		expected = sillyMul(gen, NTL::conv<ZZ>(i), e);
		result = e.pointMulMontogmery(gen, NTL::conv<ZZ>(i));
		if (result != expected) {
			cout << "FAILED" << endl;
			std::cout << "i: " << i << std::endl;
			std::cout << "Expected result: (" << expected.getX() << " , "
					<< expected.getY() << "," << expected.getZ() << ")"
					<< std::endl;
			std::cout << "Actual result: (" << result.getX() << " , "
					<< result.getY() << "," << result.getZ() << ")"
					<< std::endl;
		}
	}
	cout << "OK" << endl;

	GF2E x, y;
	ss << "[0 0 0 1 1 0 0 0 0 1]";
	ss >> x;

	ss << "[1 1 1 0 0 1 0 0 1]";
	ss >> y;
	gen = ProjectivePoint(x, y, GF2E::zero() + 1);

	cout
			<< "Checking that mul(([0 0 0 1 1 0 0 0 0 1] , [1 1 1 0 0 1 0 0 1]), i) is correct:...";
	for (int i = 1; i <= 1024; i++) {
		ProjectivePoint expected, result;
		expected = sillyMul(gen, NTL::conv<ZZ>(i), e);
		result = e.pointMulMontogmery(gen, NTL::conv<ZZ>(i));
		if (result != expected) {
			cout << "FAILED" << endl;
			std::cout << "i: " << i << std::endl;
			std::cout << "Expected result: (" << expected.getX() << " , "
					<< expected.getY() << "," << expected.getZ() << ")"
					<< std::endl;
			std::cout << "Actual result: (" << result.getX() << " , "
					<< result.getY() << "," << result.getZ() << ")"
					<< std::endl;
		}
	}
	cout << "OK" << endl << endl;

	return true;
}

bool testMulRand() {
	printf("TESTING testMulRand:\n");

	NTL::GF2X P;
	NTL::BuildIrred(P, 11);
	GF2E::init(P);
	Constants::initConstants();

	GF2E b;
	NTL::random(b);

	while (NTL::IsZero(b)) {
		NTL::random(b);
	}
	EllipticCurve e = EllipticCurve(b);

	printf("Checking that mul(p,0)=inf for all p ...");
	for (int i = 0; i < 10; i++) {
		ProjectivePoint affine = sillyGenPoint(e);

		//Sanity check:
		ProjectivePoint result = e.pointMulMontogmery(affine, ZZ::zero());
		ProjectivePoint expected = sillyMul(affine, ZZ::zero(), e);
		if (result != expected) {
			std::cout << "P: (" << affine.getX() << " , " << affine.getY()
					<< ")" << std::endl;
			std::cout << "Expected result: (" << expected.getX() << " , "
					<< expected.getY() << "," << expected.getZ() << ")"
					<< std::endl;
			std::cout << "Actual result: (" << result.toAffine().getX() << " , "
					<< result.toAffine().getY() << "," << result.getZ() << ")"
					<< std::endl;
			return false;

		}

		//And for projective:
		for (int j = 0; j < 5; j++) {
			GF2E z;
			NTL::random(z);
			if (NTL::IsZero(z)) {
				continue;
			}
			ProjectivePoint proj = ProjectivePoint(affine.getX() * z,
					affine.getY() * z * z, z);
			ProjectivePoint result = e.pointMulMontogmery(proj, ZZ::zero());
			ProjectivePoint expected = sillyMul(affine, ZZ::zero(), e);
			if (result != expected) {
				std::cout << "P: (" << affine.getX() << " , " << affine.getY()
						<< ")" << std::endl;
				std::cout << "Expected result: (" << expected.getX() << " , "
						<< expected.getY() << "," << expected.getZ() << ")"
						<< std::endl;
				std::cout << "Actual result: (" << result.toAffine().getX()
						<< " , " << result.toAffine().getY() << ","
						<< result.getZ() << ")" << std::endl;
				return false;

			}
		}

	}
	printf("OK\n");

	printf("Checking that mul(inf,k)=inf for all k ...");
	for (int i = 0; i < 10; i++) {
		ZZ k;
		k = NTL::conv<ZZ>(i);

		//Sanity check:
		ProjectivePoint result = e.pointMulMontogmery(ProjectivePoint::getINF(),
				k);
		ProjectivePoint expected = sillyMul(ProjectivePoint::getINF(), k, e);
		if (result != expected) {
			std::cout << "K: " << k << std::endl;
			std::cout << "Expected result: (" << expected.getX() << " , "
					<< expected.getY() << "," << expected.getZ() << ")"
					<< std::endl;
			std::cout << "Actual result: (" << result.getX() << " , "
					<< result.getY() << "," << result.getZ() << ")"
					<< std::endl;
			return false;

		}

		//And for projective:
		for (int j = 0; j < 5; j++) {
			GF2E x;
			NTL::random(x);
			if (NTL::IsZero(x)) {
				continue;
			}
			ProjectivePoint proj = ProjectivePoint(x, GF2E::zero(),
					GF2E::zero());
			ProjectivePoint result = e.pointMulMontogmery(proj, k);
			ProjectivePoint expected = sillyMul(proj, k, e);
			if (result != expected) {
				std::cout << "K: " << k << std::endl;
				std::cout << "Expected result: (" << expected.getX() << " , "
						<< expected.getY() << "," << expected.getZ() << ")"
						<< std::endl;
				std::cout << "Actual result: (" << result.toAffine().getX()
						<< " , " << result.toAffine().getY() << ","
						<< result.getZ() << ")" << std::endl;
				return false;

			}
		}

	}
	printf("OK\n");

	cout << "Checking mul(p,k) for random p,k ...";
	for (int i = 0; i < 10; i++) {
		ProjectivePoint point = e.genPoint();
		ZZ k = NTL::RandomBnd(ZZ(INIT_VAL, "50"));

		ProjectivePoint result = e.pointMulMontogmery(point, k);
		ProjectivePoint expected = sillyMul(point, k, e);
		if (result != expected) {
			std::cout << "K: " << k << std::endl;
			std::cout << "Expected result: (" << expected.getX() << " , "
					<< expected.getY() << "," << expected.getZ() << ")"
					<< std::endl;
			std::cout << "Actual result: (" << result.toAffine().getX() << " , "
					<< result.toAffine().getY() << "," << result.getZ() << ")"
					<< std::endl;
			return false;

		}

	}
	cout << "OK" << endl;
	printf("\n");
	return true;
}

bool testSatoh() {
	cout << "Testing Satoh point counting:" << endl;

	NTL::GF2X P;
	std::stringstream ss;
	NTL::SetSeed(ZZ(INIT_VAL, clock()));

	NTL::SetCoeff(P, 11);
	NTL::SetCoeff(P, 2);
	NTL::SetCoeff(P, 0);
	GF2E::init(P);
	Constants::initConstants();

	GF2E b;
	cout << "Checking for 11 bits..." << std::flush;
	ss.clear();
	ss << "[1 1 1 1 1 1 1 1 1 1 1]";
	ss >> b;

	StrongEC e = StrongEC(EllipticCurve(b));
	if (e.getOrder() != ZZ(INIT_VAL, "2088")) {
		cout << "Failed" << endl;
		return false;
	}
	cout << "OK" << endl;

	P = GF2X::zero();

	NTL::SetCoeff(P, 283);
	NTL::SetCoeff(P, 12);
	NTL::SetCoeff(P, 7);
	NTL::SetCoeff(P, 5);
	NTL::SetCoeff(P, 0);
	GF2E::init(P);
	Constants::initConstants();

	b = GF2E::zero();

	cout << "Checking for 283 bits. #1..." << std::flush;

	ss.clear();
	ss
			<< "[0 1 0 0 1 1 0 1 0 1 0 1 0 0 0 0 1 1 1 0 1 1 0 0 1 1 1 0 1 0 1 0 0 1 0 1 0 0 0 1 0 0 0 0 1 1 1 0 1 1 0 1 1 1 1 0 1 0 0 0 1 0 1 1 0 0 1 1 0 0 0 0 1 0 0 1 1 0 1 0 1 1 0 0 1 1 1 0 0 1 0 0 0 1 1 0 1 1 1 1 1 0 1 0 1 1 0 1 0 1 0 1 1 0 1 1 0 1 1 0 1 1 1 0 1 0 1 0 0 1 0 0 0 0 1 1 1 1 0 0 0 1 1 1 1 1 1 1 1 1 0 1 1 0 1 0 0 0 0 0 0 1 0 1 0 0 1 0 0 0 1 1 1 1 1 1 0 0 0 0 1 1 1 0 0 1 1 0 1 0 0 0 1 0 0 1 1 0 1 0 0 1 1 1 1 0 1 0 1 0 1 1 1 0 0 0 0 1 1 1 1 0 0 1 0 1 0 1 0 0 0 1 1 1 0 1 1 1 0 1 1 1 0 0 1 0 1 1 0 1 0 0 1 1 1 0 0 0 1 0 0 0 1 0 0 1 1 0 0 1 0 0 1 1 0 1 0 1 1 0]";
	ss >> b;

	e = StrongEC(EllipticCurve(b));
	if (e.getOrder()
			!= ZZ(INIT_VAL,
					"15541351137805832567355695254588151253139252155490566670891172872358422298184284915720")) {
		cout << "Failed" << endl;
		return false;
	}
	cout << "OK" << endl;

	cout << "Checking for 283 bits. #2..." << std::flush;

	ss.clear();
	ss
			<< "[0 1 0 1 0 1 0 0 1 1 1 0 0 0 0 1 1 1 1 0 0 0 0 0 0 1 0 1 0 0 1 0 1 1 0 0 1 0 0 1 0 1 0 1 1 1 0 0 0 0 0 1 0 0 0 1 1 0 1 1 0 1 0 0 0 0 0 1 1 0 1 1 0 1 1 1 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1 1 1 1 0 0 0 1 1 1 1 0 1 0 0 0 1 0 0 0 0 1 1 0 1 1 1 0 0 0 1 0 0 0 0 1 1 1 1 1 1 1 0 0 1 0 1 0 1 1 1 0 0 1 0 1 1 0 1 1 1 0 1 0 0 1 0 1 0 1 1 0 0 0 1 0 0 0 1 0 1 1 0 0 0 0 0 1 0 1 0 1 1 0 1 1 1 0 1 0 1 0 0 1 1 1 0 0 1 1 1 1 0 0 0 1 1 0 0 0 1 1 0 1 1 1 1 0 1 0 1 0 1 1 1 1 0 1 0 0 0 1 1 1 0 1 1 1 0 1 1 0 1 0 1 0 1 1 1 1 1 1 1 0 1 1 0 0 0 0 0 0 1 1 0 0 0 1 1 1 0 1 1 0 0 1 1 1 1]";
	ss >> b;

	e = StrongEC(EllipticCurve(b));
	if (e.getOrder()
			!= ZZ(INIT_VAL,
					"15541351137805832567355695254588151253139260135319672559871374792149956092824258339180")) {
		cout << "Failed" << endl;
		return false;
	}
	cout << "OK" << endl;

	cout << endl;
	return true;
}

bool testFindStrongEC() {
	cout << "Testing Finding strong EC:" << endl;
	NTL::SetSeed(ZZ(INIT_VAL, clock()));

	cout << "Checking 127 bits..." << std::flush;
	NTL::GF2X P;
	NTL::BuildIrred(P, 127);
	GF2E::init(P);
	Constants::initConstants();

	StrongEC e;
	if (!e.isStrong()) {
		cout << "Error" << endl;
		return false;
	}
	cout << "OK" << endl;

	cout << "Checking 160 bits..." << std::flush;
	NTL::BuildIrred(P, 160);
	GF2E::init(P);
	Constants::initConstants();

	e = StrongEC();
	if (!e.isStrong()) {
		cout << "Error" << endl;
		return false;
	}
	cout << "OK" << endl;

	cout << endl;
	return true;
}

bool testECDH() {
	cout << "Testing ECDH:" << endl;
	std::stringstream ss;

	NTL::GF2X P;
	NTL::BuildIrred(P, 127);
	GF2E::init(P);
	Constants::initConstants();

	cout << "Building 127 bit curve" << endl;
	ss
			<< "[0 0 1 0 1 1 0 0 1 1 0 1 0 0 1 0 0 0 0 1 0 0 0 0 0 0 0 0 1 1 1 0 1 1 0 1 0 0 0 1 1 0 0 0 0 0 0 0 0 0 1 0 1 1 1 0 0 0 1 1 0 1 1 0 0 1 1 1 0 0 0 1 0 1 0 0 1 1 0 1 1 0 1 0 0 0 1 1 0 0 0 1 1 0 0 1 0 1 0 0 0 1 1 0 1 0 1 0 1 0 1 1 0 0 1 1 0 1 0 1 1 1 0 0 0 1]";
	ss >> P;
	StrongEC e(EllipticCurve(NTL::conv<GF2E>(P)));
	e.setGenerator();

	ss.clear();
	ss
			<< "[1 0 1 0 1 1 0 0 1 1 0 1 0 0 1 0 0 0 0 1 0 0 0 0 0 0 0 0 1 1 1 0 1 1 0 1 0 0 0 1 1 0 0 0 0 0 0 0 0 0 1 0 1 1 1 0 0 0 1 1 0 1 1 0 0 1 1 1 0 0 0 1 0 1 0 0 1 1 0 1 1 0 1 0 0 0 1 1 0 0 0 1 1 0 0 1 0 1 0 0 0 1 1 0 1 0 1 0 1 0 1 1 0 0 1 1 0 1 0 1 1 1 0 0 0 1]";
	ss >> P;
	StrongEC e2(EllipticCurve(NTL::conv<GF2E>(P)));
	e2.setGenerator();

	cout << "Checking 127 bits constructor..." << std::flush;

	ZZ k1(INIT_VAL, "123456789");
	ECDHClient c1(e, k1);
	if (c1.getPvKey() != k1
			|| c1.getPubKey() != e.pointMul(e.getGenerator(), k1)) {
		cout << "Error" << endl;
		return false;
	}
	cout << "OK" << endl;

	cout << "Checking 127 bits ecdh with known pks...";
	ZZ k2(INIT_VAL, "9876543210");
	ECDHClient c2(e, k2);
	ecdh(c1, c2);
	if (c1.getSharedSecret() != e.pointMul(e.getGenerator(), k1 * k2)
			|| c2.getSharedSecret() != e.pointMul(e.getGenerator(), k1 * k2)) {
		cout << "Error" << endl;
		return false;
	}
	cout << "OK" << endl;

	cout << "Checking 127 bits ecdh with random pks...";
	for (int i = 0; i < 10; i++) {
		c1 = ECDHClient(e);
		c2 = ECDHClient(e);
		ecdh(c1, c2);
		if (c1.getSharedSecret()
				!= e.pointMul(e.getGenerator(), c1.getPvKey() * c2.getPvKey())
				|| c2.getSharedSecret()
						!= e.pointMul(e.getGenerator(),
								c1.getPvKey() * c2.getPvKey())) {
			cout << "Error" << endl;
			return false;
		}
	}
	cout << "OK" << endl;

	cout << "Checking 127 bits with different curves...";

	c1 = ECDHClient(e);
	c2 = ECDHClient(e2);
	try {
		ecdh(c1, c2);
		cout << "Error" << endl;
		return false;
	} catch (CurvesDontMatch& e) {
		cout << "OK" << endl;
	}

	cout << "Checking 127 bits with different generators...";
	c1 = ECDHClient(e);
	c2 = ECDHClient(e);
	c1.setGenerator(e.getCurve().genPoint());
	c2.setGenerator(e.getCurve().genPoint());
	try {
		ecdh(c1, c2);
		cout << "Error" << endl;
		return false;
	} catch (CurvesDontMatch& e) {
		cout << "OK" << endl;
	}

	NTL::BuildIrred(P, 160);
	GF2E::init(P);
	Constants::initConstants();

	cout << "Building 160 bit curve" << endl;
	ss.clear();
	ss
			<< "[0 1 0 0 1 1 1 1 0 1 0 0 1 0 1 0 1 0 1 0 0 0 1 1 1 1 0 0 0 1 0 1 1 0 1 0 0 0 0 0 0 1 1 1 1 0 0 0 1 0 0 1 0 0 0 0 0 0 1 0 1 1 0 0 0 1 1 0 1 0 0 0 1 1 0 0 0 0 0 1 1 1 1 0 0 1 1 1 1 0 0 1 0 0 1 1 0 0 0 1 1 0 1 0 0 0 1 0 1 1 1 1 1 0 0 0 0 1 0 0 1 0 1 1 0 1 0 1 0 1 0 0 0 1 1 0 0 0 0 1 1 0 1 1 0 1 0 0 1 1 0 1 1 0 1]";
	ss >> P;
	e = StrongEC(EllipticCurve(NTL::conv<GF2E>(P)));
	e.setGenerator();

	cout << "Checking 160 bits constructor..." << std::flush;

	c1 = ECDHClient(e, k1);
	if (c1.getPvKey() != k1
			|| c1.getPubKey() != e.pointMul(e.getGenerator(), k1)) {
		cout << "Error" << endl;
		return false;
	}
	cout << "OK" << endl;

	cout << "Checking 160 bits ecdh with known pks...";
	c2 = ECDHClient(e, k2);
	ecdh(c1, c2);
	if (c1.getSharedSecret() != e.pointMul(e.getGenerator(), k1 * k2)
			|| c2.getSharedSecret() != e.pointMul(e.getGenerator(), k1 * k2)) {
		cout << "Error" << endl;
		return false;
	}
	cout << "OK" << endl;

	cout << "Checking 160 bits ecdh with random pks...";
	for (int i = 0; i < 10; i++) {
		c1 = ECDHClient(e);
		c2 = ECDHClient(e);
		ecdh(c1, c2);
		if (c1.getSharedSecret()
				!= e.pointMul(e.getGenerator(), c1.getPvKey() * c2.getPvKey())
				|| c2.getSharedSecret()
						!= e.pointMul(e.getGenerator(),
								c1.getPvKey() * c2.getPvKey())) {
			cout << "Error" << endl;
			return false;
		}
	}
	cout << "OK" << endl;

	cout << endl;
	return true;
}

int main(int argc, const char* argv[]) {
	srand(time(NULL));
	bool all_ok = true;

	if (!testEqualityINF()) {
		all_ok = false;
		printf("We failed in equalityINF!\n\n");
	}

	if (!testEquality()) {
		all_ok = false;
		printf("We failed in equality!\n\n");
	}

	if (!testIsOn()) {
		all_ok = false;
		printf("We failed in isOn!\n\n");
	}

	if (!testToAffine()) {
		all_ok = false;
		printf("We failed in toAffine\n\n");
	}

	if (!testGenPoint()) {
		all_ok = false;
		printf("we failed in genPoint\n\n");
	}

	if (!testDouble()) {
		all_ok = false;
		printf("we failed in pointDoubling\n\n");
	}

	if (!testAdd()) {
		all_ok = false;
		printf("we failed in pointAdd\n\n");
	}

	if (!testMul()) {
		all_ok = false;
		printf("we failed in pointMul\n\n");
	}

	if (!testMulRand()) {
		all_ok = false;
		printf("we failed in pointMulRand\n\n");
	}

	if (!testSatoh()) {
		all_ok = false;
		printf("we failed in testSatoh\n\n");
	}

	if (!testFindStrongEC()) {
		all_ok = false;
		printf("we failed in testFindStrongEC\n\n");
	}

	if (!testECDH()) {
		all_ok = false;
		printf("we failed in testECDH\n\n");
	}

//	testMulTime();

	cout << "********************" << endl;
	cout << "Done testing" << endl;
	cout << (all_ok ? "All passed" : "Something went wrong") << endl;
	cout << "********************" << endl;

	return 0;
}

bool testMulTime() {
	cout << "Testing Mul Time:" << endl;
	std::stringstream ss;

	NTL::GF2X P;
	ss.clear();
	ss
			<< "0x800000000000000000000000000000000000000000000000000000000000000000010a1";
	ss >> P;
	P.normalize();
	NTL::reverse(P, P);
	GF2E::init(P);
	Constants::initConstants();

	GF2X _x, _y;
	ss.clear();
	ss
			<< "0x6893bbe18b4775bfa5e33e7543eef52fb57950d5df9e891fd5c58203fdd795991d952c8";
	ss >> _x;
	NTL::reverse(_x, _x);
	ss.clear();
	ss
			<< "0x1e550c0d00d86d37b3cc2bfd3d059b57ffc90e012c6cc20d4a81bcb8b30d94a8719e279";
	ss >> _y;
	NTL::reverse(_y, _y);
	GF2E x, y;
	x = NTL::conv<GF2E>(_x);
	y = NTL::conv<GF2E>(_y);
	ProjectivePoint p = ProjectivePoint(x, y, Constants::one);

	GF2E t1, t2, t3;
	sqr(t1, y);
	mul(t2, x, y);
	sqr(t3, x);
	mul(t3, t3, x);
	add(t1, t1, t2);
	add(t1, t1, t3);

	EllipticCurve e = EllipticCurve(t1);

	if (!e.isOn(p)) {
		cout << "something is wrong!" << endl;
	}

	double seconds = 0, time1, time2;
	int times = 0;
	ZZ k;
	set(k);
	for (int j = 0; j < 20; j++) {
		for (int i = 220; i < 270; i++) {
			ZZ n = k << i; //2^i
			int toAdd = rand();
			n = n + toAdd;
			time1 = NTL::GetTime();
			e.pointMulMontogmery(p, n);
			time2 = NTL::GetTime();
			seconds += (time2 - time1);
			times++;
		}
	}
	cout << "Montgamery:" << endl;
	cout << "Total time: " << seconds << endl;
	cout << "Avarage: " << seconds / times << endl;

	cout << endl;

	return true;
}

void statSatoh() {
	cout << "Testing Satoh counting points:" << endl;
	NTL::GF2X P;
	std::stringstream ss;
	NTL::SetSeed(ZZ(INIT_VAL, clock()));

	{
		NTL::SetCoeff(P, 283);
		NTL::SetCoeff(P, 12);
		NTL::SetCoeff(P, 7);
		NTL::SetCoeff(P, 5);
		NTL::SetCoeff(P, 0);
		P.normalize();

		GF2E::init(P);
	}

//	{
//		NTL::SetCoeff(P, 11);
//		NTL::SetCoeff(P, 2);
//		NTL::SetCoeff(P, 0);
//		P.normalize();
//		std::cout << "Irred? " << NTL::IterIrredTest(P) << std::endl;
//		GF2E::init(P);
//	}

	GF2E b;

	{
		ss.clear();
		ss
				<< "[0 1 0 0 1 1 0 1 0 1 0 1 0 0 0 0 1 1 1 0 1 1 0 0 1 1 1 0 1 0 1 0 0 1 0 1 0 0 0 1 0 0 0 0 1 1 1 0 1 1 0 1 1 1 1 0 1 0 0 0 1 0 1 1 0 0 1 1 0 0 0 0 1 0 0 1 1 0 1 0 1 1 0 0 1 1 1 0 0 1 0 0 0 1 1 0 1 1 1 1 1 0 1 0 1 1 0 1 0 1 0 1 1 0 1 1 0 1 1 0 1 1 1 0 1 0 1 0 0 1 0 0 0 0 1 1 1 1 0 0 0 1 1 1 1 1 1 1 1 1 0 1 1 0 1 0 0 0 0 0 0 1 0 1 0 0 1 0 0 0 1 1 1 1 1 1 0 0 0 0 1 1 1 0 0 1 1 0 1 0 0 0 1 0 0 1 1 0 1 0 0 1 1 1 1 0 1 0 1 0 1 1 1 0 0 0 0 1 1 1 1 0 0 1 0 1 0 1 0 0 0 1 1 1 0 1 1 1 0 1 1 1 0 0 1 0 1 1 0 1 0 0 1 1 1 0 0 0 1 0 0 0 1 0 0 1 1 0 0 1 0 0 1 1 0 1 0 1 1 0]";
		ss >> b;
	}
//	{
//		ss.clear();
//		ss
//				<< "[0 1 0 1 0 1 0 0 1 1 1 0 0 0 0 1 1 1 1 0 0 0 0 0 0 1 0 1 0 0 1 0 1 1 0 0 1 0 0 1 0 1 0 1 1 1 0 0 0 0 0 1 0 0 0 1 1 0 1 1 0 1 0 0 0 0 0 1 1 0 1 1 0 1 1 1 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1 1 1 1 0 0 0 1 1 1 1 0 1 0 0 0 1 0 0 0 0 1 1 0 1 1 1 0 0 0 1 0 0 0 0 1 1 1 1 1 1 1 0 0 1 0 1 0 1 1 1 0 0 1 0 1 1 0 1 1 1 0 1 0 0 1 0 1 0 1 1 0 0 0 1 0 0 0 1 0 1 1 0 0 0 0 0 1 0 1 0 1 1 0 1 1 1 0 1 0 1 0 0 1 1 1 0 0 1 1 1 1 0 0 0 1 1 0 0 0 1 1 0 1 1 1 1 0 1 0 1 0 1 1 1 1 0 1 0 0 0 1 1 1 0 1 1 1 0 1 1 0 1 0 1 0 1 1 1 1 1 1 1 0 1 1 0 0 0 0 0 0 1 1 0 0 0 1 1 1 0 1 1 0 0 1 1 1 1]";
//		ss >> b;
//	}
//	{
//		ss.clear();
//		ss
//				<< "[0 0 0 1 0 0 0 1 1 1 0 1 1 1 0 0 1 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 1 1 1 0 0 1 0 0 0 0 0 1 1 0 1 0 1 0 0 0 1 0 1 0 1 1 1 0 0 1 0 0 0 1 1 1 0 0 0 1 1 0 1 1 0 1 1 1 1 1 0 0 0 1 1 1 0 1 1 0 0 1 0 1 0 1 0 0 0 0 0 1 1 1 1 1 1 1 1 0 1 1 0 1 0 1 0 1 1 0 1 1 0 1 0 0 0 1 0 1 1 0 1 1 1 1 0 0 1 0 1 0 1 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 0 1 1 1 1 1 0 1 0 0 0 1 1 0 1 1 0 0 0 1 1 0 1 1 1 1 1 1 1 1 0 0 1 1 0 0 0 1 0 1 0 1 0 1 0 1 1 0 1 1 1 1 1 1 1 0 1 1 0 0 0 0 1 1 0 1 0 0 1 1 0 1 1 0 1 1 1 0 0 1 0 1 0 0 0 1 1 1 0 0 1 1 0 1 0 1 0 0 1 1 1 0 1 1 1 0 1]";
//		ss >> b;
//	}
//	{
//		ss.clear();
//		ss
//				<< "[0 1 0 1 0 1 0 0 0 1 1 0 0 1 1 0 0 1 0 1 0 0 1 1 0 1 1 1 1 0 1 1 1 0 0 1 0 1 1 0 0 1 0 1 0 1 1 1 1 1 1 1 0 0 1 0 0 1 1 1 0 1 0 1 0 0 0 1 0 0 0 0 0 1 1 1 1 0 1 1 0 0 1 1 0 1 0 0 0 1 1 1 0 1 1 1 0 1 1 0 0 0 0 0 0 1 0 0 1 0 1 1 1 1 0 1 1 1 1 1 1 0 1 1 0 1 0 0 1 1 1 1 0 1 0 1 1 0 0 1 0 0 0 1 0 1 0 0 1 0 1 0 0 1 1 1 1 0 0 1 0 0 0 1 1 1 1 1 0 0 0 1 1 0 1 0 1 0 1 1 0 1 0 0 1 1 1 0 1 1 0 0 1 0 1 0 1 1 1 0 1 1 0 1 0 1 1 1 0 0 0 1 0 1 1 1 1 1 1 1 1 1 1 0 1 0 0 0 0 0 1 1 0 0 0 1 0 1 0 0 0 1 1 0 0 1 1 1 1 1 1 0 1 0 1 0 0 0 1 0 1 1 0 0 0 0 1 1 0 0 0 1 1 1 1 1 1 0 1 1 1]";
//		ss >> b;
//	}

//	{
//		ss.clear();
//		ss << "[1 1 1 1 1 1 1 1 1 1 1]";
//		ss >> b;
//	}
	Constants::initConstants();

	EllipticCurve curve = EllipticCurve(b);

	double NTL_start, NTL_end;

	cout << "GF2E::modulus: = " << GF2E::modulus() << endl;
	cout << "EC: " << b << endl;

	NTL_start = NTL::GetTime();
	StrongEC e = StrongEC(curve);
	NTL_end = NTL::GetTime();

	cout << "*************" << endl << " DONE " << endl << "*************"
			<< endl;
	cout << "Total time: " << NTL_end - NTL_start << " seconds" << endl;

	cout << "Order (p*n): " << e.getOrder() << endl;
	cout << "p: " << e.getP() << endl;
	cout << "n: " << e.getN() << endl;
	cout << "isStrong: " << (e.isStrong() ? "Yes" : "No") << endl;

	return;

	if (!e.isStrong()) {
		return;
	}
	NTL_start = NTL::GetTime();
	ProjectivePoint g = e.getGenerator();
	NTL_end = NTL::GetTime();

	cout << "Finding generator time: " << NTL_end - NTL_start << endl;
	cout << "g =  " << g << endl;
	cout << "n*g = " << e.getCurve().pointMulMontogmery(g, e.getN()) << endl;
	cout << "p*g = " << e.getCurve().pointMulMontogmery(g, e.getP()) << endl;

	cout << "Generating random strong points" << endl;
	for (int i = 0; i < 10; i++) {
		cout << "iter #" << i << endl;
		g = e.genStrongPoint();
		cout << "g =  " << g << endl;
		cout << "n*g = " << e.getCurve().pointMulMontogmery(g, e.getN())
				<< endl;
		cout << "p*g = " << e.getCurve().pointMulMontogmery(g, e.getP())
				<< endl;
		cout << endl;
	}

	return;
}

