/*
 * Challenge.cpp
 *
 *  Created on: Oct 14, 2014
 *      Author: Itay
 */

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

#include "AffineException.hpp"
#include "ProjectivePoint.hpp"
#include "EllipticCurve.hpp"
#include "StrongEC.hpp"
#include "CryptoLib.hpp"

using NTL::GF2E;
using NTL::GF2X;
using NTL::ZZ;

using std::cout;
using std::endl;

/**
 * recieves a point on the curve (different from infinity) and returns
 * the curve that the point is on such that a=0.
 */
const EllipticCurve findParameterFromPoint(const ProjectivePoint& p);

const ZZ findLargestOddDivisor(const ZZ& order);

const EllipticCurve findParameterFromPoint(const ProjectivePoint& p) {
	GF2E b;
	b = sqr(p.getY()) + p.getX() * p.getY() + power(p.getX(), 3);
	return EllipticCurve(Constants::zero, b);
}

const ZZ findLargestOddDivisor(const ZZ& order) {
	ZZ divisor = order;
	while (NTL::divide(divisor, divisor, Constants::ZZ_2) == 1) {
	}
	return divisor;
}

int main(int argc, const char* argv[]) {

	std::stringstream ss;

	//Initializing the challenge modulus.
	NTL::GF2X P;
	NTL::SetCoeff(P, 283);
	NTL::SetCoeff(P, 12);
	NTL::SetCoeff(P, 7);
	NTL::SetCoeff(P, 5);
	NTL::SetCoeff(P, 0);
	std::cout << "P = " << P << std::endl;
	P.normalize();
	GF2E::init(P);

	ProjectivePoint::initInf();
	Constants::storeGF2EModVec(ceil(GF2E::degree() / 2.0) + 13);

	/**
	 * Challenge part 1:
	 * We use G to find the curve's b parameter, and then count the points on that curve.
	 * We assume that G's order is the large prime of the order p, so we will find it and check if its
	 * the order.
	 */

	//Creating G:
	GF2E g_x, g_y;
	ss
			<< "[0 0 0 1 0 0 1 1 0 1 0 0 1 0 1 0 1 0 0 1 1 0 1 1 1 0 0 0 1 0 0 1 1 0 0 1 1 0 1 0 1 0 0 1 1 1 1 0 1 0 1 1 1 0 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 1 0 0 0 0 0 1 1 0 1 0 0 0 1 1 1 0 1 0 1 0 1 1 1 1 1 1 1 0 0 0 1 0 0 1 0 0 0 1 0 1 1 1 1 0 0 1 1 1 1 1 1 0 1 1 1 0 1 0 1 0 1 1 0 0 0 0 1 0 1 0 1 0 0 1 1 1 1 0 1 0 1 0 1 1 0 1 1 1 1 1 0 1 0 0 1 0 1 0 1 1 1 1 0 1 1 1 0 1 1 1 1 1 0 0 0 0 1 0 1 0 1 0 1 1 1 0 0 1 1 1 1 1 0 0 1 1 0 0 0 1 1 1 1 0 1 0 0 1 0 1 1 1 1 1 1 1 0 1 1 0 1 0 1 1 1 0 1 1 1 0 0 0 1 0 1 1 0 1 0 0 0 1 1 0 0 0 0 1 1 1 1 1 0 1 1 1 0 1 1 1 0 0 1 0 0 1 0 0 0 1 0 1 1]";
	ss >> g_x;
	ss.clear();

	ss
			<< "[1 0 0 1 1 1 1 0 0 1 0 0 0 1 1 1 1 0 0 1 1 0 0 0 1 1 1 0 0 0 0 1 0 1 0 1 0 0 1 0 1 0 0 1 1 0 1 1 0 0 0 0 1 1 0 0 1 1 0 1 0 0 0 1 1 1 0 1 0 0 1 1 1 1 0 1 1 0 0 0 0 0 0 1 0 1 0 1 0 0 1 0 1 0 1 1 0 0 0 0 0 1 0 0 0 0 1 1 0 0 1 1 0 1 1 0 0 0 1 1 0 1 0 0 1 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 1 0 0 1 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 0 1 0 1 1 0 1 1 0 0 1 1 0 1 0 0 0 0 0 1 0 1 1 1 1 0 0 1 0 1 1 1 1 1 1 1 1 0 1 0 1 0 0 0 0 1 1 0 0 1 1 1 1 0 0 1 1 0 1 1 1 1 0 1 1 0 0 1 0 1 1 0 1 1 0 0 0 0 1 1 0 1 1 0 0 0 0 0 0 0 0 1 0 1 1 0 0 0 0 0 0 1 1 0 0 0 0 1 0 1 0 1 0 1 0 0 1 1 1 1]";
	ss >> g_y;
	ss.clear();
	ProjectivePoint G = ProjectivePoint(g_x, g_y);
	cout << "G = " << G << endl;

	//Finding b from G:
	EllipticCurve curve = findParameterFromPoint(G);
	cout << "curve = " << curve.getB() << endl;

	//Counting the number of points on the curve:
	StrongEC e = StrongEC(curve);
	cout << "Order: " << e.getOrder() << endl;
	cout << "P: " << e.getP() << endl;
	cout << "n: " << e.getN() << endl;
	cout << "isStrong: " << (e.isStrong() ? "Yes" : "No") << endl;

	//Checking whether G's order is p:
	ProjectivePoint test_point = curve.pointMulMontogmery(G, e.getP());
	cout << "p*G = " << test_point << endl;
	if (test_point.isINF()) {
		cout << "G's order is " << e.getP() << endl;
	}

	/**
	 * Challenge part 2:
	 * We use alice's and bob's points to find alice_b, and bob_b - the parameters of
	 * the curves they generated.
	 * Then, we count the number of points on the curves and divide them by 2 as much as we can,
	 * the numbers remaining should be the private keys.
	 */

	//Computing alice's point:
	GF2E alice_x, alice_y;

	ss
			<< "[0 1 0 0 1 0 1 0 1 0 1 0 0 1 1 0 1 1 0 0 0 0 1 1 0 0 0 0 0 0 1 0 0 1 0 0 0 1 0 1 0 1 0 1 0 1 0 1 0 0 0 1 1 0 0 1 1 0 1 1 0 0 1 0 1 1 0 0 0 0 1 1 0 1 1 0 0 0 1 0 1 0 1 1 0 0 1 1 0 1 1 1 1 1 0 1 0 1 1 1 0 0 1 0 0 1 1 1 0 0 1 0 0 1 1 0 0 1 0 1 0 0 1 1 0 0 0 1 0 1 0 1 0 0 0 1 0 1 1 0 0 1 0 0 1 1 1 0 1 1 1 1 0 1 1 0 0 0 1 0 0 1 1 0 0 1 1 0 0 1 0 0 0 1 1 0 1 1 0 0 1 1 1 0 1 0 0 0 0 1 1 1 0 0 0 0 0 1 0 1 0 0 0 0 1 1 1 0 1 0 1 0 1 0 0 0 1 1 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 1 1 1 1 1 1 0 1 0 1 1 0 1 1 0 0 0 1 0 0 0 0 0 1 0 1 0 0 0 0 1 0 1 1 1 1 1 0 1 1 1 0 1 0 0 1 1]";
	ss >> alice_x;
	ss.clear();

	ss
			<< "[0 1 0 1 1 0 0 0 1 1 0 1 0 0 0 1 0 0 0 1 0 1 1 0 1 0 0 1 0 1 1 0 0 1 0 1 1 1 1 1 0 0 1 1 1 0 0 0 1 1 0 1 1 0 1 0 0 0 1 0 1 0 1 0 0 1 0 0 0 1 1 0 1 0 0 1 1 0 0 0 1 1 1 0 0 0 1 1 0 0 1 1 0 0 1 1 0 1 0 0 0 1 0 0 0 1 0 1 1 1 0 0 1 1 0 1 1 0 1 0 1 1 1 0 1 1 1 1 0 0 0 1 1 0 0 0 0 1 1 0 1 1 0 1 1 1 1 0 0 0 1 1 0 1 1 1 0 0 1 0 1 0 1 0 0 0 1 1 1 1 0 1 1 1 0 1 0 1 0 0 1 0 0 0 1 0 0 0 1 1 0 1 1 1 0 0 1 0 0 0 0 1 0 1 1 0 0 0 1 1 1 0 1 0 1 1 0 0 1 0 0 1 0 1 0 0 0 1 1 0 0 0 0 1 0 1 1 0 1 0 0 1 1 1 0 0 0 0 1 1 0 0 0 1 0 0 1 0 0 0 0 1 1 1 1 0 1 1 0 1 0 1 0 0 1 0 1 0 0 0 1 1 1]";
	ss >> alice_y;
	ss.clear();
	ProjectivePoint P_alice = ProjectivePoint(alice_x, alice_y);
	cout << "Alice's point = " << P_alice << endl;

	//Computing bob's point:
	GF2E bob_x, bob_y;
	ss
			<< "[0 0 0 0 0 1 1 0 1 1 1 0 0 1 0 1 1 0 0 1 1 0 0 1 0 0 0 1 1 1 1 0 0 0 1 1 0 1 1 0 0 1 1 1 0 1 1 0 0 0 0 0 0 1 1 0 1 0 0 1 1 1 0 1 1 1 1 1 0 0 1 0 1 0 1 0 1 1 1 0 0 0 1 0 0 1 1 1 0 0 0 1 0 0 0 0 1 1 1 0 1 1 0 1 1 0 1 1 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 1 1 0 0 1 1 1 0 1 1 1 0 0 0 0 1 0 0 1 0 1 0 0 0 1 0 0 1 0 0 1 0 1 0 0 1 1 1 0 0 1 1 1 1 1 1 0 1 1 0 1 0 0 1 1 1 0 0 1 0 1 0 1 1 1 1 1 0 0 0 1 0 1 1 1 0 1 1 0 0 0 1 0 0 1 1 1 1 0 0 0 0 1 0 0 1 1 1 1 0 0 1 0 0 1 1 1 0 0 1 1 0 0 1 0 1 0 1 0 1 1 0 0 0 0 0 0 1 1 0 1 1 1 1 0 0 1 0 1 0 0 1 0 1 1 0 1 1 1 1 1 0 1 1 0 0 1 1]";
	ss >> bob_x;
	ss.clear();

	ss
			<< "[0 1 1 1 1 0 1 0 0 1 0 0 1 0 1 1 1 1 0 1 1 0 1 1 1 0 1 0 0 1 0 0 0 0 0 0 1 0 1 1 0 1 1 1 0 0 0 0 0 0 0 0 1 0 0 1 0 1 0 0 0 1 0 0 0 1 1 0 0 1 1 0 1 0 1 0 1 0 1 1 1 0 0 1 0 0 1 1 1 1 0 1 1 0 1 1 1 0 1 0 1 0 0 0 1 1 0 1 0 0 0 1 0 1 1 1 0 1 0 1 0 0 0 0 1 0 1 0 1 1 0 1 1 1 0 1 1 0 0 0 1 0 1 0 1 0 0 0 1 0 1 1 0 0 1 1 0 1 1 0 1 0 0 1 0 1 0 1 0 1 1 0 0 0 1 0 0 1 1 1 1 1 0 0 1 1 1 0 1 1 1 1 1 0 1 0 0 0 0 1 0 1 1 0 1 0 1 1 0 0 1 0 1 0 0 1 1 1 1 0 1 1 1 0 0 0 1 0 0 1 1 0 1 0 1 1 0 0 0 0 0 0 0 1 0 1 0 1 1 1 0 1 1 0 0 1 0 0 0 1 1 1 1 1 0 1 1 1 0 1 1 1 0 0 1 0 0 0 0 1 1 1 1]";
	ss >> bob_y;
	ss.clear();
	ProjectivePoint P_bob = ProjectivePoint(bob_x, bob_y);
	cout << "Bob's point = " << P_bob << endl;

	//Finding alice's and bob's curves:
	EllipticCurve alice_curve = findParameterFromPoint(P_alice);
	cout << "alice's curve = " << alice_curve.getB() << endl;

	EllipticCurve bob_curve = findParameterFromPoint(P_bob);
	cout << "bob's curve = " << bob_curve.getB() << endl;

	//Counting the number of points on alice's and bob's curves:
	ZZ alice_order = StrongEC(alice_curve).getOrder();
	cout << "alice's order = " << alice_order << endl;

	ZZ bob_order = StrongEC(bob_curve).getOrder();
	cout << "bob's order = " << bob_order << endl;

	//Calculating the largest odd divisors - aka the private keys:
	ZZ alice_key = findLargestOddDivisor(alice_order);
	cout << "alice's private key = " << alice_key << endl;

	ZZ bob_key = findLargestOddDivisor(bob_order);
	cout << "bob's private key = " << bob_key << endl;

	//Doing ECDH between alice and bob:
	e.setGenerator(G);
	ECDHClient alice = ECDHClient(e,alice_key);
	ECDHClient bob = ECDHClient(e,bob_key);
	ecdh(alice,bob);
	cout << "alice's and bob's shared secret is: " << alice.getSharedSecret() << endl;

	return 0;
}
