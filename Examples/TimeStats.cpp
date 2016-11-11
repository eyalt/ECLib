/*
 * Main.cpp
 *
 *  Created on: 18 баеч 2014
 *      Author: Eyal
 */
#include "StrongEC.hpp"
#include <NTL/GF2XFactoring.h>
#include <assert.h>
#include <time.h>

using std::cout;
using std::endl;

NTL::GF2X getIrred(long bits) {
	NTL::GF2X P;
	switch (bits) {
	case 113: {
		NTL::SetCoeff(P, 113);
		NTL::SetCoeff(P, 15);
		NTL::SetCoeff(P, 0);
		break;
	}
	case 163: {
		NTL::SetCoeff(P, 163);
		NTL::SetCoeff(P, 7);
		NTL::SetCoeff(P, 6);
		NTL::SetCoeff(P, 3);
		NTL::SetCoeff(P, 0);
		break;
	}
	case 233: {
		NTL::SetCoeff(P, 233);
		NTL::SetCoeff(P, 74);
		NTL::SetCoeff(P, 0);
		break;
	}
	case 283: {
		NTL::SetCoeff(P, 283);
		NTL::SetCoeff(P, 12);
		NTL::SetCoeff(P, 7);
		NTL::SetCoeff(P, 5);
		NTL::SetCoeff(P, 0);
		break;
	}
	case 409: {
		NTL::SetCoeff(P, 409);
		NTL::SetCoeff(P, 87);
		NTL::SetCoeff(P, 0);
		break;
	}
	default:
		assert(false);
	}
	assert(NTL::deg(P) == bits);
	assert(NTL::IterIrredTest(P));
	return P;
}

void testArithTime(long bits);
void testSatohTime(long bits);
void testFindStrong(long bits);

int main() {
	NTL::SetSeed(NTL::ZZ(NTL::INIT_VAL, clock()));
	long bit_arr[5] = { 113, 163, 233, 283, 409 };
	for (int j = 0; j < 3; j++) {
		cout << "******* Iter #" << j << " *******" << endl;
		for (int i = 0; i < 1; i++) {
			testFindStrong(bit_arr[i]);
		}
		cout << endl << endl;
	}
//	testArithTime(113);
//	testArithTime(163);
//	testArithTime(233);
//	testArithTime(283);
//	testArithTime(409);
	return 0;
}

void testFindStrong(long bits) {
	cout << "Testing finding strong on " << bits << " bits" << endl;
	NTL::GF2X P = getIrred(bits);
	GF2E::init(P);
	Constants::initConstants();

	double satoh_time = 0;
	long satoh_count = 0;
	long satoh_tries = 0;
	double gen_time = 0;
	long gen_count = 0;

	for (int i = 0; i < 10; i++) {
		satoh_count++;
		StrongEC sec;
		satoh_time += sec.total_time;
		satoh_tries += sec.num_tries;
		cout << "Satoh iter #" << i << " time = " << sec.total_time << " # = "
				<< sec.num_tries << endl;
		for (int j = 0; j < 1000; j++) {
			gen_count++;
			gen_time -= NTL::GetTime();
			sec.setGenerator();
			gen_time += NTL::GetTime();
		}
	}

	cout << "Satoh total time: " << satoh_time << endl;
	cout << "Satoh avg. time: " << satoh_time / satoh_count << endl;
	cout << "Satoh total tries: " << satoh_tries << endl;
	cout << "Satoh avg. tries: " << (double) satoh_tries / satoh_count << endl;
	cout << "Gen total time: " << gen_time << endl;
	cout << "Gen avg. time: " << (gen_time / gen_count) * pow(10, 6) << "us"
			<< endl;
	cout << endl;
}

void testSatohTime(long bits) {
	cout << "Testing satoh on " << bits << " bits" << endl;
	NTL::GF2X P = getIrred(bits);
	GF2E::init(P);
	Constants::initConstants();

	double satoh_time = 0;
	long satoh_count = 0;
	double gen_time = 0;
	long gen_count = 0;

	for (int i = 0; i < 3; i++) {
		satoh_count++;
		GF2E b;
		NTL::random(b);
		EllipticCurve e(b);
		StrongEC sec(e);
//		satoh_time += sec.time_satoh;
//		cout << "Satoh iter #" << i << " time = " << sec.time_satoh << endl;
		for (int j = 0; j < 1000; j++) {
			gen_count++;
			gen_time -= NTL::GetTime();
			e.genPoint();
			gen_time += NTL::GetTime();
		}
	}

	cout << "Satoh total time: " << satoh_time << endl;
	cout << "Satoh avg. time: " << satoh_time / satoh_count << endl;
	cout << "Gen total time: " << gen_time << endl;
	cout << "Gen avg. time: " << (gen_time / gen_count) * pow(10, 6) << "us"
			<< endl;
	cout << endl;
}

void testArithTime(long bits) {
	cout << "Testing arith on " << bits << " bits" << endl;
	NTL::GF2X P = getIrred(bits);
	GF2E::init(P);
	Constants::initConstants();

	GF2E x, y;
	double add_time = 0;
	double double_time = 0;
	double mul_time = 0;
	long count = 0;
	cout.precision(16);
	NTL::ZZ k;
	NTL::ZZ max_k = NTL::power(NTL::ZZ(NTL::INIT_VAL, 2), (int) (bits / 2));
	for (int i = 0; i < pow(10, 5); i++) {
		GF2E b;
		NTL::random(b);
		EllipticCurve e(b);
		NTL::RandomBnd(k, max_k);
		ProjectivePoint p1 = e.genPoint();
		ProjectivePoint p2;
		do {
			p2 = e.genPoint();
		} while (p2 == p1);

		add_time -= NTL::GetTime();
		e.pointAdd(p1, p2);
		add_time += NTL::GetTime();

		double_time -= NTL::GetTime();
		e.pointDoubling(p1);
		double_time += NTL::GetTime();

		mul_time -= NTL::GetTime();
		e.pointMulMontogmery(p1, k);
		mul_time += NTL::GetTime();

		count++;
	}

	cout << "Addition total time: " << add_time << endl;
	cout << "Addition avg. time: " << add_time / count << endl;
	cout << "Doubling total time: " << double_time << endl;
	cout << "Doubling avg. time: " << double_time / count << endl;
	cout << "Montogmery total time: " << mul_time << endl;
	cout << "Montogmery avg. time: " << mul_time / count << endl;

	cout << endl;
}
