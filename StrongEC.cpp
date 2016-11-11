/*
 * StrongEC.cpp
 *
 *  Created on: Jul 28, 2014
 *      Author: Itay
 */

#include "StrongEC.hpp"
#include <NTL/ZZ_pE.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZVec.h>
#include <assert.h>
#include <sstream>
#include "Constants.hpp"
#include <NTL/ZZ_pX.h>
#include <stack>
#include "EC_Exceptions.hpp"

using NTL::GF2E;
using NTL::ZZ;
using NTL::GF2X;
using NTL::ZZ_pE;
using NTL::ZZ_pEBak;
using NTL::ZZ_p;
using NTL::power;
using NTL::conv;
using NTL::INIT_VAL;
using NTL::ZZX;
using NTL::ZZ_pX;
using NTL::MulMod;
using NTL::SqrMod;
using NTL::rem;

using std::cout;
using std::endl;

int isFirst = 0;
double LIFTED_INV_TIME[2] = { 0, 0 };
double LIFTED_INV_NUMBER[2] = { 0, 0 };
double INV_MOD_NUMBER[2] = { 0, 0 };
double time_newton[2] = { 0, 0 };
double time_ABC[2] = { 0, 0 };
double time_ND[2] = { 0, 0 };
double time_ND_core_calc[2] = { 0, 0 };
double time_liftPrev[2] = { 0, 0 };
double time_applyRem[2] = { 0, 0 };
double time_applyRem_remainder[2] = { 0, 0 };
double time_main_calcs = 0;
double time_liftedInv_core[2] = { 0, 0 };
double time_liftedInv_calc[2] = { 0, 0 };

double time_calc_newton_mul[2] = { 0, 0 };
double number_calc_newton[2] = { 0, 0 };

void reset_timers() {
	isFirst = 0;
	LIFTED_INV_TIME[0] = LIFTED_INV_TIME[1] = 0;
	LIFTED_INV_NUMBER[0] = LIFTED_INV_NUMBER[1] = 0;
	INV_MOD_NUMBER[0] = INV_MOD_NUMBER[1] = 0;
	time_newton[0] = time_newton[1] = 0;
	time_ABC[0] = time_ABC[1] = 0;
	time_ND[0] = time_ND[1] = 0;
	time_ND_core_calc[0] = time_ND_core_calc[1] = 0;
	time_liftPrev[0] = time_liftPrev[1] = 0;
	time_applyRem[0] = time_applyRem[1] = 0;
	time_applyRem_remainder[0] = time_applyRem_remainder[1] = 0;
	time_main_calcs = 0;
	time_liftedInv_core[0] = time_liftedInv_core[1] = 0;
	time_liftedInv_calc[0] = time_liftedInv_calc[1] = 0;

	time_calc_newton_mul[0] = time_calc_newton_mul[1] = 0;
	number_calc_newton[0] = number_calc_newton[1] = 0;
}

ZZX liftPrevJInvRec(const ZZX& A, const ZZX& B, const ZZX& C, const ZZX& jNext,
		const int N);
ZZX calcNewton(const ZZX& A, const ZZX& B, const ZZX& C, const ZZX& j,
		const int N);

ZZX liftPrevJInv(const ZZX& _jNext, const int N);
ZZX liftFirstJInv(const GF2E& j0, const int N);
ZZX liftedInverse(const ZZX& a, const int N);
ZZX liftedInverseSqrt(const ZZX& a, const ZZX& z0, const int N);

void applyModP(ZZX& a);
void applyRemF(ZZX& a, long N);
void applyRemF(ZZX& a, ZZ_pX& temp, long N);
ZZX getRemF(const ZZX& a, long N);

void applyRemFMod(ZZX& a, long N);

void StrongEC::calcStrength() {
	ZZ i = Constants::ZZ_2;
	p = order;
	NTL::set(n);
	double div_time = -NTL::GetTime();
	for (; i < MAX_ALLOWED_PRIMES; i++) {
		while ((NTL::divide(p, p, i) == 1)) {
			mul(n, n, i);
		}
	}
	div_time += NTL::GetTime();

#ifdef SHOW_MR_STATS
	std::cout << "calcStrength division time: " << div_time << std::endl;
#endif

	double mr_time = -NTL::GetTime();
	strong = NTL::ProbPrime(p, MILLER_RABIN_TRIALS);
	mr_time += NTL::GetTime();

#ifdef SHOW_MR_STATS
	std::cout << "calcStrength miller rabin time: " << mr_time << std::endl;
#endif
}

ProjectivePoint StrongEC::pointAdd(const ProjectivePoint& p1,
		const ProjectivePoint& p2) const {
	return curve.pointAdd(p1, p2);
}
ProjectivePoint StrongEC::pointDoubling(const ProjectivePoint& p) const {
	return curve.pointDoubling(p);
}
ProjectivePoint StrongEC::pointMul(const ProjectivePoint& p,
		const NTL::ZZ& k) const {
	return curve.pointMulMontogmery(p, k);
}

ProjectivePoint StrongEC::getGenerator() const {
	if (!generator.isAffine()) {
		throw GeneratorNotSetException();
	}
	return generator;
}

void StrongEC::setGenerator() { //will replace the current one!
	if (!isStrong()) {
		generator = curve.genPoint();
		return;
	}
	bool found = false;
	while (!found) {
		generator = curve.genPoint();
		generator = curve.pointMulMontogmery(generator, n);
		if (!generator.isINF()) {
			found = true;
		}
	}
}

bool StrongEC::hasGenerator() const {
	return generator.isAffine();
}

ProjectivePoint StrongEC::genStrongPoint() const {
	ZZ k;
	RandomBnd(k, p - 1);
	NTL::add(k, k, 1);
	return curve.pointMulMontogmery(getGenerator(), k);
}

void StrongEC::countPointsSatoh() {
	reset_timers();

	int N = ceil(double(GF2E::degree()) / 2.0) + 13;
	int M = N - 10;

	ZZ two_pow_N = Constants::getPowerOf2(N);
	ZZ two_pow_M = Constants::getPowerOf2(M);

	ZZX j, jPrime;

	double start = NTL::GetTime();

	j = liftFirstJInv(curve.jInvariant(), N);

	double timeFirst = NTL::GetTime() - start;
	double timeJMain = 0;

	ZZX cn, cd, z, t, toInv;
	NTL::set(cn);
	NTL::set(cd);

	ZZ_p ZZp_N_195120, ZZp_N_4095, ZZp_N_660960000, ZZp_N_563760, ZZp_N_512,
			ZZp_N_372735, ZZp_N_8981280000;

	conv(ZZp_N_195120, Constants::ZZ_195120);
	conv(ZZp_N_4095, Constants::ZZ_4095);
	conv(ZZp_N_660960000, Constants::ZZ_660960000);
	conv(ZZp_N_563760, Constants::ZZ_563760);
	conv(ZZp_N_512, Constants::ZZ_512);
	conv(ZZp_N_372735, Constants::ZZ_372735);
	conv(ZZp_N_8981280000, Constants::ZZ_8981280000);

	const ZZ_pXModulus& ZZpX_N_f = Constants::loadGF2EModVec(N);

	ZZ_p::init(two_pow_M);
	ZZ_p ZZp_M_12, ZZp_M_1728, ZZp_M_36, ZZp_M_504, ZZp_M_12096, ZZp_M_240;

	conv(ZZp_M_12, Constants::ZZ_12);
	conv(ZZp_M_1728, Constants::ZZ_1728);
	conv(ZZp_M_36, Constants::ZZ_36);
	conv(ZZp_M_504, Constants::ZZ_504);
	conv(ZZp_M_12096, Constants::ZZ_12096);
	conv(ZZp_M_240, Constants::ZZ_240);

	ZZ_pX cn_pX, cd_pX;
	NTL::set(cn_pX);
	NTL::set(cd_pX);

	const ZZ_pXModulus& ZZpX_M_f = Constants::loadGF2EModVec(M);

	isFirst = 1;

	start = NTL::GetTime();
	for (int i = 0; i <= GF2E::degree() - 1; i++) {
		ZZ_p::init(two_pow_N);

		timeJMain -= NTL::GetTime();
		jPrime = liftPrevJInv(j, N);
		timeJMain += NTL::GetTime();

		time_main_calcs -= NTL::GetTime();
		ZZ_pX z_pX, j_pX, jPrime_pX, t_pX;
		ZZ_pX j2_pX;
		conv(j_pX, j);
		conv(jPrime_pX, jPrime);

		SqrMod(j2_pX, j_pX, ZZpX_N_f);

		z_pX = -(j2_pX + ZZp_N_195120 * j_pX + ZZp_N_4095 * jPrime_pX
				+ ZZp_N_660960000);
		conv(z, z_pX);
		NTL::div(z, z, 4096);
		conv(z_pX, z);

		ZZ_pX toInv_pX = j2_pX
				+ MulMod(j_pX, (ZZp_N_563760 - ZZp_N_512 * jPrime_pX), ZZpX_N_f)
				+ ZZp_N_372735 * jPrime_pX + ZZp_N_8981280000;
		conv(toInv, toInv_pX);
		NTL::div(toInv, toInv, 512);
		conv(toInv_pX, toInv);

		LIFTED_INV_TIME[isFirst] -= NTL::GetTime();
		conv(toInv_pX, liftedInverse(toInv, N));
		LIFTED_INV_TIME[isFirst] += NTL::GetTime();
		LIFTED_INV_NUMBER[isFirst]++;

		MulMod(z_pX, z_pX, toInv_pX, ZZpX_N_f);
		conv(z, z_pX);

		ZZ_p::init(two_pow_M);
		conv(j_pX, j);
		conv(jPrime_pX, jPrime);
		conv(z_pX, z);

		t_pX = MulMod(ZZp_M_12 * SqrMod(z_pX, ZZpX_M_f) + z_pX,
				(jPrime_pX - ZZp_M_1728), ZZpX_M_f) - ZZp_M_36;

		MulMod(cn_pX, cn_pX,
				jPrime_pX
						- MulMod(ZZp_M_504 + ZZp_M_12096 * z_pX, t_pX,
								ZZpX_M_f), ZZpX_M_f);

		MulMod(cd_pX, cd_pX, ZZp_M_240 * t_pX + jPrime_pX, ZZpX_M_f);

		j = jPrime;
		time_main_calcs += NTL::GetTime();
	}
	conv(cn, cn_pX);
	conv(cd, cd_pX);

	LIFTED_INV_TIME[isFirst] -= NTL::GetTime();
	ZZX cd_inv = liftedInverse(cd, M);
	LIFTED_INV_TIME[isFirst] += NTL::GetTime();
	LIFTED_INV_NUMBER[isFirst]++;

	ZZX res = cn * cd_inv;
	applyRemFMod(res, M);

	res = res * liftedInverseSqrt(res, Constants::ZZX_1, M);
	applyRemFMod(res, M - 1);

	ZZ ZZ_res = NTL::coeff(res, 0);
	if (sqr(ZZ_res) > Constants::getPowerOf2(GF2E::degree() + 2)) {
		ZZ_res -= Constants::getPowerOf2(M - 1);
	}

#ifdef SHOW_SATOH_STATS
	std::cout << "Total time in liftedInv: " << LIFTED_INV_TIME[0] << " "
	<< LIFTED_INV_TIME[1] << std::endl;
	cout << "Total time in lifted core: " << time_liftedInv_core[0] << " "
	<< time_liftedInv_core[1] << endl;
	cout << "Total time in lifted calc: " << time_liftedInv_calc[0] << " "
	<< time_liftedInv_calc[1] << endl;
	cout << "Total time in main loop j Inv: " << timeJMain << endl;

	cout << "Total time in newton: " << time_newton[0] << " " << time_newton[1]
	<< endl;
	cout << "Total time in ND: " << time_ND[0] << " " << time_ND[1] << endl;
	cout << "Total time in ND core: " << time_ND_core_calc[0] << " "
	<< time_ND_core_calc[1] << endl;
	cout << "Total time in calc newton mul: " << time_calc_newton_mul[0] << " "
	<< time_calc_newton_mul[1] << endl;
	cout << "# of times in calc newton mul: " << number_calc_newton[0] << " "
	<< number_calc_newton[1] << endl;
	cout << "Total time in ABC: " << time_ABC[0] << " " << time_ABC[1] << endl;
	cout << "Total time in liftPrevJInvRec: " << time_liftPrev[0] << " "
	<< time_liftPrev[1] << endl;
	cout << "Total time in applyRem: " << time_applyRem[0] << " "
	<< time_applyRem[1] << endl;
	cout << "Total time in applyRem remainder: " << time_applyRem_remainder[0]
	<< " " << time_applyRem_remainder[1] << endl;
	cout << "Total time in main calculations: " << time_main_calcs << endl;

	std::cout << "Total time in main loops: " << timeFirst << " "
	<< NTL::GetTime() - start << std::endl;

	cout << "Final t: " << ZZ_res << endl;
	cout << "#E: " << order << endl;
#endif

	order = Constants::getPowerOf2(GF2E::degree()) + 1 - ZZ_res;
}

void calcND(const ZZX& j, const ZZX& A, const ZZX& B, const ZZX& C, ZZ_pX& n,
		ZZ_pX& d, long N) {
	time_ND[isFirst] -= NTL::GetTime();
	ZZ_pX pxJ1, pxJ2, pxJ3, pxA, pxB, pxC;
	const ZZ_pXModulus& f = Constants::loadGF2EModVec(N);

	conv(pxA, A);
	conv(pxB, B);
	conv(pxC, C);

	conv(pxJ1, j);
	NTL::rem(pxJ1, pxJ1, f);

	SqrMod(pxJ2, pxJ1, f);
//	mul(pxJ3, pxJ1, pxJ2);

	time_ND_core_calc[isFirst] -= NTL::GetTime();

	n = MulMod(pxJ2, 2 * pxJ1 + pxA, f) - pxC;
	d = 3 * pxJ2 + 2 * MulMod(pxA, pxJ1, f) + pxB;
//	rem(n, n, f);
//	rem(d, d, f);

//	n = pxJ3 + pxA * pxJ2 + pxB * pxJ1 + pxC;
//	rem(n, n, f);
//	d = 3 * pxJ2 + 2 * pxA * pxJ1 + pxB;
//	rem(d, d, f);

//	n = pxJ3 + MulMod(pxA, pxJ2, f) + MulMod(pxB, pxJ1, f) + pxC;
//	d = 3 * pxJ2 + 2 * MulMod(pxA, pxJ1, f) + pxB;

	time_ND_core_calc[isFirst] += NTL::GetTime();
	time_ND[isFirst] += NTL::GetTime();
}

ZZX calcNewton(const ZZX& A, const ZZX& B, const ZZX& C, const ZZX& j,
		const int N) {
	time_newton[isFirst] -= NTL::GetTime();
	ZZ_pX n, d;
	calcND(j, A, B, C, n, d, N);

	ZZX ZZX_d, ret;
	conv(ZZX_d, d);

	LIFTED_INV_TIME[isFirst] -= NTL::GetTime();
	time_calc_newton_mul[isFirst] -= NTL::GetTime();
	conv(d, liftedInverse(ZZX_d, N));
	time_calc_newton_mul[isFirst] += NTL::GetTime();
	LIFTED_INV_TIME[isFirst] += NTL::GetTime();
	LIFTED_INV_NUMBER[isFirst]++;

	ZZ_pX tmpRet;

	conv(tmpRet, j);

	conv(ret, MulMod(n, d, Constants::loadGF2EModVec(N)));
	number_calc_newton[isFirst]++;

	time_newton[isFirst] += NTL::GetTime();
	return ret;
}

void calcABC(const ZZX& jNext, ZZX& A, ZZX& B, ZZX& C, long N) {
	time_ABC[isFirst] -= NTL::GetTime();
	ZZ_pX pxJNext1, pxJNext2, pxJNext3, pxA, pxB, pxC;
	const ZZ_pXModulus& f = Constants::loadGF2EModVec(N);

	conv(pxJNext1, jNext);
	NTL::rem(pxJNext1, pxJNext1, f);
	SqrMod(pxJNext2, pxJNext1, f);
	MulMod(pxJNext3, pxJNext1, pxJNext2, f);

	ZZ_p a0, a1, a2;

	conv(a0, Constants::ZZ_162000);
	conv(a1, Constants::ZZ_1488);
	pxA = -pxJNext2 + a1 * pxJNext1 - a0;

	conv(a2, Constants::ZZ_1488);
	conv(a1, Constants::ZZ_40773375);
	conv(a0, Constants::ZZ_8748000000);
	pxB = a2 * pxJNext2 + a1 * pxJNext1 + a0;

	conv(a2, Constants::ZZ_162000);
	conv(a1, Constants::ZZ_8748000000);
	conv(a0, Constants::ZZ_157464000000000);
	pxC = pxJNext3 - a2 * pxJNext2 + a1 * pxJNext1 - a0;

	conv(A, pxA);
	conv(B, pxB);
	conv(C, pxC);

	time_ABC[isFirst] += NTL::GetTime();
}

/***
 * @param jNext - mod 2^(N-1)
 * @return - j (calced from A,B,C,jNext) mod 2^N
 */
ZZX liftPrevJInvRec(const ZZX& A, const ZZX& B, const ZZX& C, const ZZX& jNext,
		const int N) {

	ZZX ret;

	if (N == 1) {
		GF2X tempRet;
		conv(tempRet, jNext);
		SqrMod(tempRet, tempRet, GF2E::modulus());
		conv(ret, tempRet);

		return ret;
	}

	int NPrime = ceil(double(N) / 2);

	NTL::ZZ_pBak bak;
	bak.save();
	ZZ_p::init(Constants::getPowerOf2(NPrime));
	ZZX temp = liftPrevJInvRec(A, B, C, jNext, NPrime);
	bak.restore();

	temp = calcNewton(A, B, C, temp, N);
	return temp;
}

/***
 * @return - j mod 2^N
 */
ZZX liftPrevJInv(const ZZX& jNext, const int N) {
	ZZX A, B, C;

//	assert(ZZ_p::modulus() == power(two, N));

	calcABC(jNext, A, B, C, N);
	time_liftPrev[isFirst] -= NTL::GetTime();
	ZZX res = liftPrevJInvRec(A, B, C, jNext, N);
	time_liftPrev[isFirst] += NTL::GetTime();
	return res;
}

ZZX liftFirstJInv(const GF2E& j0, const int N) {
	std::stringstream ss;

	GF2E tempJ0 = j0;
	long deg = GF2E::degree();

	for (int i = 0; i < deg - N + 1; i++) {
		NTL::sqr(tempJ0, tempJ0);
	}

	ZZX j;
	conv(j, conv<GF2X>(j0));

	ZZ_pX tempJ, f;

	for (int i = 2; i <= N; i++) {
		ZZ_p::init(Constants::getPowerOf2(i));
		j = liftPrevJInv(j, i);
		applyRemF(j, N);
	}

	return j;
}

/***
 * @return - a^(-1) mod 2^N
 */
ZZX liftedInverse(const ZZX& a, const int N) {
	ZZX ret, z;
	const ZZ_pXModulus& f = Constants::loadGF2EModVec(N);
	if (N == 1) {
		GF2X tempRet;

		conv(tempRet, a);

		time_liftedInv_core[isFirst] -= NTL::GetTime();
		NTL::rem(tempRet, tempRet, GF2E::modulus());
		NTL::InvMod(tempRet, tempRet, GF2E::modulus());
		time_liftedInv_core[isFirst] += NTL::GetTime();

		INV_MOD_NUMBER[isFirst]++;

		conv(ret, tempRet);
		return ret;
	}

	int n = ceil(double(N) / 2);

	NTL::ZZ_pBak bak;
	bak.save();
	ZZ_p::init(Constants::getPowerOf2(n));
	z = liftedInverse(a, n);
	bak.restore();

	ZZ_pX pxZ, pxA;
	conv(pxZ, z);
	conv(pxA, a);

//ret = z + z (1 - a * z)
	ZZ_pX pxRet;

	{
		time_liftedInv_calc[isFirst] -= NTL::GetTime();
		SqrMod(pxRet, pxZ, f);
		MulMod(pxRet, pxRet, pxA, f);
		sub(pxRet, 2 * pxZ, pxRet);
		time_liftedInv_calc[isFirst] += NTL::GetTime();
	}

	conv(ret, pxRet);
	return ret;

}

/**
 * @return sqrt(a)^(-1) mod 2^N
 */
ZZX liftedInverseSqrt(const ZZX& a, const ZZX& z0, const int N) {
	ZZX ret, z;

	if (N <= 2) {
		ret = z0;
		applyRemFMod(ret, N);
		return ret;
	}

	int n = ceil(double(N + 1) / 2);

	NTL::ZZ_pBak bak;
	bak.save();
	ZZ_p::init(Constants::getPowerOf2(n));
	z = liftedInverseSqrt(a, z0, n);
	bak.restore();

	ret = z + z * (1 - a * z * z) / 2;
	applyRemFMod(ret, N);

	return ret;
}

void applyModP(ZZX& a) {
	static ZZ_pX temp;
	conv(temp, a);
	conv(a, temp);
}

void applyRemF(ZZX& a, long N) {
	ZZ_pX temp;
	applyRemF(a, temp, N);
}

void applyRemF(ZZX& a, ZZ_pX& temp, long N) {
	time_applyRem[isFirst] -= NTL::GetTime();
	conv(temp, a);

	time_applyRem_remainder[isFirst] -= NTL::GetTime();
	NTL::rem(temp, temp, Constants::loadGF2EModVec(N));
	time_applyRem_remainder[isFirst] += NTL::GetTime();

	conv(a, temp);
	time_applyRem[isFirst] += NTL::GetTime();

}

ZZX getRemF(const ZZX& a, long N) {
	ZZX temp = a;
	applyRemF(temp, N);
	return temp;
}

void applyRemFMod(ZZX& a, long N) {
	NTL::ZZ_pBak bak;
	bak.save();
	ZZ_p::init(Constants::getPowerOf2(N));
	ZZ_pX temp;
	applyRemF(a, temp, N);
	bak.restore();
}

StrongEC::StrongEC() :
		curve(GF2E::zero()) {
	GF2E b;
	int i = 0;
	while (true) {
		i++;
//		cout << ">> Searching strong EC. iter #" << i << endl;
		NTL::random(b);
		while (NTL::IsZero(b)) {
			NTL::random(b);
		}

		curve = EllipticCurve(b);
		countPointsSatoh();
		calcStrength();
		if (isStrong()) {
			break;
		}
//		cout << endl;
	}
	setGenerator();

//	cout << endl << "took a total of " << i << " iterations" << endl;

}

const EllipticCurve& StrongEC::getCurve() const {
	return curve;
}

const ZZ& StrongEC::getN() const {
	return n;
}

const ZZ& StrongEC::getP() const {
	return p;
}

const ZZ& StrongEC::getOrder() const {
	return order;
}

bool StrongEC::isStrong() const {
	return strong;
}

void StrongEC::setGenerator(const ProjectivePoint& g) {
	if (g.isINF() || !curve.isOn(g)) {
		throw SetGeneratorException();
	}
	if (!g.isAffine()) {
		throw AffineException();
	}
	generator = g;
}

bool StrongEC::isStrongPoint(const ProjectivePoint& g) const {
	return (!g.isINF()) && curve.pointMulMontogmery(g, p).isINF();
}

bool operator==(const StrongEC& e1, const StrongEC& e2) {
	return e1.getCurve() == e2.getCurve()
			&& e1.getGenerator() == e2.getGenerator();
}
bool operator!=(const StrongEC& e1, const StrongEC& e2) {
	return !(e1 == e2);
}

bool StrongEC::isOn(const ProjectivePoint& p) const {
	return curve.isOn(p);
}

ProjectivePoint StrongEC::genPoint() const {
	return curve.genPoint();
}

GF2E StrongEC::jInvariant() const {
	return curve.jInvariant();
}
