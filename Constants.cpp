/*
 * Constants.cpp
 *
 *  Created on: Jul 28, 2014
 *      Author: Itay
 */

#include "Constants.hpp"
#include <NTL/ZZ_p.h>
#include "ProjectivePoint.hpp"

NTL::GF2E_NoAlloc_type a;
NTL::GF2E Constants::zero = NTL::GF2E(a);
NTL::GF2E Constants::one = NTL::GF2E(a);
NTL::ZZX Constants::GF2E_modulus = NTL::ZZX(NTL::INIT_SIZE, 1);

vector<NTL::ZZ_pXModulus> Constants::vec_GF2E_modulus(0);
vector<NTL::ZZ> Constants::powersOf2(0);

ZZX Constants::ZZX_1;

ZZ Constants::ZZ_2;
ZZ Constants::ZZ_12;
ZZ Constants::ZZ_36;
ZZ Constants::ZZ_240;
ZZ Constants::ZZ_504;
ZZ Constants::ZZ_512;
ZZ Constants::ZZ_1488;
ZZ Constants::ZZ_1728;
ZZ Constants::ZZ_4095;
ZZ Constants::ZZ_12096;
ZZ Constants::ZZ_162000;
ZZ Constants::ZZ_195120;
ZZ Constants::ZZ_372735;
ZZ Constants::ZZ_563760;
ZZ Constants::ZZ_40773375;
ZZ Constants::ZZ_660960000;
ZZ Constants::ZZ_8748000000;
ZZ Constants::ZZ_8981280000;
ZZ Constants::ZZ_157464000000000;

void Constants::initConstants() {
	zero = NTL::GF2E::zero();
	one = zero + 1;

	long N = 300 > (GF2E::degree() + 2) ? 300 : (GF2E::degree() + 2);

	storeGF2EMod();
	storeGF2EModVec(N);

	ZZX_1 = ZZX::zero() + 1;

	ZZ_2 = ZZ(INIT_VAL, "2");
	ZZ_12 = ZZ(INIT_VAL, "12");
	ZZ_36 = ZZ(INIT_VAL, "36");
	ZZ_240 = ZZ(INIT_VAL, "240");
	ZZ_504 = ZZ(INIT_VAL, "504");
	ZZ_512 = ZZ(INIT_VAL, "512");
	ZZ_1488 = ZZ(INIT_VAL, "1488");
	ZZ_1728 = ZZ(INIT_VAL, "1728");
	ZZ_4095 = ZZ(INIT_VAL, "4095");
	ZZ_12096 = ZZ(INIT_VAL, "12096");
	ZZ_162000 = ZZ(INIT_VAL, "162000");
	ZZ_195120 = ZZ(INIT_VAL, "195120");
	ZZ_372735 = ZZ(INIT_VAL, "372735");
	ZZ_563760 = ZZ(INIT_VAL, "563760");
	ZZ_40773375 = ZZ(INIT_VAL, "40773375");
	ZZ_660960000 = ZZ(INIT_VAL, "660960000");
	ZZ_8748000000 = ZZ(INIT_VAL, "8748000000");
	ZZ_8981280000 = ZZ(INIT_VAL, "8981280000");
	ZZ_157464000000000 = ZZ(INIT_VAL, "157464000000000");

	storePowersOf2(N);

	ProjectivePoint::initInf();
}

void Constants::storeGF2EMod() {
	conv(GF2E_modulus, NTL::GF2E::modulus());
}

/*stores whatever we want to the GF2E_modulus... use only for testing!*/
void Constants::storeTestMod(const NTL::ZZ_pX& mod) {
	conv(GF2E_modulus, mod);
}

const NTL::ZZX& Constants::loadGF2EMod() {
	return GF2E_modulus;
}

void Constants::storeGF2EModVec(long maxDeg) {
	vec_GF2E_modulus = vector<ZZ_pXModulus>(0);
	vec_GF2E_modulus.reserve(maxDeg);
	NTL::ZZ_pBak bak;
	bak.save();
	ZZ mod = ZZ(INIT_VAL, "2");
	ZZX f;
	NTL::ZZ_p::init(mod);
	NTL::ZZ_pX ZZ_pX_f;
	conv(f, NTL::GF2E::modulus());
	for (int i = 1; i <= maxDeg; i++) {
		NTL::ZZ_p::init(mod);
		conv(ZZ_pX_f, f);
		vec_GF2E_modulus.push_back(ZZ_pXModulus(ZZ_pX_f));
		mod <<= 1;
	}

	bak.restore();
}

const ZZ_pXModulus& Constants::loadGF2EModVec(long deg) {
	if (deg - 1 >= (long) vec_GF2E_modulus.size()) {
		std::cout
				<< "Error: Tried to access GF2E_mod vector at invalid precision "
				<< deg << ". size = " << vec_GF2E_modulus.size() << std::endl;
		exit(1);
	}
	if (deg <= 0) {
		std::cout
				<< "Error: tried to access GF2E_mod vector at precision 0 or less"
				<< std::endl;
		exit(1);
	}
	return vec_GF2E_modulus[deg - 1];
}

void Constants::storePowersOf2(long n) {
	powersOf2 = vector<ZZ>(0);
	powersOf2.reserve(n + 1);
	ZZ currPow(NTL::INIT_VAL, 1);

	for (int i = 0; i <= n; i++) {
		powersOf2.push_back(currPow);
		NTL::mul(currPow, currPow, 2);
	}
}

const ZZ& Constants::getPowerOf2(long n) {
	if (n >= (long) powersOf2.size()) {
		std::cout
				<< "Error: Tried to access powers_of_2 vector at invalid power "
				<< n << ". size = " << powersOf2.size() << std::endl;
		exit(1);
	}
	if (n < 0) {
		std::cout
				<< "Error: tried to access GF2E_mod vector at precision 0 or less"
				<< std::endl;
		exit(1);
	}
	return powersOf2[n];
}

