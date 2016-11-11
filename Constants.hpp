/*
 * Constants.hpp
 *
 *  Created on: Jul 28, 2014
 *      Author: Itay
 */

#ifndef CONSTANTS_HPP_
#define CONSTANTS_HPP_

#include <vector>

#include <NTL/GF2E.h>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/vec_ZZ.h>

using std::vector;

using NTL::ZZ;
using NTL::INIT_VAL;
using NTL::ZZX;
using NTL::vec_ZZ;
using NTL::ZZ_pXModulus;

class Constants {
private:
	static NTL::ZZX GF2E_modulus;
	static vector<ZZ_pXModulus> vec_GF2E_modulus;
	static vector<ZZ> powersOf2;
public:
	static NTL::GF2E zero;
	static NTL::GF2E one;

	static ZZX ZZX_1;

	static ZZ ZZ_2;
	static ZZ ZZ_12;
	static ZZ ZZ_36;
	static ZZ ZZ_240;
	static ZZ ZZ_504;
	static ZZ ZZ_512;
	static ZZ ZZ_1488;
	static ZZ ZZ_1728;
	static ZZ ZZ_4095;
	static ZZ ZZ_12096;
	static ZZ ZZ_162000;
	static ZZ ZZ_195120;
	static ZZ ZZ_372735;
	static ZZ ZZ_563760;
	static ZZ ZZ_40773375;
	static ZZ ZZ_660960000;
	static ZZ ZZ_8748000000;
	static ZZ ZZ_8981280000;
	static ZZ ZZ_157464000000000;

	static void initConstants();

	static void storeGF2EMod();
	static const NTL::ZZX& loadGF2EMod();
	static void storeGF2EModVec(long maxDeg);
	static const ZZ_pXModulus& loadGF2EModVec(long deg);

	static void storeTestMod(const NTL::ZZ_pX& mod);

	static NTL::GF2X findIrred(int p);

	static void storePowersOf2(long n);
	static const ZZ& getPowerOf2(long n);

};

#endif /* CONSTANTS_HPP_ */
