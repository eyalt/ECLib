/*
 * CryptoLib.hpp
 *
 *  Created on: Jul 28, 2014
 *      Author: Itay
 */

#ifndef CRYPTOLIB_HPP_
#define CRYPTOLIB_HPP_

#include "StrongEC.hpp"
#include <NTL/ZZ_p.h>

class ECDHClient {
private:
	StrongEC ec;
	ZZ pvKey;
	ProjectivePoint pubKey; /* G*pvKey */
	ProjectivePoint sharedSecret;

public:
	ECDHClient(const StrongEC& _ec, const ZZ& _pvKey = NTL::conv<ZZ>(0));

	const ZZ& getPvKey();
	const ProjectivePoint& getPubKey();
	const ProjectivePoint& getSharedSecret();
	const StrongEC& getCurve();

	void calcSharedSecret(const ProjectivePoint& pubKey2);

	void setCurve(const StrongEC _ec);
	void setPvKey(const ZZ& _pvKey);
	void setGenerator(const ProjectivePoint& g);

};

/**
 * Performs ECDH between two clients and initialize them accordingly.
 * Both clients must be built in advance and agree on the same StrongEC and generator.
 */
void ecdh(ECDHClient& client1, ECDHClient& client2);

#endif /* CRYPTOLIB_HPP_ */
