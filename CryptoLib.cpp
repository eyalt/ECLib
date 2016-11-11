/*
 * CryptoLib.cpp
 *
 *  Created on: Jul 28, 2014
 *      Author: Itay
 */

#include "CryptoLib.hpp"
#include "EC_Exceptions.hpp"
#include <assert.h>

ECDHClient::ECDHClient(const StrongEC& _ec, const ZZ& _pvKey) :
		ec(_ec), pvKey(_pvKey) {
	if (NTL::IsZero(_pvKey)) {
		NTL::ZZ_p rand;

		NTL::ZZ_pBak temp;
		temp.save();
		NTL::ZZ_p::init(ec.getP());

		NTL::random(rand);
		conv(pvKey, rand);

		temp.restore();
	}

	pubKey = ec.pointMul(ec.getGenerator(), pvKey);
}

const ZZ& ECDHClient::getPvKey() {
	return pvKey;
}

const ProjectivePoint& ECDHClient::getPubKey() {
	return pubKey;
}

const ProjectivePoint& ECDHClient::getSharedSecret() {
	return sharedSecret;
}

const StrongEC& ECDHClient::getCurve() {
	return ec;
}

void ECDHClient::calcSharedSecret(const ProjectivePoint& pubKey2) {
	sharedSecret = ec.pointMul(pubKey2, pvKey);
}

void ECDHClient::setCurve(const StrongEC _ec) {
	ec = _ec;
	if (!ec.hasGenerator()) {
		ec.setGenerator();
	}
	pubKey = ec.pointMul(ec.getGenerator(), pvKey);
}

void ECDHClient::setPvKey(const ZZ& _pvKey) {
	pvKey = _pvKey;
	pubKey = ec.pointMul(ec.getGenerator(), pvKey);
}

void ECDHClient::setGenerator(const ProjectivePoint& g) {
	ec.setGenerator(g);
	pubKey = ec.pointMul(ec.getGenerator(), pvKey);
}

void ecdh(ECDHClient& client1, ECDHClient& client2) {
	if (client1.getCurve() != client2.getCurve()) {
		throw CurvesDontMatch();
	}

	client1.calcSharedSecret(client2.getPubKey());
	client2.calcSharedSecret(client1.getPubKey());

	assert(client1.getSharedSecret() == client2.getSharedSecret());
}

