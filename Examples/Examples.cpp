#include <CryptoLib.hpp>

#include <NTL/GF2XFactoring.h>
#include <sstream>
#include <time.h>

int main() {
	//initialize NTL seed to ensure randomness
	NTL::SetSeed(ZZ(INIT_VAL, clock()));

	//Build irreducible polynomial
	NTL::GF2X P;
	NTL::BuildIrred(P, 127);
	//Initialize GF2E modulus
	GF2E::init(P);
	//Initialize constants to fit the new modulus
	Constants::initConstants();

	//Creating GF2E polynomial from a string
	std::stringstream ss;
	GF2E b;
	ss << "[0 1 1 0 1]";
	ss >> b;
	//Creating an elliptic curve with b
	EllipticCurve e(b);
	//Printing the b parameter of the curve
	std::cout << "b = " << e.getB() << std::endl;

	//Creating a custom strong EC
	StrongEC sec1(e);
	std::cout << "order = " << sec1.getOrder() << std::endl;

	//Creating a random (truly Strong) strong EC
	StrongEC sec2;

	//Printing the curve's details
	std::cout << "b = " << sec2.getCurve().getB() << std::endl;
	std::cout << "order = " << sec2.getOrder() << std::endl;
	std::cout << "strong = " << (sec2.isStrong() ? "True" : "False")
			<< std::endl;
	std::cout << "p = " << sec2.getP() << std::endl;

	//Generate a strong (of order p) point from a strong EC
	StrongEC e1;
	ProjectivePoint p1 = e1.genStrongPoint();
	ProjectivePoint p2 = e1.genStrongPoint();

	//Point arithmetics
	ProjectivePoint p_res;
	p_res = e1.pointAdd(p1, p2); // p_res = p1+p2
	if (!p_res.isINF()) {
		p_res = p_res.toAffine(); // convert p_res to affine
	}
	p_res = e1.pointDoubling(p1); // p_res = 2p1
	p_res = e1.pointMul(p1, ZZ(INIT_VAL, 6)); // p_res = 6p1
	if (!e1.getCurve().isOn(p_res)) {
		std::cout << "result is not on the curve" << std::endl;
	}

	//A generator is chosen in random by default, thus
	//both clients agree on a the curve and generator.
	StrongEC e2;
	//deciding on both private keys.
	NTL::ZZ pv_key_alice = ZZ(INIT_VAL, 1337);
	NTL::ZZ pv_key_bob = ZZ(INIT_VAL, 9001);

	ECDHClient alice = ECDHClient(e2, pv_key_alice);
	ECDHClient bob = ECDHClient(e2, pv_key_bob);
	//agreeing on the shared secret (1337*9001*generator).
	ecdh(alice, bob);
	std::cout << "Shared secret: " << alice.getSharedSecret() << std::endl;

	//switching the generator for ecdh. (notice that the curve is still e)
	ProjectivePoint new_generator = e2.genStrongPoint();
	bob.setGenerator(new_generator);
	alice.setGenerator(new_generator);
	ecdh(alice, bob);
	std::cout << "Shared secret: " << alice.getSharedSecret() << std::endl;

	return 0;
}
