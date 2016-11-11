/*
 * AffineException.hpp
 *
 *  Created on: 8 αιεπ 2014
 *      Author: Eyal
 */

#ifndef AFFINEEXCEPTION_HPP_
#define AFFINEEXCEPTION_HPP_

class AffineException: public std::exception {
	const char* what() {
		return "Must use affine point";
	}
};

class SetGeneratorException: public std::exception {
	const char* what() {
		return "The point isn't a strong point";
	}
};

class GeneratorNotSetException: public std::exception {
	const char* what() {
		return "A generator was yet to be set";
	}
};

class CurvesDontMatch: public std::exception {
	const char* what() {
		return "the curves of the clients don't match";
	}
};

#endif /* AFFINEEXCEPTION_HPP_ */
