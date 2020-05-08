/*
 * Exceptions.h
 *
 *  Created on: Apr 21, 2020
 *      Author: shailesh
 */

#include "CentreIException.h"
#include <cstdlib>

using namespace std;
class NCAttibuteTextError: public IException {
public:
	NCAttibuteTextError(const string &message) {
		cerr << message << endl;
	}
};

class NCAttibuteLengthError: public IException {
public:
	NCAttibuteLengthError(const string &message) {
		cerr << message << endl;
	}
};

class NCConventionReadError: public IException {
public:
	NCConventionReadError(const string &message) {
		cerr << message << endl;
	}
};

class NCUnknownConventionError: public IException {
public:
	NCUnknownConventionError(const string &message) {
		cerr << message << endl;
	}
};

class NCDimensionIDReadError: public IException {
public:
	NCDimensionIDReadError(const string &message) {
		cerr << message << endl;
	}
};

class NCDimensionLengthReadError: public IException {
public:
	NCDimensionLengthReadError(const string &message) {
		cerr << message << endl;
	}
};

class NCGenericLibraryError: public IException {
public:
	NCGenericLibraryError(const string &message) {
		cerr << message << endl;
	}
};

class NCOpenForReadError: public IException {
public:
	NCOpenForReadError(const string &message) {
		cerr << message << endl;
	}
};

class NCOpenForWriteError: public IException {
public:
	NCOpenForWriteError(const string &message) {
		cerr << message << endl;
	}
};

class NCDataReadError: public IException {
public:
	NCDataReadError(const string &message) {
		cerr << message << endl;
	}
};

class NCDataWriteError: public IException {
public:
	NCDataWriteError(const string &message) {
		cerr << message << endl;
	}
};

class NCAttributeCreateError: public IException {
public:
	NCAttributeCreateError(const string &message) {
		cerr << message << endl;
	}
};

class NCDimensionCreateError: public IException {
public:
	NCDimensionCreateError(const string &message) {
		cerr << message << endl;
	}
};

class NCVariableCreateError: public IException {
public:
	NCVariableCreateError(const string &message) {
		cerr << message << endl;
	}
};
