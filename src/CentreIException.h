/*
 * CentreIException.h
 *
 *  Created on: Apr 21, 2020
 *      Author: shailesh
 */

#ifndef SRC_CENTREIEXCEPTION_H_
#define SRC_CENTREIEXCEPTION_H_

#include <iostream>
#include <sstream>

using namespace std;
extern ostringstream cwarn;

/**
 * class CException defines the format of exceptions shown in the log file.
 *
 * Exception is a severe problem that the program does not know how to correct it. Program terminates
 * immediately with a message when an exception is thrown.
 */
class IException {
public:
	IException() {
		cerr << "[FATAL ERROR] ";
	}

	~IException() {
		cerr << "\n";
	}
};

/**
 * class CWarning defines the format of warnings shown in the log file.
 *
 * Warning is a problem that the program can correct (usually reset the variable to its default value)
 * for the user without terminating the current run. Instead, a message is shown to remind the user
 * what has been corrected.
 */

#ifndef OMP_MPI
class IWarning // Program continues its run with warnings
{
private:

public:
	static int warning_id;

	IWarning() {
		warning_id++;

		cwarn << "[WARNING #" << warning_id << "] ";
	}

	~IWarning() {
		cwarn << "\n";
	}
};
#endif

#ifdef OMP_MPI
class IWarning // Program continues its run with warnings
{
   private:

   public:
      static int warning_id;
      IWarning()
      {
        warning_id++;

         cwarn << "   " << "\033[1;34m" << "[WARNING #" << warning_id << "] ";
      }

      ~IWarning()
      {
         cwarn << "\033[0m";
      }
};
#endif

#endif /* SRC_CENTREIEXCEPTION_H_ */
