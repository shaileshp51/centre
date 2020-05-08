/*
 * CentreIException.cpp
 *
 *  Created on: Apr 21, 2020
 *      Author: shailesh
 */

#include "CentreIException.h"

/*
 * Declaring a static member variable is not enough, you need to also define it somewhere.
 * This should be in a .cpp file that includes the header (.h) file that declared them.
 */

int IWarning::warning_id = 0;
std::ostringstream cwarn;

