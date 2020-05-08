/*
 * EntropyScorer.h
 *
 *  Created on: Sep 5, 2019
 *      Author: shailesh
 */

#ifndef SRC_ENTROPYSCORER_H_
#define SRC_ENTROPYSCORER_H_

//#include "Inputs.h"
#include "main_common.h"

class EntropyScorer {
	Inputs &inputs;

public:
	EntropyScorer(Inputs &inps) :
			inputs(inps) {
	}
	void run();
};

#endif /* SRC_ENTROPYSCORER_H_ */
