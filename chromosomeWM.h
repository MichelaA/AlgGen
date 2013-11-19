/*
 * chromosomeWM.h
 *
 *  Created on: May 27, 2010
 *      Author: miki
 */

#ifndef CHROMOSOMEWM_H_
#define CHROMOSOMEWM_H_

#include "chromosome.h"

class chromosomeWM: protected chromosome
{
	struct gen{
		int index;
		int* vectAnt;
	};

	gen* vettR;
	int dimWM; //number of rules of the WM rule base


public:
	chromosomeWM(int,unsigned int**);
	void elimina();
	virtual ~chromosomeWM();
};

#endif /* CHROMOSOMEWM_H_ */
