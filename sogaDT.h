/*
 * sogaDT.h
 *
 *  Created on: May 20, 2010
 *      Author: miki
 */

#ifndef SOGADT_H_
#define SOGADT_H_

/*
    Copyright (c) <year>, <copyright holder>
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:
        * Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.
        * Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.
        * Neither the name of the <organization> nor the
        names of its contributors may be used to endorse or promote products
        derived from this software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY <copyright holder> ''AS IS'' AND ANY
    EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL <copyright holder> BE LIABLE FOR ANY
    DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
    ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef SOGADT_H
#define SOGADT_H

#include "evolutionDT.h"
#include "frbs.h"


class sogaDT:public EvolutionDT
{
	int numElit;///< Number of elitist solutions
	//dataset* m;///< Bi-dimensional vector of chromosomes. It is used to store the offspring solutions
	int dimPopol;

	dataset* parent_pop;	///< chromosome vector for the parent population
    dataset* child_pop; 	///< chromosome vector for the child  population
    dataset* mixed_pop;	///< chromosome vector for the merged population

public:
	sogaDT(int,int);
	~sogaDT();
	void evaluate_pop(chromosome* arcV, int archlen,frbs&,int );
	void evol(chromosome* arcV, int archlen, frbs&,int);
	void extractPopulation ();///< Updates the current population substituting the worst numElit solution with the best numElit offspring
	int tournament(dataset&, dataset&);
	void generate_child();
	//dataset& getBest(){return parent_pop[0];}
	int getDimPopol(){return dimPopol;}
	dataset& getDataset(int i){return parent_pop[i];}
	dataset* getDt(int i){return &parent_pop[i];}
	void writePop(int);
};

#endif // SOGADT_H


#endif /* SOGADT_H_ */
