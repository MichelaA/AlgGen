/*
 * crom.h
 *
 *  Created on: Oct 12, 2012
 *      Author: miki
 */

#ifndef CROM_H_
#define CROM_H_

class crom {
	int* gen;
	double fitness;
	int num;
	double nregole;
	void copiaCr(const crom& );


public:
	crom(int);
	crom(){};
	crom(const crom&);
	~crom();
	void genera(int );
	crom& operator=(const crom&);

	int* getGen(){return gen;}
	double getFitness(){return fitness;}
	void setGen(int* g){for (int i=0;i<num;i++) gen[i]=g[i];}
	void setFitness(double f){fitness=f;}
	void setNregole(int n){nregole=n;}
	int getNregole(){return nregole;}
};

#endif /* CROM_H_ */
