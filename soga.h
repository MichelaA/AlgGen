/*
 * soga.h
 *
 *  Created on: Sep 28, 2012
 *      Author: miki
 */

#ifndef SOGA_H_
#define SOGA_H_


#include "supportchrom.h"
#include <sstream>
#include "frbs.h"
#include <map>
#include "crom.h"

class soga {

	struct matGran
	{	int** mat;
		int nRighe;
		int ncolonne;
		double fitness;
	};

	crom* parent_pop;
	crom* child_pop;



	int numGen;
	int numPop;
	int oldnumVar;


	static map<string,matGran> matSal;

	static	double** oldTr;
	static	double** oldTs;
	static int oldMaxterms,oldNumVar,indMat;
	static int* oldNumParts;
	static double* oldMinVal,* oldMaxVal;

	string faiStringa(int*);
public:

	soga();
	soga(const soga&);
	~soga();
	void inizializePop(crom* p,int n);
	void evolution();
	void evaluatePop(crom*);
	int tournament(crom& ind1, crom& ind2);
	void generate_child();
	void ordinaPopol(crom* popol, int n); ///< Orders a vector in ascendent way
	void extractPopulation ();
	void setc45Best();
	void mutation(crom& par);
	void onePointCross(const crom& par1, const crom& par2, crom& suc1,crom& suc2);

	double valutamatrice(int**, double**, int);
	void copiaValori(double**, double**, int, int, int*, double*, double*);
	void ripristinaValori();


};

#endif /* SOGA_H_ */
