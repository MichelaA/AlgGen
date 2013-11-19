/*
 * crom.cpp
 *
 *  Created on: Oct 12, 2012
 *      Author: miki
 */

#include "crom.h"

crom::crom(int num)
{
	genera(num);
}

void crom::genera(int nume)
{	gen=new int[nume];
	num=nume;
	fitness=0;
	nregole=0;

}

crom::~crom() {
	delete[] gen;
}

void crom::copiaCr(const crom& cr)
{	for (int i=0;i<cr.num;i++)
		gen[i]=cr.gen[i];
	fitness=cr.fitness;
	nregole=cr.nregole;
}


crom::crom(const crom& cr)
{	gen=new int[cr.num];
	copiaCr(cr);
}

crom& crom::operator=(const crom& cr)
{	if (gen!=0)
		delete[] gen;
	gen=new int[cr.num];
	copiaCr(cr);
	return *this;
}


