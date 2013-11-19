/*
 * dataset.h
 *
 *  Created on: May 20, 2010
 *      Author: miki
 */

#ifndef DATASET_H_
#define DATASET_H_

/**
	@author miki <michela.antonelli@iet.unipi.it>
*/
#include "chromosome.h"
#include "supportchrom.h"
#include "frbs.h"
#include <fstream>


using namespace std;




class dataset
{
	int* Blocks;
	double* obj;
	int dimChrom;
	int nObj;

	void Elimina();


public:
	int conta();
	dataset(){};
    ~dataset();
	dataset(const dataset&);
	dataset& operator=(const dataset& d);

	bool operator==(const dataset& );
	bool operator!=(const dataset& );

	void inizialize(int,int);
	//double evaluateDataset(double**, int,chromosome&, frbs&);
	friend void cross(dataset&,dataset&,dataset&,dataset&,int);
	void Mutacion(int);
	fstream& writeDT(fstream&,int);

	double* getObj(){return obj;}
	double getObj(int i){return obj[i];}
	int getnObj(){return nObj;}
	void copyDataset(dataset&);
	int getBlock(int j){return Blocks[j];}
	int* getBlocks(){return Blocks;}
	void setObj(double errors,int ind){obj[ind]=errors;}



	int getDim(){return dimChrom;}




};

int compareChromBlock(const  void * _a, const  void * _b);

#endif //DATASET_H_
