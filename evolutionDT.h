/*
 * evolutionDT.h
 *
 *  Created on: May 20, 2010
 *      Author: miki
 */

#ifndef EVOLUTIONDT_H_
#define EVOLUTIONDT_H_

#include "dataset.h"
#include "supportchrom.h"

class EvolutionDT
{

  protected:

  struct block{
		int* setpoint; //indexes of point belonging to the block
		int dimBlock;
		double errBlock;//error of the block
	};

	block* indTr;     //contains the current indices for the reduced training set

	int numBlock;

	int dimChrom;




public:

	EvolutionDT(int,int);
	~EvolutionDT();
	double** buildMatrix(dataset&,int&);
	double** buildMatrix(int*,int);
	//evolution();
	void crossover(dataset&,dataset&,dataset&,dataset&);
	void inizializePop(dataset* dt, int dimChrom,int dimP,int); ///< Inizializes the chromosome population
	void mutation(dataset& ,double); ///< permforms the mutation operations


	/*int getNumRule(chromosome chrom[], int index){ return chrom[index].getNumRule();} ///< Returns the number of rules of the index-th chromosome
	double getObj(chromosome chrom[], int index, int num){return chrom[index].getObj(num);} ///< Returns the num-th objective of the index-th chromosome
	double getECM_ts(chromosome chrom[], int index, bool conc){double d; return chrom[index].ECM_ts(fis,d,conc);} ///< Returns the error on the testset of the index-th chromosome
	*/

/**
 	*  \param chrom vector of chromosomes to be saved
    *  \param popsize dimension of the population
 	*  \param conc if equal to 1 save the results of the concrete RB  otherwise save the results of the virtual RB
	*
*/
	//void saveresults(chromosome* chrom, int popsize,bool conc);  ///< Save the results of the population on files

/**
 	*  \param arc0 chromosome with the best accuracy
    *  \param it number of iteration
	*  \param arclength archive dimension
	*  \param conc if equal to 1 save the results of the concrete RB  otherwise save the results of the virtual RB
	*
*/
//void stampaBest(chromosome& arc0,int it,int arclength,bool conc); ///< Print the solution with best accuracy


/**
 	*  \param chrom vector of chromosomes to be saved
    *  \param popsize dimension of the population
	*  \param conc if equal to 1 save the results of the concrete RB  otherwise save the results of the virtual RB
	*  \param numEval number of evaluation
	*
*/
//void savestep(chromosome* chrom, int popsize,int numEval,bool conc); ///< Save the population after numEval evaluation


/**
	*  \param chrom population of chromosomes
    *  \param dimPop dimension of the population
*/
//void ordinaPopChrom(chromosome* chrom,int dimPop); ///< Sort the chromosome population by ascendent accuracy




};

#endif // EVOLUTIONDT_H
