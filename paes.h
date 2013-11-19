/*
 * paes.h
 *
 *  Created on: May 20, 2010
 *      Author: miki
 */

#ifndef PAES_H_
#define PAES_H_

#include "evolution.h"
#include "dataset.h"

#include "sogaDT.h"
#include "frbs.h"
#include "chromosomeWM.h"
/**
 * Class for Modified (2+2)PAES algorithm
 */

class paes:evolution
{
	int depth;///< this is the number of recursive subdivisions of the objective space carried out in order to divide the objective space into a grid for the purposes of diversity maintenance. Values between 3 and 6 are useful, depending on number of objectives.
	int archive; ///<archive dimension
	int iteration;///<number of iterations
	int arclengthV; ///< current length of the virtual archive
	int arclengthC; ///< current length of the concrete archive
	double* gl_offset; ///< vector of the minumum values for each objective
	double* gl_range;///< vector of the range values for each objective
	double* gl_largest;///< vector of the maximum values for each objective
	int dimGrid;///<grid dimension
	int *grid_pop;   ///< the array holding the numbero of  residing in each grid location
	chromosome* c;///< vector of the current solutions (2)
	chromosome* m;///< vector of new generated solutions (2)
	chromosome* arcV;///Archive: vector of non-dominatedchromosomes.
	int nsolution; ///< number of current solutions

public:

	void evaluatePop(double**, int);
	void evolvPop(double** inOutRed, int numEvol, int numPatterRed,int);
	paes(int depth,int nsolution);///< Inizializes the grid structures, the current solutions and the archive(contains the current solutions)

	~paes();
	//void evol();///< Performs the genetic evolution
	void evolRB_GR(); ///<Performs the generic evolution using 2 archive, one for the virtual RB and one for the concrete RB
	void evol2();
	void evolvRed(); //Evolution of paes for reduced dataset
	void evolApprox();///< Performs the genetic evolution exploiting some datareduction approaches
	void updateArc(double** inOutRed, int numPatterRed,int numObj);
	void updateArcDaChrom(chromosome* pop,int& dimpop,double** inOutRed, int numPatterRed);


	//void write(chromosome*);
  /**
 	*  \param m chromosome to be compared with the solutions in the archive
	*  \param arc archive to be considered
	*  \param arclength number of solution of the archive
	*  \return -1 if dominated by any member, 1 if dominates any member, and 0 otherwise
*/
	int compare_to_archive(chromosome& m,chromosome* arc,int arclength); ///< compares a solution to every member of the archive.

/**
 	*  \param m chromosome be updated compared with the solutions in the archive
	*  \return -1 if dominated by any member, 1 if dominates any member, and 0 otherwise
*/


/**
 	*  \param s chromosome candidate to be insert into the archive
	*  \param arc archive of solutions
	*  \param arclength number of solution of the archive
*/
	void update_grid(chromosome& s, chromosome* arc, int arclength); ///< Update (if needed) the grid structures  and the location of the archive solutions. Calculate the location of the chromosome s

/**
 * Given a solution s, add it to the archive if
   	a) the archive is empty
   	b) the archive is not full and s is not dominated or equal to anything currently in the archive
   c) s dominates anything in the archive
   d) the archive is full but s is nondominated and is in a no more crowded square than at least one solution.
	 \param s chromosome candidate to be insert into the archive
	\param arc archive to be considered
	\param arclength number of solution of the archive
*/
	void archive_soln(chromosome& s, chromosome* arc, int& arclength,double** matrPoints,int numPatt); ///< Add a solution to the archive following the specific paes strategy

/**
 	*  \param s chromosome candidate to be insert into the archive
	*  \param arc archive of solutions
	*  \param arclength number of solution of the archive
*/
	void add_to_archive(chromosome& s,chromosome* arc, int& arclength); ///< add a solution to the archive updating the length

/**
 	*  \param eval objective vector
*/
	int find_loc(double * eval, int numObj);///< finds the grid location of a solution given a vector of its objective values

	//void save(chromosome* appo,int dim, int num);

	chromosome& getChrom(int k){return arcV[k];}
	int getRealNumRule(int index, chromosome* arc){return arc[index].getRealNumRule();}
	//int getNumRule(int index,chromosome* arc){ return arc[index].getNumRule();}
    //double getObjTot(int index, int num,chromosome* arc){return arc[index].getObjTot(num);}
    //double getECM_ts(int index, frbs& fis, chromosome* arc,bool conc){double* d=0; int comp ; return *arc[index].ECM(comp,fis,inOutTs,numPatterTs,d,true);}




	int getIteration(){return iteration;}
	int getnsolution(){return nsolution;}
	chromosome* getarcV(){return arcV;}
	int getarclengthV(){return arclengthV;}
	frbs& getFis(){return fis;}
	chromosome& getchrom(int i){return arcV[i];}


};



#endif /* PAES_H_ */
