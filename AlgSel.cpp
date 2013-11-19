
/*
 * AlgGenNew.cpp
 *
 *  Created on: Dec 17, 2010
 *      Author: miki
 *
 */

//============================================================================
// Name        : AlgGenRed.cpp
// Author      : Miki
// Version     :
// Copyright   : Your copyright notice
// Description : Ansi-style
//============================================================================

/***************************************************************************
 *   Copyright (C) 2008 by miki   *
 *   michela.antonelli@iet.unipi.it   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <fstream>

#include "paes.h"
#include "sogaDT.h"
#include "dataset.h"
#include "ctime"
#include "soga.h"


using namespace std;
extern int* numPart;
extern int** matC45;
extern int dimmatC45;


int leggi(vector<int>* matrice)
{	ifstream infile;
	infile.open ("../Indicifeat.txt");
	if (!infile)
	{		cout<<"Errore"<<endl;
		exit(1);
	}
	int numero;
	int linea=0;
	int colonna=0;
	string s;

	while(getline(infile,s))
	{	colonna=0;
		istringstream iss(s);
	    while (iss>>numero)
	    {	matrice[linea].resize(matrice[linea].size()+1);
	    	matrice[linea][colonna++]=numero;
	    }
		linea++;
	}

}

void leggiGran(vector<int>* gran,int numfold,double** miV, double** maV)
{
	ifstream f("gran.txt");
	if (!f)
	{	cout<<"errore"<<endl;
		exit(1);
	}
	string stringa,sub,sub1;

	for (int i=0;i<numfold;i++)
	{	getline(f,stringa);
		cout<<stringa<<endl;
		istringstream iss(stringa);
		int j=0;
		while(iss>>sub)
		{	gran[i].resize(gran[i].size()+1);
		  	gran[i][j]=atoi(sub.c_str());
		  	j++;
		}
	/*	getline(f,stringa);
		istringstream iss1(stringa);
		miV[i]=new double[j];
		maV[i]=new double[j];
		j=0;
		while(iss1>>sub)
		{	iss1>>sub1;
		    istringstream iss2(sub);
		    istringstream iss3(sub1);
			iss2>>miV[i][j];
			iss3>>maV[i][j];
			j++;
		}*/
		miV[i]=0;
		maV[i]=0;
	}

	/*for (int i=0;i<numfold;i++)
			  {	  for (int j=0;j<gran[i].size();j++)
					  cout<<gran[i][j]<<' ';
			  	  cout<<endl;
		  }*/
}

int main ( int argc, char *argv[] )
{
	  time_t start,end;
	  double dif;
	  time (&start);

	 // vector<int>* matFeat=0;
	  int numfold=10;

	//  double** miV, **maV;

	/*  if (!CLASSIFICATION)
	  {	  //inizializeVar ( "GLAT11.ENT","pippo", 0 );
		  // inizializeVar(argv[1],argv[2],0);
		  vector<int>*  gran;
		  gran = new vector<int>[numfold];
		  miV=new double*[numfold];
		  maV=new double*[numfold];

		 // leggiGran(gran,numfold,miV,maV);



	/*	  for (int i=0;i<numfold;i++)
		  {	  for (int j=0;j<gran[i].size();j++)
				  cout<<gran[i][j]<<' ';
		  	  cout<<endl;
	  }*/
/*		  for (int i=0;i<numfold;i++)
		  	  {	 miV[i]=maV[i]=0;
		  	  //inizializeVar ( "GLAT11.ENT","pippo", i); //, gran[i],miV[i],maV[i]);//

		  	  	  inizializeVar(argv[1],argv[2],i); //,gran[i],0,0);
		  	  /*if (matFeat==0)
		  	  {		matFeat=new vector<int>[numfold];
	  				for (int j=0;j<numfold;j++)
	  					matFeat[j].reserve(numVar);
	  				leggi(matFeat);
		  	  }
		  	  calcolaWMrandom(matFeat[i]);*/

/*		  seleziona();
	  	  }
		   	   exit(1);
	  	  }

	  	  if (CLASSIFICATION)
	  	  {	soga sg;
	  	  	sg.evolution();
	  	  	sg.setc45Best();
	  	  }
*/
//bool fatto=false;
//for (int i=0;i<numfold;i++)
	//inizializeVar( "GLAT11.ENT","pippo", i);
	int i=0;
	inizializeVar(argv[1],argv[2],i);
	if (CLASSIFICATION)
	{
		int numAtt=0;
	 	frbs fis(numVar);
	 	int* indici=0;
		matC45=calcolaMatC45(dimmatC45,numAtt,indici);
		cambiaVariabili(numAtt,indici,numParts);
		delete[] indici;
	 }
/*	for (int i=0;i<5;i++)
	{	//inizializeVar(argv[1],argv[2],i);
		inizializeVar ( "GLAT11.ENT","pippo",i);
*/
	paes popol(4,2);

	if (Perc==100)
		 popol.evolvPop(inOutTr,0,numPatterTr,1);
	 else
	 {	 char nomefile[100]="EvolDataset.txt";
		 fstream fpObj,fpObj1;
		 fpObj.open(nomefile, ios::out | ios::app);
		 fpObj<<"0 0 0 0 0 0 0 ";
		 popol.evolvRed();
	 }

	 time (&end);
	 dif = difftime (end,start);
	 printf ("It took  %.6lf seconds.\n", dif );

	 char nomefile1[]="tempo.txt";
	 fstream fptempo;
	 fptempo.open(nomefile1,ios::out | ios::app);
	 if (!fptempo)
	 {	cout<<("Error opening file\n");
		 exit(1);
	 }
	 fptempo<<dif<<endl;

	/* delete[] minVal;
	 delete[] maxVal;*/

	 for ( int j=0;j<numPatterTr;j++ )
	 	 delete[] inOutTr[j];

	 for ( int j=0;j<numPatterTs;j++ )
		 delete[] inOutTs[j];

	 delete[] inOutTr;
	 delete[] inOutTs;
}
	 //return EXIT_SUCCESS;




