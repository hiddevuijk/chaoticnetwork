// nr_headers (nr_headers/)
#include "nr_headers/nr3.h"
#include "nr_headers/ran.h"

// std headers
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <math.h>
#include <time.h>

//other headers (headers/)
#include "headers/generate_matrix.h"
#include "headers/write_matrix.h"

using namespace::std;


typedef std::vector<std::vector<int> > Matint;
int main()
{
	clock_t time_secs = clock();
	int time_1, time_2;
	ofstream log("log.txt");

	int N = 100;
	int Ne = N;
	int Ni = N;
	int No = N;

	int K = 10;

	double Jeo = 1.0;
	double Jee = 1.0;
	double Jie = 1.0;
	double Jio = 0.8;
	double Jei = -0.2;
	double Jii = -1.8;

	double se = 1.0;
	double si = 0.8;
	double D = 0.3;

	int tmax = 10;
	double tau = 0.9;

	int seed = 123456789;
	Ranq1 r(seed);

	
	time_1 = clock();
	Matint EE = connectivity_matrix(Ne,Ne,K,r);
	Matint EI = connectivity_matrix(Ne,Ni,K,r);
	Matint IE = connectivity_matrix(Ni,Ne,K,r);
	Matint II = connectivity_matrix(Ni,Ni,K,r);
	Matint EO = connectivity_matrix(Ne,No,K,r);
	Matint IO = connectivity_matrix(Ni,No,K,r);
	time_2 = clock();
	log << "initializen connection matrices took: " << (time_2-time_1)/(double)CLOCKS_PER_SEC << endl;
	vector<double> the(Ne,0.0);
	vector<double> thi(Ni,0.0);
	double sqK = sqrt(1.*K);
	for(int i=0;i<Ne;i++) the[i] = (se + r.doub()*D)*sqK;
	for(int i=0;i<Ni;i++) thi[i] = (si + r.doub()*D)*sqK;


	vector<int> nwe(Ne,0);
	vector<int> nwi(Ni,0);
	vector<int> nwo(No,0);
//	for(int i=0;i<Ne;i++) if(0.5<r.doub()) nwe[i] = 1;
//	for(int i=0;i<Ni;i++) if(0.5<r.doub()) nwi[i] = 1;
	for(int i=0;i<No;i++) if(0.08<r.doub()) nwo[i] = 1;

	vector<int> empty_vecint;
	vector<vector<int> > nwet(tmax,empty_vecint);
	vector<vector<int> > nwit(tmax,empty_vecint);
	nwet[0] = nwe;
	nwit[0] = nwi;

	vector<double> re(Ne,0.0);
	vector<vector<double> > currentEE(tmax,re);
	vector<vector<double> > currentEI(tmax,re);
	vector<double> ri(Ni,0.0);
	vector<vector<double> > currentIE(tmax,ri);
	vector<vector<double> > currentII(tmax,ri);	
	
	// start simulation
	time_1 = clock();
	double input, inputEE, inputEI,inputIE, inputII, inputEO, inputIO;
	int ie, ii;
	for(int t=1;t<tmax;t++) {
		//update exitatory population
//		for(int ti=0;ti<Ne;ti++) {
			ie = r.int64() % (Ne-1);
			inputEE = Jee*dotproduct(EE[ie],nwe,Ne);
			inputEI = Jei*dotproduct(EI[ie],nwi,Ni);
			inputEO = Jeo*dotproduct(EO[ie],nwo,No);
			currentEE[t][ie] = inputEE+inputEO;
			currentEI[t][ie] = inputEI;
			input = inputEE+inputEO+inputEI;
			if(input>the[ie]) nwe[ie] = 1;
			else nwe[ie] = 0;
//		}
		// update inhibitory population
//		for(int ti=0;ti<(Ni/tau);ti++) {
			ii = r.int64() % (Ni-1);
			inputIE = Jie*dotproduct(IE[ii],nwe,Ne);
			inputII = Jii*dotproduct(II[ii],nwi,Ni);
			inputIO = Jio*dotproduct(IO[ii],nwo,No);
			currentIE[t][ii] = inputIE+inputIO;
			currentII[t][ii] = inputII;
			input = inputIE+inputII+inputIO;
			if(input>thi[ii]) nwi[ii] = 1;
			else nwi[ii] = 0;
//		}

		nwet[t] = nwe;
		nwit[t] = nwi;
	}
	time_2 = clock();
	log << "simulation took: " << (time_2-time_1)/(double)CLOCKS_PER_SEC << endl;

	write_matrix(nwet,tmax,Ne,"nwet.csv");
	write_matrix(nwit,tmax,Ni,"nwit.csv");
	write_matrix(nwo,No,"nwo.csv");
	write_matrix(currentEE,tmax,Ne,"currentEE.csv");
	write_matrix(currentEI,tmax,Ne,"currentEI.csv");
	write_matrix(currentIE,tmax,Ni,"currentIE.csv");
	write_matrix(currentII,tmax,Ni,"currentII.csv");
	write_matrix(EE,Ne,Ne,"EE.csv");
	write_matrix(EI,Ne,Ni,"EI.csv");
	write_matrix(IE,Ni,Ne,"IE.csv");
	write_matrix(II,Ni,Ni,"II.csv");
	write_matrix(EO,Ne,No,"EO.csv");
	write_matrix(IO,Ni,No,"IO.csv");
	write_matrix(the,Ne,"the.csv");
	write_matrix(thi,Ni,"thi.csv");
	return 0;
}







