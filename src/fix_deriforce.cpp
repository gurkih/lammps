#include "string.h"
#include "stdlib.h"
#include "math.h"
#include "fix_deriforce.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "input.h"
#include "variable.h"
#include "group.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "error.h"
#include "memory.h"

#include "mpi.h"
#include "assert.h"
#include "float.h" // needed for dbl_max

using namespace LAMMPS_NS;
using namespace FixConst;

ForceDerivative::ForceDerivative(LAMMPS *lmp, int narg, char **arg) 
: Fix(lmp, narg, arg) {
  if (narg < 4) error->all(FLERR,"Illegal fix print command");

  nevery = atoi(arg[3]);
  if (nevery <= 0) error->all(FLERR,"Illegal fix print command");
	memory->create(this->lastf, atom->nmax, 3, "ForceDerivative:lastf");
	atom->add_callback(0);
}

ForceDerivative::~ForceDerivative() {
	memory->destroy(lastf);
	atom->delete_callback(id,0);
}

int ForceDerivative::setmask() {
	int mask = 0;
	mask |= FixConst::END_OF_STEP;
	return mask;
}

void ForceDerivative::sort(int myindex) {
	//this doesnt work as sort is unable to access indicesofblah
	printf("sorting ... \n");
	double blah = indicesofclosestatoms[0][0];
//	double blab = distancesofclosestatoms[0][0];
	printf("i did blah, myindex is %d %f \n", myindex, blah);
	if (distancesofclosestatoms [myindex][2] < distancesofclosestatoms[myindex][0]) {

		printf("1 \n");
		int itmp = indicesofclosestatoms[myindex][2];
		double dtmp = distancesofclosestatoms[myindex][2];

		printf("2 \n");
		indicesofclosestatoms[myindex][2] = indicesofclosestatoms[myindex][1];
		distancesofclosestatoms[myindex][2]= distancesofclosestatoms[myindex][1];
		
		printf("3 \n");
		indicesofclosestatoms[myindex][1] = indicesofclosestatoms[myindex][0];
		distancesofclosestatoms[myindex][1]=  distancesofclosestatoms[myindex][0];

		printf("4 \n");
		indicesofclosestatoms[myindex][0] = itmp;
		distancesofclosestatoms[myindex][0]= dtmp;
	}
	
	else if (distancesofclosestatoms [myindex][2] < distancesofclosestatoms[myindex][1]) {
		
		int itmp = indicesofclosestatoms[myindex][2];
		double dtmp = distancesofclosestatoms[myindex][2];

		indicesofclosestatoms[myindex][2] = indicesofclosestatoms[myindex][1];
		distancesofclosestatoms[myindex][2]= indicesofclosestatoms [myindex][1];

		indicesofclosestatoms[myindex][1] = itmp;
		distancesofclosestatoms[myindex][1]= dtmp;

	}
	printf("finished sorting! \n");
}

double ForceDerivative::euclideandistance(double* firstatom, double* secondatom) {
				printf("measuring ... \n");
				if (firstatom[0] == secondatom[0] && firstatom[1] == secondatom[1] && firstatom[2] == secondatom[2]) {
					printf("done measuring! \n");
						return DBL_MAX; //return inf if we got handed over the same atom twice
				}
				//printf("%f", tmp[0]);
				/*
				double xdistance = pow((poscopy[firstatom][0]-poscopy[secondatom][0]),2);
				double ydistance = pow((poscopy[firstatom][1]-poscopy[secondatom][1]),2);
				double zdistance = pow((poscopy[firstatom][2]-poscopy[secondatom][2]),2);
				*/
				double xdistance = pow(firstatom[0]-secondatom[0],2);
				double ydistance = pow(firstatom[1]-secondatom[1],2);
				double zdistance = pow(firstatom[2]-secondatom[2],2);
				double mydistance = sqrt(xdistance+ydistance+zdistance);
				printf("done measuring! \n");
				return mydistance;
}



double **indicesofclosestatoms;
double **distancesofclostesatoms;

void ForceDerivative::end_of_step() {


	int nlocal = atom->nlocal;
	//indicesofclosestatoms = new double[9999][3]; // just looking for the three closest atoms right now. this might change.
	//distancesofclosestatoms = new double[9999][3];
	

//	extern double** indicesofclosestatoms;
//	extern double** distancesofclostesatoms;
	indicesofclosestatoms = new double*[3];
	distancesofclosestatoms = new double*[3];
	double foo[3] = {0,0,0};
	for (int i = 0; i < nlocal; i++) {
		indicesofclosestatoms[i] = foo;
		distancesofclosestatoms[i] = foo;
	}
//	indicesofclosestatoms[0][0] = 42;
	printf("i just did stuff\n");

	double deriforceoutput[nlocal][3];
	double forceoutput[nlocal][3];
	double speedoutput[nlocal][3];

	float averagedenominator[3] = {0,0,0};

	float averagederif[3] = {0,0,0};
//	printf("i am being run!");


//	printf("; \n");
	for (int indexOfParticle = 0; indexOfParticle < nlocal; ++indexOfParticle) {
	
	double **poscopy = atom->x;
	for(int i = 0; i < 10; i++) {
		for (int j = 0; j < 3; j++) {
			printf("poscopy %d/%d = %f ",i, j, poscopy[i][j]);
		}
		printf("\n");
	}
	double **forcecopy = atom->f;
	double **speedcopy = atom->v;
	int** specialcopy = atom->bond_type;
	bool ifoundsomething = false;
//	printf("i managed to copy!\n");
//	printf("my bond type is: %d \n",atom->bond_type[2][2]);

	indicesofclosestatoms[indexOfParticle][0];
		if (atom->mask[indexOfParticle] & groupbit) {
		//	printf("Hello World");
			
			for (int i = 0; i < 3; i++) {
				indicesofclosestatoms [indexOfParticle][i] = i;
				if (indexOfParticle <=2) {
//					indicesofclosestatoms [indexOfParticle][i]+=3;
				}	
				printf("i managed to init! \n");

// we have to create local copys to avoid some scope trouble.

				double first [3];
				double second [3];
				for (int x = 0; x < 3; x++) {
					first[x] = poscopy[indexOfParticle][x];
					int foo = indicesofclosestatoms[indexOfParticle][i];
					second[x] = poscopy[foo][x];
				}
				//distancesofclosestatoms[indexOfParticle][i] = euclideandistance(indexOfParticle,indicesofclosestatoms [indexOfParticle][i]);i
				distancesofclosestatoms[indexOfParticle][i] = euclideandistance(first, second);
				
			}
			printf("my distance is %f \n", distancesofclosestatoms[indexOfParticle][0]);

			sort(indexOfParticle);

			for (int i = 0; i < nlocal; i++) {
				ifoundsomething = false;
				if (i == indexOfParticle) {
					break;
				}

// we have to create local copys to avoid some scope trouble.

				double first [3];
				double second [3];
				for (int x = 0; x < 3; x++) {
					first[x] = poscopy[indexOfParticle][x];
					second[x] = poscopy[i][x];
				}

				double mydistance = euclideandistance(first, second);

				if (mydistance < distancesofclosestatoms[indexOfParticle][2]) {
					indicesofclosestatoms[indexOfParticle][2] = i;
					distancesofclosestatoms[indexOfParticle][2] = mydistance;
					ifoundsomething = true;
				}
				if (ifoundsomething) {				
					sort(indexOfParticle);
				}
			}
//			printf("%d, %d, %d; \n",poscopy[indexOfParticle][0],poscopy[indexOfParticle][1],poscopy[indexOfParticle][2]);
		}
	}
}

double ForceDerivative::memory_usage() {
	int nmax = atom->nmax;
	double bytes = 0.0;
	if (this->lastf !=0) {
		bytes +=nmax * 3 *sizeof(double);
	}
	return bytes;	
}

void ForceDerivative::grow_arrays(int nmax) {
	if (this->lastf != 0) {
		memory->grow(this->lastf, nmax, 3, "ForceDerivative:lastf");
	}
}

void ForceDerivative::copy_arrays(int i, int j) {
	if (this->lastf != 0) {
		memcpy(this->lastf[j], this->lastf[i],sizeof(double)*3);
	}
}

void ForceDerivative::set_arrays(int i) {
	if (this->lastf != 0) {
		memset(this->lastf[i], 0, sizeof(double) * 3);
	}
}
