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

using namespace LAMMPS_NS;
using namespace FixConst;

bool ForceDerivative::beenhere = false;
double ** ForceDerivative::lastf;
//double** ForceDerivative::lastf;

ForceDerivative::ForceDerivative(LAMMPS *lmp, int narg, char **arg) 
: Fix(lmp, narg, arg) {
	int nlocal = atom->nlocal;
	//if(update->ntimestep == 1) {
		lastf = new double*[atom->nlocal];
		for (int i = 0; i < nlocal; i++) {
			lastf[i] = new double[3];
		}
		printf("i just initialized lastf \n");
	//}

  if (narg < 4) error->all(FLERR,"Illegal fix print command");

  nevery = atoi(arg[3]);
  if (nevery <= 0) error->all(FLERR,"Illegal fix print command");
	memory->create(this->lastf, atom->nmax, 3, "ForceDerivative:lastf");
	atom->add_callback(0);

}

ForceDerivative::~ForceDerivative() {
//	memory->destroy(lastf);
	atom->delete_callback(id,0);
}

int ForceDerivative::setmask() {
	int mask = 0;
	mask |= FixConst::END_OF_STEP;
	return mask;
}

void ForceDerivative::end_of_step() {
	int nlocal = atom->nlocal;
	double deriforceoutput[nlocal][3];
	double forceoutput[nlocal][3];
	double speedoutput[nlocal][3];

	float averagedenominator[3] = {0,0,0};

//	printf("timestep is %lu, my timestep is %lu \n",update->ntimestep, this->mytimestep);
	float averagederif[3] = {0,0,0};
//	printf("i am being run!");

	double **forcecopy = atom->f;
	double **speedcopy = atom->v;
	if(beenhere == true) { //act, this should be false. but it works this way. TODO: find out why
		printf("i am writing lastfs \n");

		printf("false. my value is %d \n ",lastf[nlocal-2][0]);
		for (int indexOfParticle = 0; indexOfParticle < nlocal; ++indexOfParticle) {
			for (int i = 0; i < 3; i++) {
				lastf[indexOfParticle][i] = forcecopy[indexOfParticle][i];
			}
		}
		printf("i just wrote the value %f \n", lastf[nlocal-2][0]);
		beenhere = false; //this should be true
	} else {
//		printf("plah \n");
		beenhere = true; //this should be false
		printf("true. my value is %d \n ",lastf[nlocal-2][0]);
//		this->mytimestep++;
		int tmp = 0;
		double add;
		double average = 0;
		for (int indexOfParticle = 0; indexOfParticle < nlocal; ++indexOfParticle) {	
			if (atom->mask[indexOfParticle] & groupbit) {
				if (lastf[indexOfParticle][2] == 0)
					tmp++;
				}
				for (int i = 0; i < 3; i++) {
//					printf("lastf = %f, actualf = %f \n", lastf[indexOfParticle][i],forcecopy[indexOfParticle][i]);
					add =this->lastf[indexOfParticle][i] - forcecopy[indexOfParticle][i];
					add /=nlocal;
					add /=1.01;
				}
				average+=add;
			}
		double deri = lastf[1][0]-forcecopy[1][0];
		printf("the number of 0 forces is %d, the derivation is %f. \n",tmp, average);
		for (int i = 0 ; i < nlocal; ++i) {
			//printf("%f, %f, %f last %f, %f, %f \n", forcecopy[i][0],forcecopy[i][1],forcecopy[i][2],lastf[i][0],lastf[i][1],lastf[i][2]);
		}
		printf("deri is %f \n",deri);
		
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

