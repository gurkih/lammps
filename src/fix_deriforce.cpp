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
	//memory->destroy(indicesofclosestatoms);
	//memory->destroy(distancesofclosestatoms);
}

int ForceDerivative::setmask() {
	int mask = 0;
	mask |= FixConst::END_OF_STEP;
	return mask;
}

double ForceDerivative::euclideandistance(double* firstatom, double* secondatom) {
//				printf("measuring ... \n");
				if (firstatom[0] == secondatom[0] && firstatom[1] == secondatom[1] && firstatom[2] == secondatom[2]) {
//					printf("done measuring! \n");
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
//				printf("done measuring! \n");
				return mydistance;
}

void ForceDerivative::end_of_step() {
	int nlocal = atom->nlocal;
	double deriforce[nlocal][3];
//	double forceoutput[nlocal][3];
//	double speedoutput[nlocal][3];
//	float averagedenominator[3] = {0,0,0};
//	float averagederif[3] = {0,0,0};
	double **poscopy = atom->x;
	double **forcecopy = atom->f;
/*
	for(int i = 0; i < 10; i++) {
		for (int j = 0; j < 3; j++) {
			printf("poscopy %d/%d = %f ",i, j, poscopy[i][j]);
		}
		printf("\n");
	}
*/
//	double **forcecopy = atom->f;
//	double **speedcopy = atom->v;
//	int** specialcopy = atom->bond_type;
//	bool ifoundsomething = false;
	printf("force on atom is: \n");
	for (int indexOfParticle = 0; indexOfParticle < nlocal; indexOfParticle++) { //was: ++indexOfParticle++
//	printf("my bond type is: %d \n",atom->bond_type[2][2]);

//	printf("currentatom is %i \n", indexOfParticle);

		if (atom->mask[indexOfParticle] & groupbit) {
			printf("x: %f, y: %f, z: %f \n", forcecopy[indexOfParticle][0], forcecopy[indexOfParticle][1], forcecopy[indexOfParticle][2]);
			if (nlocal - indexOfParticle >= 3) {
				double distancederi [3][3];
				double forcederi [3][3];
				double distances[3];
				double sumofdistances = 0;
				for(int i = 0; i < 3; i++) { //calculating stuff for atom i
					double* a = poscopy[indexOfParticle];
					double* b = poscopy[indexOfParticle+i+1];
					distances[i] = euclideandistance(a,b);
					if (distances[i] < 0) {
						distances[i] *= -1;
					}
					sumofdistances+=distances[i];
					for (int x = 0; x < 3; x++) { //calculating stuff for direction x_x
						distancederi[i][x]=poscopy[indexOfParticle][x]-poscopy[indexOfParticle+i+1][x];
						forcederi[i][x]=forcecopy[indexOfParticle][x]-forcecopy[indexOfParticle+i+1][x];
					}
				}
//				printf("\n");
				for(int x = 0; x < 3; x++) {
					deriforce[indexOfParticle][x] = 0;
//					printf("\n");
					for (int i = 0; i < 3; i++) {
						double myweighting = distances[i]/sumofdistances;
//						printf("my weighting is %f, the added derif is %f ", myweighting, forcederi[i][x]);
						deriforce[indexOfParticle][x] +=  myweighting * forcederi[i][x];
					}
				}
//			printf("distances: %f %f %f \n", distances[0], distances[1], distances[2]);
			printf("derif: ");
			for (int i = 0; i < 3; i++) {
				printf(" %f",deriforce[indexOfParticle][i]);
			}
			printf("\n");
			}
		}
	}

//	printf("\n");
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
