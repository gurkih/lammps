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

void ForceDerivative::end_of_step() {
	int nlocal = atom->nlocal;
	double deriforceoutput[nlocal][3];
	double forceoutput[nlocal][3];
	double speedoutput[nlocal][3];

	float averagedenominator[3] = {0,0,0};

	float averagederif[3] = {0,0,0};
	printf("i am being run!");


	//here: calculate the DeriForce for each atom. Note that this will just calculate the DeriForce for the atoms this thread owns (like stated in atom->mask)

	for (int indexOfParticle = 0; indexOfParticle < nlocal; ++indexOfParticle) {	

	double **forcecopy = atom->f;
	double **speedcopy = atom->v;
	int** specialcopy = atom->bond_type;
	printf("i managed to copy!\n");
//	printf("my bond type is: %d \n",atom->bond_type[2][2]);
		if (atom->mask[indexOfParticle] & groupbit) {
		//	printf("Hello World");
			if (this->lastf != 0) {
				//printf("nspecial = %d \n",specialcopy[0][0]);
				
				//printf(" %d \n,",nlocal);
				for(int i = 0; i < 3; i++) {

					// Note: forceoutput and speedoutput are not needed for the deriforce computations. They are still being calculated here in order to save a for-loop
					forceoutput[indexOfParticle][i] = forcecopy[indexOfParticle][i]; //this calcs the f for each atom
					speedoutput[indexOfParticle][i] = speedcopy[indexOfParticle][i]; //this calcs the v for each atom

					deriforceoutput[indexOfParticle][i] = lastf[indexOfParticle][i]-forcecopy[indexOfParticle][i]; //this calcs the derivative force
					if (deriforceoutput[indexOfParticle][i] != deriforceoutput[indexOfParticle][i]) { //this checks for nan and overwrites with 0 if nan
						deriforceoutput[indexOfParticle][i] = 0; 
						//this might be nonsense; TODO: check for plausibility
					} else {
						deriforceoutput[indexOfParticle][i] = deriforceoutput[indexOfParticle][i]*deriforceoutput[indexOfParticle][i];
					}
				}	
			} 
		}
	}
	
	for (int indexOfParticle = 0; indexOfParticle < nlocal; ++indexOfParticle) {	
		if (atom->mask[indexOfParticle] & groupbit) {
			//printf("%f is gonna be added \n", averagederif[0]); //deriforceoutput[indexOfParticle][0]);	
			

	//now: calculate -F_i / v_i
		for(int i = 0; i < 3; i++) {
			if(speedoutput[indexOfParticle][i] != 0) {
				averagedenominator[i] += -forceoutput[indexOfParticle][i] / speedoutput[indexOfParticle][i];
			}
		
	//now: calculate the average DeriForce over threads atoms. Note that - again - we have to just use our threads atoms - so compare with atom->mask
			averagederif[i] = averagederif[i] + (float)deriforceoutput[indexOfParticle][i];	
			//printf("%f is gonna be added \n", averagederif[0]);	

			//averagederif[i] /= (float)nlocal;
			//averagedenominator[i] /= (float)nlocal;

			printf("averagederif is %f , averagedenominator = %f ",averagederif[i], averagedenominator[i]);
		}
		printf("\n");
		}
	}

	for (int i = 0; i < 3; i++) {
		averagederif[i] /= (float)nlocal;
		averagedenominator[i] /= (float)nlocal;
	}


	float globalderif[3];
	float globaldenominator[3];
	// mpi stuff is taken from http://mpitutorial.com/tutorials/mpi-reduce-and-allreduce/	
	//int num_elements_per_proc = atoi(argv[1]);
	int num_elements_per_proc = 1;
	
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	MPI_Reduce(&averagederif[0], &globalderif[0], 3, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&averagedenominator[0], &globaldenominator[0], 3, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
	
	for (int i = 0; i < 3; i++) {

		globalderif[i] /= (world_size * num_elements_per_proc);
		globaldenominator[i] /= (world_size * num_elements_per_proc);

	}
	//printf("worldsize: %d \n",world_size);

	if(world_rank == 0)
	{
	float sumderif = 0;
	float sumdenom = 0;
	for (int i = 0; i < 3; i++) {
		sumderif += globalderif[i];
		sumdenom += globaldenominator[i];
	}
	float result = sumderif/(sumdenom*3);
	printf("Temperature: %f \n", result);
	
	//	printf("x_f_deriv = %f, y_f_deriv = %f, z_f_deriv = %f, globaldenom = %f,  nlocal = %d, worldsize = %d \n", globalderif[0], globalderif[1], globalderif[2], globaldenominator[0], nlocal, world_size);
	}

	//MPI_Reduce TODO

	double **forcecopy = atom->f;
	for (int indexOfParticle = 0; indexOfParticle < nlocal; ++indexOfParticle) 
	{
		for (int i = 0; i < 3; i++) {
			this->lastf[indexOfParticle][i] = forcecopy[indexOfParticle] [i];
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
