#ifdef FIX_CLASS
	FixStyle(force/deri,ForceDerivative)
#else

#ifndef LMP_FIX_DERIFORCE_H
#define LMP_FIX_DERIFORCE_H

#include "fix.h"

namespace LAMMPS_NS {

class ForceDerivative : public Fix {
	public:
		ForceDerivative(class LAMMPS *, int, char **);
		~ForceDerivative();
		int setmask();
		void end_of_step();
		void sort(int);
		double euclideandistance(double *, double *);
		double memory_usage();
		void grow_arrays(int);
		void copy_arrays(int, int);
		void set_arrays(int);

	private:
		double **lastf;
		double **poscopy;
		double **forcecopy;



};

}

#endif
#endif
