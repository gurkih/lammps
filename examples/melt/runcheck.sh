echo "1 mpi process"
mpirun -np 1 ../../src/lmp_ompi_g++ -in in.melt | grep x_f_deriv
echo "2 mpi processes"
mpirun -np 2 ../../src/lmp_ompi_g++ -in in.melt | grep x_f_deriv
echo "4 mpi processes"
mpirun -np 4 ../../src/lmp_ompi_g++ -in in.melt | grep x_f_deriv
echo "5 mpi processes"
mpirun -np 5 ../../src/lmp_ompi_g++ -in in.melt | grep x_f_deriv
echo "10 mpi processes"
mpirun -np 10 ../../src/lmp_ompi_g++ -in in.melt | grep x_f_deriv
echo "20 mpi processes"
mpirun -np 20 ../../src/lmp_ompi_g++ -in in.melt | grep x_f_deriv
echo "50 mpi processes"
mpirun -np 50 ../../src/lmp_ompi_g++ -in in.melt | grep x_f_deriv
