#include "pch.h"
#include <iostream>
#include <utility>
#include <exception>
#include "GPE_Tools.h"
#include <fstream>

int main()
{
	using namespace std;
	//Basic simulation of GPE


	//We can first "non-dimensionalize" the GPE, getting something like eq. 3.1 in Weizhu Bao, Dieter Jaksch, Peter Markowich "Numerical solution of the GPE for BEC":
	// i*e*PHI_T = -e^2/2 * PHI_XX + x^2/2 * PHI + K*|PHI|^2 * PHI
	// For a harmonic trap, with frequency w, mass m, and system size x. Where e = hbar/(w*m*x^2) and K = sgn(a)/2*(a0/xh * a0/x)^2, where a is the scattering length, a0 is the harmonic length sqrt(hbar/w*m), and xh is the healing length (8*Pi*|a|*N/x^3)^(-1/2).

	// We enforce periodic boundary conditions, and divide our region into units of size h = (a-b)/M = Length/M, for M an even number.  We pick time step k.  Let x_j = a + jh,t_n = nk.  Let Psi^n_j be the approximation of Psi(x_j,t_n), and Psi^n the vector of these points, Psi^n_j.
	// We define a time-splitting spectral method (TSSP) with the following:
	// solve:  i*e*dPsi/dt = -e^2/2 * d^2Psi/dx^2 for time step k, using spectral methods (eq. 3.4 loc. cit).  
	//Then solve:  i*e*dPsi/dt = x^2/2*Psi(x,t)+K*|Psi(x,t)|^2 * Psi(x,t) for the same time step (eq. 3.5 loc. cit).
	// It turns out that eq. 3.5 leaves |Psi| invariant, so we can replace |Psi|^2 in the above with the initial value for timestep n: i*e*dPsi/dt = x^2/2*Psi(x,t)+K*|Psi(x_n,t)|^2 * Psi(x,t), so this is now exactly integrable.
	// For the simulated step between timestep n and timestep n+1, we combine the two equations above with Strang Splitting.
	// Psi*_j = exp(-i(x_j^2/2 + K*|Psi^n_j|^2)*k/(2e))*Psi^n_j
	// Psi**_j = 1/M * SUM(l = -M/2, l = M/2-1, exp(-i*e*k*u_l^2/2) Psi-hat*_l exp(i*u_l*(x_j-a))
	// Psi^n+1_j = exp(-i*(x_j/2+K|Psi**_j|^2)*k/(2*e) * Psi**_j
	// where Psi-hat*_l, the fourier components of Psi*, are defined as: u_l = 2*Pi*l/(b-a)
	// and Psi-hat*_l = SUM(j = 0, j = M-1, Psi*_j exp(-i*u_l(x_j - a))) and l is in the range [-M/2, M/2 - 1]
	// Due to Strang splitting, the time discretization error is second order in k. The space-discretization error is machine limited (as the spectral method is only limited by the fourier transform precision).

	//Physical Constants
	double hbar = 1.0545718e-34; //Joule-seconds
	double mass_of_potassium39_u = 38.96370668; //in u
	double mass_of_potassium39_kg = 1.66054e-27*mass_of_potassium39_u; //in kg
	double bohr_radius = 5.2917721067e-11; //Bohr Radius, in meters
	double scattering_length_of_potassium39 = 64.61; //in units of Bohr Radii

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//Set parameters (height, width, mass, trap, time, etc.)
//in the future, replace this with a file read to keep from having to recompile to change parameters
	//physical parameters
		double height = 1; // height in meters
		double width = 1; //width in meters
		double mass = mass_of_potassium39_kg; //mass of atom in kg
		double trap_frequency = 100; //frequency of harmonic trap, in Hertz
		double N = 100000; //number of atoms in trap
		double time = 10; //length of simulation in seconds
		double characteristic_size = 0.001; //characteristic length of the condensate, in meters
		double epsilon = hbar/(mass*trap_frequency*characteristic_size); // a non-dimensionalized energy scale, equals hbar/(w*m*x^2)
		double a = scattering_length_of_potassium39 *bohr_radius; //scattering length, in meters
		double a0 = sqrt(hbar/(trap_frequency*mass)); //harmonic lengthscale, sqrt(hbar/w*m)
		double healing_length = pow(8 * pi*abs(a)*N / pow(characteristic_size, 3), -0.5);//the healing length, (8 * Pi* | a | *N / x ^ 3) ^ (-1 / 2)
		double K; //Interaction parameter, sgn(a) / 2 * (a0 / xh * a0 / x) ^ 2
		if (a > 0) { K = 0.5*pow(a0 / healing_length * a0 / characteristic_size, 2); }
		else K = -0.5*pow(a0 / healing_length * a0 / characteristic_size, 2);

	//simulation parameters
		int nx = 32; //number of space steps in x-axis, must be a power of 2;
		int ny = 32; //number of space steps in y-axis, must be a power of 2;
		int nt = 1000; //number of time steps

		double x_stepsize = height / nx; //physical size of single step in x-axis, in meters
		double y_stepsize = width / ny;//physical size of single step in y-axis, in meters
		double t_stepsize = time / nt;//physical duration of single step in time, in seconds
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//open files to store all data

/*		ofstream output_file;
		output_file.open("C:\\Users\\Grant\\Desktop\\CPPGPE\\TestOutput\\test1.txt");
		if (output_file.is_open()) { std::cout << "opening file test1.txt\n"; }
		else {
			std::cout << "Failed to open file.\n";
			return 0;
		}*/

	//Initialize an array of Complex Numbers to store the field;
		Complex* Psi = new Complex[nx*ny];

	//Start in some initial state

		//for now, just start in a gaussian distribution
		for (int x = 0; x < nx; x++) {
			for (int y = 0; y < ny; y++) {
				Psi[x*ny + y] = exp(-pow(x-nx/2, 2)/nx - pow(y-ny/2, 2)/ny);
			}
		}

		WritePsi_Abs_SquaredToFile(Psi, nx, ny, "C:\\Users\\Grant\\Desktop\\CPPGPE\\TestOutput\\step_initial.csv");

    //Main loop
		for (int t = 0; t < nt; t++) {
			//Strang Split
				//Solve for  Psi*_j as a function of Psi^n_j
			for (int x = 0; x < nx; x++) {
				for (int y = 0; y < ny; y++) {
					Psi[x*ny + y] = Psi[x*ny + y] *
						((Complex(0, -0.5 * t_stepsize / epsilon * (pow(x*x_stepsize, 2) / 2 + pow(y*y_stepsize, 2) / 2 + K * (Psi[x*ny + y].Abs_squared())))).Exp_pure_imaginary());
				}
			}

			//Solve for  Psi**_j as a function of Psi*_j
			FFT_2D(Psi, nx, ny);

			for (int x = 0; x < nx; x++) {
				for (int y = 0; y < ny; y++) {
					Psi[x*ny + y] = Psi[x*ny + y] *
						(Complex(0, -0.5*epsilon*t_stepsize*(pow(2.0*pi*x / height, 2) + pow(2.0*pi*y / width, 2))).Exp_pure_imaginary());
				}
			}
			IFT_2D(Psi, nx, ny);

			//Solve for Psi^n+1_j as a function of Psi**_j
			for (int x = 0; x < nx; x++) {
				for (int y = 0; y < ny; y++) {
					Psi[x*ny + y] = Psi[x*ny + y] *
						((Complex(0, -0.5 * t_stepsize / epsilon * (pow(x*x_stepsize, 2) / 2 + pow(y*y_stepsize, 2) / 2 + K * (Psi[x*ny + y].Abs_squared())))).Exp_pure_imaginary());
				}
			}
			//Write the current state to a file every so often
			WritePsi_Abs_SquaredToFile(Psi, nx, ny, "C:\\Users\\Grant\\Desktop\\CPPGPE\\TestOutput\\Abs\\step_"+to_string(t)+".csv");
			WritePsi_ArgToFile(Psi, nx, ny, "C:\\Users\\Grant\\Desktop\\CPPGPE\\TestOutput\\Arg\\step_" + to_string(t) + ".csv");

		}

	//Check for conserved quantities and write summary file

	//close files
	//output_file.close();
}


/*
//Strang Split Helper functions
void step_one(Complex* data, int nx, int ny, double t) {
	for (int x = 0; x < nx; x++) {
		for (int y = 0; y < ny; y++) {
			data[x*ny + y] = (Complex(0, -1.0 * t / epsilon * (pow(x*x_stepsize, 2) / 2 + K * (data[x*ny + y].Abs_squared()))));
		}
	}
}*/