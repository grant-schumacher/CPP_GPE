#pragma once


//Library file for larger functions

//standard includes so VS doesn't complain

//#include "pch.h"
#include <iostream>
#include <utility>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

const double pi = 3.14159265358979323846;



//basic class for complex numbers and arithmetic
class Complex {
public:
	double real_part;
	double imaginary_part;


	//constructors
	Complex() {
		real_part = 0;
		imaginary_part = 0;
	}

	Complex(double x, double y) {
		real_part = x;
		imaginary_part = y;
	}

	Complex(double x) : Complex(x, 0) {};

	Complex(int x) : Complex((double)x) {};


	//Arithmetic, Comparison, stream overloads
	friend std::ostream& operator<<(std::ostream& os, const Complex& z)
	{
		os << "(" << z.real_part << ", " << z.imaginary_part << ")";
		return os;
	}


	void operator=(const Complex& rhs) {
		this->real_part = rhs.real_part;
		this->imaginary_part = rhs.imaginary_part;
	}

	void operator=(double rhs) {
		this->real_part = rhs;
		this->imaginary_part = 0;
	}

	friend Complex operator*(const Complex& lhs, const Complex& rhs) {
		double new_real_part = lhs.real_part * rhs.real_part - lhs.imaginary_part * rhs.imaginary_part;
		double new_imaginary_part = lhs.imaginary_part * rhs.real_part + lhs.real_part * rhs.imaginary_part;
		return Complex(new_real_part, new_imaginary_part);
	}

	friend Complex operator*(const double& lhs, const Complex& rhs) {
		double new_real_part = lhs * rhs.real_part;
		double new_imaginary_part = lhs * rhs.imaginary_part;
		return Complex(new_real_part, new_imaginary_part);
	}

	friend Complex operator*(const Complex& lhs, const double& rhs) {
		double new_real_part = lhs.real_part * rhs;
		double new_imaginary_part = lhs.imaginary_part * rhs;
		return Complex(new_real_part, new_imaginary_part);
	}

	friend Complex operator+(const Complex& lhs, const Complex& rhs) {
		double new_real_part = lhs.real_part + rhs.real_part;
		double new_imaginary_part = lhs.imaginary_part + rhs.imaginary_part;
		return Complex(new_real_part, new_imaginary_part);
	}

	friend Complex operator+(const double& lhs, const Complex& rhs) {
		double new_real_part = lhs + rhs.real_part;
		double new_imaginary_part = rhs.imaginary_part;
		return Complex(new_real_part, new_imaginary_part);
	}

	friend Complex operator+(const Complex& lhs, const double& rhs) {
		double new_real_part = lhs.real_part + rhs;
		double new_imaginary_part = lhs.imaginary_part;
		return Complex(new_real_part, new_imaginary_part);
	}

	friend Complex operator-(const Complex& lhs, const Complex& rhs) {
		double new_real_part = lhs.real_part - rhs.real_part;
		double new_imaginary_part = lhs.imaginary_part - rhs.imaginary_part;
		return Complex(new_real_part, new_imaginary_part);
	}

	friend Complex operator-(const double& lhs, const Complex& rhs) {
		double new_real_part = lhs - rhs.real_part;
		double new_imaginary_part = (-1)*rhs.imaginary_part;
		return Complex(new_real_part, new_imaginary_part);
	}

	friend Complex operator-(const Complex& lhs, const double& rhs) {
		double new_real_part = lhs.real_part - rhs;
		double new_imaginary_part = lhs.imaginary_part;
		return Complex(new_real_part, new_imaginary_part);
	}

	friend bool operator==(const Complex& lhs, const Complex& rhs) {
		return (lhs.real_part == rhs.real_part && lhs.imaginary_part == rhs.imaginary_part);
	}

	void operator+=(const Complex& rhs) {
		this->real_part = this->real_part + rhs.real_part;
		this->imaginary_part = this->imaginary_part + rhs.imaginary_part;
	}

	//swap overload for Complex (to allow for Cooley Tukey ordering for FFT)
	friend void swap(Complex& a, Complex& b)
	{
		using std::swap; // bring in swap for built-in types

		swap(a.real_part, b.real_part);
		swap(a.imaginary_part, b.imaginary_part);
	}

	//More math
	double Abs() {
		return (sqrt(pow(real_part, 2) + pow(imaginary_part, 2)));
	}

	double Abs_squared() {
		return (pow(real_part, 2) + pow(imaginary_part, 2));
	}

	Complex Exp() {
		return exp(real_part)*Complex(cos(imaginary_part), sin(imaginary_part));
	}

	Complex Exp_pure_imaginary() {
		return Complex(cos(imaginary_part), sin(imaginary_part));
	}

	double Arg() {
		return (atan2(imaginary_part, real_part));
	}
};

void Print2D(Complex* data, int n_x, int n_y) {
	for (int i = 0; i < n_x; i++) {
		for (int j = 0; j < n_y; j++) {
			std::cout << data[i*n_y + j];
			std::cout << '\t';
		}
		std::cout << '\n';
	}
}

//basic implementation of FFT, uses extra memory but puts the result back in the same memory (in order to simulate an in-place operation).
void FFT_naive(Complex* data, unsigned long n)\
{
	Complex* transform = new Complex[n];
	for (int k = 0; k < n; k++) {
		for (int j = 0; j < n; j++){
			transform[k] += (data[j] * (Complex(0, -2 * pi*k*j / n).Exp_pure_imaginary()));
		}
	}

	for (int k = 0; k < n; k++) {
		data[k] = transform[k];
	}

	delete[] transform;	
}


//FFT taken from Numerical recipes in C++, as quoted in "A Simple and Efficient FFT Implementation in C++", Vlodymyr Myrnyy, Dr. Dobbs Journal, 2007
//Note that this requires length to be a power of 2. 
//The initial signal is stored in the array data of length 2*nn, where each even element corresponds to the real part and each odd element to the imaginary part of a complex number.
//nn is limited to <= 2^14
void FFT_reference(double* data, unsigned long nn)
{
	using namespace std;
	unsigned long n, mmax, m, j, istep, i;
	double wtemp, wr, wpr, wpi, wi, theta;
	double tempr, tempi;

	// reverse-binary reindexing, to implement Cooley-Tukey
	n = nn << 1;
	j = 1;
	for (i = 1; i < n; i += 2) {
		if (j > i) {
			swap(data[j - 1], data[i - 1]);
			swap(data[j], data[i]);
		}
		m = nn;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	};

	// here begins the Danielson-Lanczos section
	mmax = 2;
	while (n > mmax) {
		istep = mmax << 1;
		theta = -(2 * 3.1415926535 / mmax);
		wtemp = sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi = sin(theta);
		wr = 1.0;
		wi = 0.0;
		for (m = 1; m < mmax; m += 2) {
			for (i = m; i <= n; i += istep) {
				j = i + mmax;
				tempr = wr * data[j - 1] - wi * data[j];
				tempi = wr * data[j] + wi * data[j - 1];

				data[j - 1] = data[i - 1] - tempr;
				data[j] = data[i] - tempi;
				data[i - 1] += tempr;
				data[i] += tempi;
			}
			wtemp = wr;
			wr += wr * wpr - wi * wpi;
			wi += wi * wpr + wtemp * wpi;
		}
		mmax = istep;
	}
}


//FFT designed to work with Complex number array of length nn.  This just unpacks the array to make FFT_reference work
void FFT_Complex(Complex* data_complex, unsigned long nn)
{
	double* data = new double[2 * nn];
	for (int index = 0; index < nn; index++) {
		data[2 * index] = data_complex[index].real_part;
		data[2 * index + 1] = data_complex[index].imaginary_part;
	}

	using namespace std;
	unsigned long n, mmax, m, j, istep, i;
	double wtemp, wr, wpr, wpi, wi, theta;
	double tempr, tempi;

	// reverse-binary reindexing, to implement Cooley-Tukey
	n = nn << 1;
	j = 1;
	for (i = 1; i < n; i += 2) {
		if (j > i) {
			swap(data[j - 1], data[i - 1]);
			swap(data[j], data[i]);
		}
		m = nn;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	};

	// here begins the Danielson-Lanczos section
	mmax = 2;
	while (n > mmax) {
		istep = mmax << 1;
		theta = -(2 * 3.1415926535 / mmax);
		wtemp = sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi = sin(theta);
		wr = 1.0;
		wi = 0.0;
		for (m = 1; m < mmax; m += 2) {
			for (i = m; i <= n; i += istep) {
				j = i + mmax;
				tempr = wr * data[j - 1] - wi * data[j];
				tempi = wr * data[j] + wi * data[j - 1];

				data[j - 1] = data[i - 1] - tempr;
				data[j] = data[i] - tempi;
				data[i - 1] += tempr;
				data[i] += tempi;
			}
			wtemp = wr;
			wr += wr * wpr - wi * wpi;
			wi += wi * wpr + wtemp * wpi;
		}
		mmax = istep;
	}

	for (int index = 0; index < nn; index++) {
		data_complex[index].real_part = data[2 * index];
		data_complex[index].imaginary_part = data[2 * index + 1];
	}
}

//Inverse Fourier Transform for 1 Dimension
void IFT_1D(Complex* data, const int n_x) {
	Complex* transform = new Complex[n_x];
	for (int j = 0; j < n_x; j++) {
		for (int x = 0; x < n_x; x++) {
				transform[j] += 1.0 / (n_x)*(data[x] * (Complex(0, 2 * (pi*j*x / n_x )).Exp_pure_imaginary()));
		}
	}


	for (int x = 0; x < n_x; x++) {
		data[x] = transform[x];
	}

	delete[] transform;
}

//Naive FT for 2D
void FFT_2D(Complex* data, const unsigned long n_x, const unsigned long n_y)\
{
	Complex* transform = new Complex[n_x*n_y];
	for (int j = 0; j < n_x; j++) {
		for (int k = 0; k < n_y; k++) {
			for (int x = 0; x < n_x; x++) {
				for (int y = 0; y < n_y; y++) {
					transform[j*n_y + k] += (data[x*n_y + y] * (Complex(0, -2 *(pi*j*x / n_x + pi*k * y / n_y)).Exp_pure_imaginary()));
				}
			}
		}
	}


	for (int x = 0; x < n_x; x++) {
		for (int y = 0; y < n_y; y++) {
			data[x*n_y + y] = transform[x*n_y + y];
		}
	}

	delete[] transform;
}

void IFT_2D(Complex* data, const int n_x, const int n_y) {
	Complex* transform = new Complex[n_x*n_y];
	for (int j = 0; j < n_x; j++) {
		for (int k = 0; k < n_y; k++) {
			for (int x = 0; x < n_x; x++) {
				for (int y = 0; y < n_y; y++) {
					transform[j*n_y + k] += 1.0/(n_x*n_y)*(data[x*n_y + y] * (Complex(0, 2 * (pi*j*x / n_x + pi * k * y / n_y)).Exp_pure_imaginary()));
				}
			}
		}
	}


	for (int x = 0; x < n_x; x++) {
		for (int y = 0; y < n_y; y++) {
			data[x*n_y + y] = transform[x*n_y + y];
		}
	}

	delete[] transform;
}

std::string WritePsiToFile(Complex* data, const int nx, const int ny, std::string filename) {
	std::ofstream file;
	file.open(filename);
	if (!file.is_open()) {
		return "Error Opening File";
	}

	for (int i = 0; i < nx; i++){

		file << data[i*ny + 0];
		for (int j = 1; j < ny; j++) {
			file << ',';
			file << data[i*ny + j];
		}
		file << '\n';
	}
	file.close();
	return "Success";
}

std::string WritePsi_Abs_SquaredToFile(Complex* data, const int nx, const int ny, std::string filename) {
	std::ofstream file;
	file.open(filename);
	if (!file.is_open()) {
		return "Error Opening File";
	}

	for (int i = 0; i < nx; i++) {
		file << (data[i*ny + 0].Abs_squared());

		for (int j = 1; j < ny; j++) {
			file << ',';
			file << (data[i*ny + j].Abs_squared());
		}
		file << '\n';
	}
	file.close();
	return "Success";
}

std::string WritePsi_ArgToFile(Complex* data, const int nx, const int ny, std::string filename) {
	std::ofstream file;
	file.open(filename);
	if (!file.is_open()) {
		return "Error Opening File";
	}

	for (int i = 0; i < nx; i++) {

		file << (data[i*ny + 0].Arg());
		for (int j = 1; j < ny; j++) {
			file << ',';
			file << (data[i*ny + j].Arg());
		}
		file << '\n';
	}
	file.close();
	return "Success";
}

namespace D_1
{
std::string WritePsiToFile(Complex* data, const int nx, std::string filename) {
	std::ofstream file;
	file.open(filename);
	if (!file.is_open()) {
		return "Error Opening File";
	}

	file << data[0];
	for (int j = 1; j < nx; j++) {
		file << ',';
		file << data[ j];
	}
	file << '\n';
	file.close();
	return "Success";
}

std::string WritePsi_Abs_SquaredToFile(Complex* data, const int nx, std::string filename) {
	std::ofstream file;
	file.open(filename);
	if (!file.is_open()) {
		return "Error Opening File";
	}

	
	file << (data[0].Abs_squared());
	for (int j = 1; j < nx; j++) {
		file << ',';
		file << (data[j].Abs_squared());
	}
	file << '\n';

	file.close();
	return "Success";
}

std::string WritePsi_ArgToFile(Complex* data, const int nx, std::string filename) {
	std::ofstream file;
	file.open(filename);
	if (!file.is_open()) {
		return "Error Opening File";
	}


	file << (data[0].Arg());
	for (int j = 1; j < nx; j++) {
		file << ',';
		file << (data[j].Arg());
	}
	file << '\n';

	file.close();
	return "Success";
}
}