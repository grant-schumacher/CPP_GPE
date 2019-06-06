#include "stdafx.h"
#include "CppUnitTest.h"
#include <locale>
#include <codecvt>
#include <string>

#include "../CPP_GPE/GPE_Tools.h"

//In order to allow Assert::AreEqual to work with Complex class
namespace Microsoft {
	namespace VisualStudio{
		namespace CppUnitTestFramework {
			template<> static std::wstring ToString<Complex>(const Complex& z) {
				std::string str = "(";
				str += std::to_string(z.real_part);
				str += ", ";
				str += std::to_string(z.imaginary_part);
				str += ")";
				std::wstring_convert<std::codecvt_utf8_utf16<wchar_t>> converter;
				return converter.from_bytes(str);
			}
		}
	} 
}

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace CPP_GPE_test
{		
	TEST_CLASS(FourierTransformTests)
	{
	public:
		/*
		TEST_METHOD(FT_ref_of_zero_size2)
		{
			const int size = 2;
			double data[2*size];
			double initial_data[2 * size];
			for (int i = 0; i < 2*size; i++) {
				initial_data[i] = data[i] = 0;
			}

			FFT_reference(data, size);

			for (int i = 0; i < 2 * size; i++) {
				Assert::AreEqual(initial_data[i], data[i]);
			}
		}

		TEST_METHOD(FT_ref_of_zero_size64)
		{
			const int size = 64;
			double data[2 * size];
			double initial_data[2 * size];
			for (int i = 0; i < 2 * size; i++) {
				initial_data[i] = data[i] = 0;
			}

			FFT_reference(data, size);

			for (int i = 0; i < 2 * size; i++) {
				Assert::AreEqual(initial_data[i], data[i]);
			}
		}

		TEST_METHOD(FT_ref_of_zero_size16384)
		{
			const int size = 16384;
			double data[2 * size];
			double initial_data[2 * size];
			for (int i = 0; i < 2 * size; i++) {
				initial_data[i] = data[i] = 0;
			}

			FFT_reference(data, size);

			for (int i = 0; i < 2 * size; i++) {
				Assert::AreEqual(initial_data[i], data[i]);
			}
		}
		*/
		TEST_METHOD(FT_size256)
		{
			const int size = 256;
			Complex data[size];
			Complex initial_data[size];
			for (int i = 0; i < size; i++) {
				initial_data[i] = 0.0;
				data[i] = 0.0;
			}

			FFT_Complex(data, size);

			for (int i = 0; i < size; i++) {
				Assert::AreEqual(initial_data[i], data[i]);
			}
		}

		TEST_METHOD(FT_naive_size256)
		{
			const int size = 256;
			Complex data[size];
			Complex initial_data[size];
			for (int i = 0; i < size; i++) {
				initial_data[i] = 0.0;
				data[i] = 0.0;
			}

			FFT_naive(data, size);

			for (int i = 0; i < size; i++) {
				Assert::AreEqual(initial_data[i], data[i]);
			}
		}

		TEST_METHOD(FT_naive_vs_FT_size256)
		{
			const int size = 256;
			Complex data[size];
			Complex data2[size];
			for (int i = 0; i < size; i++) {
				data2[i] = i+0.0;
				data[i] = i+0.0;
			}

			FFT_naive(data, size);
			FFT_Complex(data2, size);

			for (int i = 0; i < size; i++) {
				Assert::IsTrue((data[i]-data2[i]).Abs() < 0.0001, L"Disagreement between FFT methods.\n");
			}
		}

		TEST_METHOD(FT_8)
		{
			Complex data[8] = { 0, 0, 2, 3, 4, 0, 0, 0 };
			FFT_naive(data,8);
			Complex answer[8] = { sqrt(8)* Complex(3.18,0),
						 sqrt(8)*Complex(-2.16, -1.46),
						 sqrt(8)*Complex(0.71,1.06),
						 sqrt(8)*Complex(-0.66,-0.04),
						 sqrt(8)*Complex(1.06, 0.0),
						 sqrt(8)*Complex(-0.66, 0.04),
						 sqrt(8)*Complex(0.71, -1.06),
						 sqrt(8)*Complex(-2.16, 1.46) };
			for (int i = 0; i < 8; i++) {
				
				Assert::IsTrue(((data[i] - answer[i]).Abs()) < 0.1);
			}
		}

		TEST_METHOD(FT_2D_4x4)
		{
			Complex data[16];
			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < 4; j++) {
					data[i * 4 + j] = (Complex(0, 2 * pi*i / 4).Exp_pure_imaginary());
				}
			}
			FFT_2D(data, 4, 4);
			Complex answer[16] = {0,0,0,0,
								  16,0,0,0,
								  0,0,0,0,
								  0,0,0,0};
			for (int i = 0; i < 16; i++) {

				Assert::IsTrue(((data[i] - answer[i]).Abs()) < 0.0001);
			}
		}

		TEST_METHOD(FT_2D_16x16)
		{
			Complex data[256];
			for (int i = 0; i < 16; i++) {
				for (int j = 0; j < 16; j++) {
					data[i * 16 + j] = (Complex(0, 2 * (pi*i / 16 +3*pi*j/16)).Exp_pure_imaginary());
				}
			}
			FFT_2D(data, 16, 16);
			Complex answer[256];
			for (int i = 0; i < 16; i++) {
				for (int j = 0; j < 16; j++) {
					answer[i * 16 + j] = 0;
				}
			}
			answer[1 * 16 + 3] = 256;

			for (int i = 0; i < 256; i++) {

				Assert::IsTrue(((data[i] - answer[i]).Abs()) < 0.0001);
			}
		}

		TEST_METHOD(FT_2D_64x64)
		{
			Complex data[64*64];
			for (int i = 0; i < 64; i++) {
				for (int j = 0; j < 64; j++) {
					data[i * 64 + j] = (Complex(0, 2 * (pi*4*i / 64 + 8 * pi*j / 64)).Exp_pure_imaginary());
				}
			}
			FFT_2D(data, 64, 64);
			Complex answer[64*64];
			for (int i = 0; i < 64; i++) {
				for (int j = 0; j < 64; j++) {
					answer[i * 64 + j] = 0;
				}
			}
			answer[4 * 64 + 8] = 64*64;

			for (int i = 0; i < 64*64; i++) {

				Assert::IsTrue(((data[i] - answer[i]).Abs()) < 0.0001);
			}
		}

		TEST_METHOD(FT_IFT_2D_16x16)
		{
			Complex data1[16*16];
			Complex data2[16 * 16];
			for (int i = 0; i < 16*16; i++) {
				data1[i] = i;
				data2[i] = i;
			}
			FFT_2D(data1, 16, 16);
			IFT_2D(data1, 16, 16);

			for (int i = 0; i < 16*16; i++) {

				Assert::IsTrue(((data1[i] - data2[i]).Abs()) < 0.0001);
			}
		}

		TEST_METHOD(FT_IFT_1D_16)
		{
			Complex data1[16 ];
			Complex data2[16 ];
			for (int i = 0; i < 16 ; i++) {
				data1[i] = i;
				data2[i] = i;
			}
			FFT_naive(data1, 16);
			IFT_1D(data1, 16);

			for (int i = 0; i < 16 ; i++) {

				Assert::IsTrue(((data1[i] - data2[i]).Abs()) < 0.0001);
			}
		}


	};

	TEST_CLASS(ComplexTests) 
	{
	public:
		TEST_METHOD(MultiplyRealBy2)
		{
			Complex x(1, 0);
			Complex y(2, 0);			

			x = x * y;

			Assert::AreEqual<Complex>(x, Complex(2,0));
		}

		TEST_METHOD(MultiplyImaginaryBy2)
		{
			Complex x(0, 1);
			Complex y(2, 0);

			x = x * y;

			Assert::AreEqual<Complex>(x, Complex(0, 2));
		}

		TEST_METHOD(MultiplyRealBy2i)
		{
			Complex x(1, 0);
			Complex y(0, 2);

			x = x * y;

			Assert::AreEqual<Complex>(x, Complex(0, 2));
		}

		TEST_METHOD(MultiplyImaginaryBy2i)
		{
			Complex x(0, 1);
			Complex y(0, 2);

			x = x * y;

			Assert::AreEqual<Complex>(x, Complex(-2, 0));
		}

		TEST_METHOD(Magnitude) {
			Complex z(1, 1);
			Assert::AreEqual(z.Abs(), sqrt(2));
		}

		TEST_METHOD(EXPof1) {
			Complex z(1, 0);
			Assert::AreEqual<Complex>(z.Exp(), Complex(exp(1),0));
		}

		TEST_METHOD(EXPofI) {
			Complex z(0, 1);
			Assert::AreEqual<Complex>(z.Exp(), Complex(cos(1), sin(1)));
		}

		TEST_METHOD(EXP_imaginary_ofI) {
			Complex z(0, 1);
			Assert::AreEqual<Complex>(z.Exp_pure_imaginary(), Complex(cos(1), sin(1)));
		}

		TEST_METHOD(EXP) {
			Complex z(2, 3);
			Assert::AreEqual<Complex>(z.Exp(), exp(2)*Complex(cos(3), sin(3)));
		}

	};
}