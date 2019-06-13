#pragma once
#include <iostream>
#include <complex>
class FT
{
private:
	
public:
	FT();
	void DiscreteFourierTransform(int** InputImage, int** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w);
	void DFT(double** pFreqReal, double** pFreqImag, int** InputImage, int h, int w, int u, int v);

	void InverseDiscreteFourierTransform(int** InputImage, int** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w);
	void InverseDFT(double ** InverseReal, double ** InverseImag, double ** pFreqReal, double ** pFreqImag, int h, int w, int x, int y);

	void FastFourierTransform(int ** InputImage, int ** OutputImage,std::complex<float>** CImage, int h, int w);
	void FFT(double** pFreqReal, double** pFreqImag, int** InputImage, int h, int w, int u, int v);

	void InverseFastFourierTransform( int** OutputImage, std::complex<float>** CImage, int h, int w);
	void InverseFFT(double ** InverseReal, double ** InverseImag, double ** pFreqReal, double ** pFreqImag, int h, int w, int x, int y);

	void LowpassFilter(int ** OutputImage, std::complex<float> ** CImage, int h, int w, float Dcutoff, float n);
	void HighpassFilter(int ** OutputImage, std::complex<float> ** CImage, int h, int w, float Dcutoff, float n);

private:

};



