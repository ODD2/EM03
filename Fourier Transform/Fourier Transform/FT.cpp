#include "FT.h"
#include <complex>
#include <vector>
#define �k 3.1415926
using namespace std;


void CalcFFT(vector<complex<float>>& x, int N, bool inv = false) {
	/* bit-reversal permutation */
	

	/* dynamic programming */
	//for (int k = 2; k <= N; k <<= 1)
	//{
	//	float �c = (inv?-1:1)* -2.0 * �k / k;
	//	complex<float> d�s(cos(�c), sin(�c));
	//	// �Ck�Ӱ��@��FFT
	//	for (int j = 0; j < N; j += k)
	//	{
	//		// �ek/2�ӻP��k/2���T����ƭȫ�n��١A
	//		// �]������٪��@�_���C
	//		complex<float> �s(1, 0);
	//		for (int i = j; i < j + k / 2; i++)
	//		{
	//			if (i < N)
	//			{
	//				complex<float> a = target[i];
	//				complex<float> b;
	//				if (i + k / 2 < N) {
	//					b = target[i + k / 2] * �s;
	//					target[i + k / 2] = a - b;
	//				}
	//				else {
	//					b = 0;
	//				}
	//				target[i] = a + b;
	//			}
	//			�s *= d�s;
	//		}
	//	}
	//}


	for (int k = N; k >= 2; k >>= 1)
	{
		float �c = (inv ? -1 : 1)* 2.0 * �k / k;
		complex<float> d�s(cos(�c), sin(�c));
		complex<float> �s(1, 0);
		for (int i = 0; i < k / 2; i++)
		{
			for (int j = i; j < N; j += k)
			{
				// �s�ܦ��۴��~���W�C
				complex<float> a = x[j] ;
				complex<float> b;
				if (j + k / 2 < N)
				{
					b = x[j + k / 2];
					x[j + k / 2] = (a - b)* �s;
				}
				else {
					b = complex<float>(0,0);
				}
				x[j] = a + b;
			}
			�s *= d�s;
		}
	}
	
	for (int i = 1, j = 0; i < N; ++i)
	{
		for (int k = N >> 1; !((j ^= k)&k); k >>= 1);
		//      for (int k=N>>1; k>(j^=k); k>>=1) ;
		if (i > j) swap(x[i], x[j]);
		//      if (i<j) swap(target[i], target[j]);
	}
}

FT::FT()
{
}

void FT::DiscreteFourierTransform(int** InputImage, int** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w)
{
	int M = h;
	int N = w;

	double** pFreq = new double*[M];

	for (int newcnt = 0; newcnt<M; newcnt++)
	{
		pFreq[newcnt] = new double[N]; // �ť߸��W�v�}�C
	}

	for (int forzero_i = 0; forzero_i<M; forzero_i++)
	{
		for (int forzero_j = 0; forzero_j<N; forzero_j++)
		{
			pFreq[forzero_i][forzero_j] = 0.0f;
		}
	}
	//-------------------------------------------
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			DFT(FreqReal, FreqImag, InputImage,M, N, j, i);
		}
	}

	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			// �N�p��n���ť߸���ƻP��Ƴ����@���X 
			pFreq[i][j] = sqrt(pow(FreqReal[i][j], (double) 2.0) + pow(FreqImag[i][j], (double) 2.0));
			// ���X�ᤧ�W�v���J�v���}�C����� 
			OutputImage[i][j] = pFreq[i][j];
		}
	}

	//-------------------------------------------
	for (int delcnt = 0; delcnt < M; delcnt++)
	{
		delete[] pFreq[delcnt];
	}
	delete[] pFreq;
}

void FT::DFT(double ** pFreqReal, double ** pFreqImag, int ** InputImage, int h, int w, int u, int v)
{
	// M = N �����O��}
	int M = h;
	int N = w;

	for (int y = 0; y < M; y++)
	{
		for (int x = 0; x < N; x++)
		{
			// �i���p��Eular's equation e^{j*theta} = cos{theta}+j*sin{theta}			
			double angleDFT = (-1.0f * 2.0f * 3.14159 * (double)(u*x + v*y) / (double)M);
			double c = cos(angleDFT);
			double s = sin(angleDFT);

			// �Q��Eular's equation�p��ť߸������Ƴ���
			pFreqReal[u][v] += (double)InputImage[y][x] * c;
			pFreqImag[u][v] -= (double)InputImage[y][x] * s;
		}
	}

	pFreqReal[u][v] = pFreqReal[u][v] / (double)(M);
	pFreqImag[u][v] = pFreqImag[u][v] / (double)(M);
}

void FT::InverseDiscreteFourierTransform(int ** InputImage, int ** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w)
{
	int M = h;
	int N = w;

	double** InverseReal = new double*[M];
	double** InverseImag = new double*[M];
	double** pFreq = new double*[M];

	for (int i = 0; i<M; i++)
	{
		InverseReal[i] = new double[N];
		InverseImag[i] = new double[N];
		pFreq[i] = new double[N]; // �ť߸��W�v�}�C
	}

	for (int i = 0; i<M; i++)
	{
		for (int j = 0; j<N; j++)
		{
			InverseReal[i][j] = 0.0f;
			InverseImag[i][j] = 0.0f;
			pFreq[i][j] = 0.0f;
		}
	}
	//-------------------------------------------
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			InverseDFT(InverseReal, InverseImag,FreqReal, FreqImag, M, N, j, i);
		}
	}
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			// �N�p��n���ť߸���ƻP��Ƴ����@���X 
			pFreq[i][j] = sqrt(pow(InverseReal[i][j], (double) 2.0) + pow(InverseImag[i][j], (double) 2.0));
			// ���X�ᤧ�W�v���J�v���}�C����� 
			OutputImage[i][j] = pFreq[i][j];
			//�s�U�ϳť߸���ƻP��Ƴ���
			FreqReal[i][j] = InverseReal[i][j];
			FreqImag[i][j] = InverseImag[i][j];

		}
	}
	//-------------------------------------------
	for (int i = 0; i < M; i++)
	{
		delete[] pFreq[i];
		delete[] InverseReal[i];
		delete[] InverseImag[i];

	}
	delete[] pFreq;
	delete[] InverseReal;
	delete[] InverseImag;

}

void FT::InverseDFT(double ** InverseReal, double ** InverseImag, double ** pFreqReal, double ** pFreqImag, int h, int w, int x, int y)
{
	// M = N �����O��}
	int M = h;
	int N = w;

	for (int v = 0; v < M; v++)
	{
		for (int u = 0; u < N; u++)
		{
			// �i���p��Eular's equation e^{j*theta} = cos{theta}+j*sin{theta}			
			double angleIDFT = (2.0f * 3.14159 * (double)(u*x + v*y) / (double)M);
			double c = cos(angleIDFT);
			double s = sin(angleIDFT);

			// �Q��Eular's equation�p��ť߸������Ƴ���
			InverseReal[x][y] += (pFreqReal[v][u] * c - pFreqImag[v][u] * s);
			InverseImag[x][y] += (pFreqReal[v][u] * s + pFreqImag[v][u] * c);
		}
	}
	InverseReal[x][y] = InverseReal[x][y] / (float)M;
	InverseImag[x][y] = InverseImag[x][y] / (float)M;
}

void FT::FastFourierTransform(int ** InputImage, int ** OutputImage,complex<float>** CImage, int h, int w)
{
	vector<complex<float>> target;
	int M = h;
	int N = w;

	for (int i = 0; i < M; ++i) {
		for (int j = 0; j < N; ++j) {
			CImage[i][j] =complex<double>(pow(-1, i + j)*InputImage[i][j], 0);
		}
	}
	//-------------------------------------------
	target.resize(N);
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			target[j] = CImage[i][j];
		}
		CalcFFT(target,N);
		for (int j = 0; j < N; j++)
		{
			CImage[i][j] = target[j];
		}
	}

	target.resize(M);
	for (int j = 0; j < N; j++)
	{
		for (int i = 0; i < M; i++)
		{
			target[i] = CImage[i][j];
		}
		CalcFFT(target,M);
		for (int i = 0; i < M; i++)
		{
			CImage[i][j] = target[i];
		}
	}

	double scaler = pow(M * N,0.5);
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			// �N�p��n���ť߸���ƻP��Ƴ����@���X 
			// ���X�ᤧ�W�v���J�v���}�C����� 
			OutputImage[i][j] = pow(pow(CImage[i][j].real(), 2) + pow(CImage[i][j].imag() , 2), 0.5)/scaler;
		}
	}
	//-------------------------------------------
}

void FT::FFT(double** pFreqReal, double** pFreqImag, int** InputImage, int h, int w, int u, int v)
{

}

void FT::InverseFastFourierTransform(int** OutputImage, std::complex<float>** CImage, int h, int w)
{
	vector<complex<float>> target;
	int M = h;
	int N = w;
	//-------------------------------------------
	target.resize(N);
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			target[j] = CImage[i][j];
		}
		CalcFFT(target, N,true);
		for (int j = 0; j < N; j++)
		{
			CImage[i][j] = target[j];
		}
	}

	target.resize(M);
	for (int j = 0; j < N; j++)
	{
		for (int i = 0; i < M; i++)
		{
			target[i] = CImage[i][j];
		}
		CalcFFT(target, M,true);
		for (int i = 0; i < M; i++)
		{
			CImage[i][j] = target[i];
		}
	}

	float scaler = M * N;
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{	float scaler = pow(M*N, 0.5);
			CImage[i][j]._Val[0] /=scaler;
			CImage[i][j]._Val[1]/=scaler;
			// �N�p��n���ť߸���ƻP��Ƴ����@���X 
			// ���X�ᤧ�W�v���J�v���}�C����� 
			OutputImage[i][j] = pow(pow(CImage[i][j].real(), 2) + pow(CImage[i][j].imag(), 2), 0.5)/scaler;
		}
	}
}

void FT::InverseFFT(double ** InverseReal, double ** InverseImag, double ** pFreqReal, double ** pFreqImag, int h, int w, int x, int y)
{

}

void FT::LowpassFilter(int ** OutputImage,complex<float> ** CImage, int h, int w,float Dcutoff,float n)
{
	double scaler = h * w;
	double cw=w/2.0, ch=h/2.0, u, v;
	for (int i = 0; i < h; ++i) {
		for (int j = 0; j < w; ++j) {
			u = j - cw;
			v = i - ch;
			CImage[i][j] *= 1.0 /(1+ pow(pow(u*u + v* v, 0.5) / Dcutoff, 2 * n));
			OutputImage[i][j] = pow(pow(CImage[i][j].real(), 2) + pow(CImage[i][j].imag(), 2), 0.5);
		}
	}
}

void FT::HighpassFilter( int ** OutputImage, complex<float> ** CImage, int h, int w, float Dcutoff, float n)
{
	double scaler = h * w;
	double cw = w / 2.0, ch = h / 2.0, u, v;
	for (int i = 0; i < h; ++i) {
		for (int j = 0; j < w; ++j) {
			u = j - cw;
			v = i - ch;
			CImage[i][j] *= 1- 1.0 / (1 + pow(pow(u*u + v * v, 0.5) / Dcutoff, 2 * n));
			OutputImage[i][j] = pow(pow(CImage[i][j].real(), 2) + pow(CImage[i][j].imag(), 2), 0.5);
		}
	}
}


