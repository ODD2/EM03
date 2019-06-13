#pragma once
#include <complex>


class DataManager
{
private:
	int** InputImage; //��J�v��
	int** OutputImage; //��X�v��
	double ** FreqReal; // �ť߸���Ƴ���
	double ** FreqImag; // �ť߸���Ƴ���
	int ImageHeight;
	int ImageWidth;
	std::complex<float>** CImage=nullptr;
public:
	DataManager(int h, int w);
	~DataManager();
	//�]�w�W�v��T
	void SetPixel(int x, int y, int pixelValue);
	void SetFreqReal(int x, int y, double value);
	void SetFreqImag(int x, int y, double value);
	void DataManager::SetCImage(std::complex<float> ** c);
	int GetImageHeight();
	int GetImageWidth();

	//���o�ܼƦ�������
	int** GetInputImage();
	int** GetOutputImage();
	double** GetFreqReal();
	double** GetFreqImag();
	std::complex<float> ** GetCImage();
};


