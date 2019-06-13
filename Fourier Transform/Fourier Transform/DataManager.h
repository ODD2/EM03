#pragma once
#include <complex>


class DataManager
{
private:
	int** InputImage; //輸入影像
	int** OutputImage; //輸出影像
	double ** FreqReal; // 傅立葉實數部分
	double ** FreqImag; // 傅立葉虛數部分
	int ImageHeight;
	int ImageWidth;
	std::complex<float>** CImage=nullptr;
public:
	DataManager(int h, int w);
	~DataManager();
	//設定頻率資訊
	void SetPixel(int x, int y, int pixelValue);
	void SetFreqReal(int x, int y, double value);
	void SetFreqImag(int x, int y, double value);
	void DataManager::SetCImage(std::complex<float> ** c);
	int GetImageHeight();
	int GetImageWidth();

	//取得變數成員指標
	int** GetInputImage();
	int** GetOutputImage();
	double** GetFreqReal();
	double** GetFreqImag();
	std::complex<float> ** GetCImage();
};


