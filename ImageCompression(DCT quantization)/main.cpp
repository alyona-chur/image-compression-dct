#include <opencv2/opencv.hpp>
#include <opencv/highgui.h>
#include <cstdlib>
#include <vector>
#include <iostream>
#include "dct.h"
#include "mtrx_functions.h"
#include "quantization.h"

using namespace std;
using namespace cv; //Open CV 2.4.13.3
typedef long double ld;
typedef vector < vector <ld> > vvld;

void getChanels(const Mat &img, vvld &chanel0, vvld &chanel1, vvld &chanel2) {
	int rows = img.rows;
	chanel0.resize(rows);
	chanel1.resize(rows);
	chanel2.resize(rows);
	for (int i = 0; i < img.rows; ++i) {
		for (int j = 0; j < img.cols; ++j) {
			chanel0[i].push_back((ld)img.at<Vec3b>(i, j)[0]);
			chanel1[i].push_back((ld)img.at<Vec3b>(i, j)[1]);
			chanel2[i].push_back((ld)img.at<Vec3b>(i, j)[2]);
		}
	}
}

void setChanels(Mat &img, const vvld &chanel0, const vvld &chanel1, const vvld &chanel2) {
	for (int i = 0; i < img.rows; ++i) {
		for (int j = 0; j < img.cols; ++j) {
			img.at<Vec3b>(i, j)[0] = (unsigned char)chanel0[i][j];
			img.at<Vec3b>(i, j)[1] = (unsigned char)chanel1[i][j];
			img.at<Vec3b>(i, j)[2] = (unsigned char)chanel2[i][j];
		}
	}
}

int main() {
	//Read
	long double q;
	cout << "Enter quantization quality (0 - no changes)" << endl;
	cin >> q;
	Mat imgIn = imread("test.jpg");
	int cols = imgIn.cols, rows = imgIn.rows;
	//Show input img
	namedWindow("in", WINDOW_AUTOSIZE);
	imshow("in", imgIn);

	//Get chanels from img
	vvld chanel0, chanel1, chanel2;
	getChanels(imgIn, chanel0, chanel1, chanel2); 
	//Do DCT and quantization
	vvld chanel0DCT = chanel0, chanel1DCT = chanel1, chanel2DCT = chanel2;
	DCT(chanel0, chanel0DCT);
	DCT(chanel1, chanel1DCT);
	DCT(chanel2, chanel2DCT); 
	quant(chanel0DCT, chanel0, q);
	quant(chanel1DCT, chanel1, q);
	quant(chanel2DCT, chanel2, q);
	//Do it back
	quant_back(chanel0, chanel0DCT, q);
	quant_back(chanel1, chanel1DCT, q);
	quant_back(chanel2, chanel2DCT, q);
	DCT_back(chanel0DCT, chanel0);
	DCT_back(chanel1DCT, chanel1);
	DCT_back(chanel2DCT, chanel2); 
	
	//Show undo img
	Mat imgOut = imgIn;
	setChanels(imgOut, chanel0, chanel1, chanel2);
	namedWindow("out", WINDOW_AUTOSIZE);
	imshow("out", imgOut);
	waitKey(0);
}