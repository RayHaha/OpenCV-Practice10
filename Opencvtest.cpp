// Opencvtest.cpp : 定義主控台應用程式的進入點。
//

#include "stdafx.h"
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/core/core.hpp>
#include <iostream>
#include <fstream>
#include <string.h>
#include <time.h>
#include <math.h>

using namespace cv;
using namespace std;

void LaplacianMask1(Mat originimage, Mat resultimage, double threshhold);
void LaplacianMask2(Mat originimage, Mat resultimage, double threshold);
void MinimumVarianceLaplacian(Mat originimage, Mat resultimage, double threshold);
void LaplacianOfGaussian(Mat originimage, Mat resultimage, double threshold);
void DifferenceOfGaussian(Mat originimage, Mat resultimage, double threshold);

int _tmain(int argc, _TCHAR* argv[])
{
	Mat image = imread("lena.bmp");
	Mat grayimage;
	int pixel = 0;
	cvtColor(image,grayimage, CV_BGR2GRAY);
	
	Mat laplacian1image = grayimage.clone();
	LaplacianMask1(grayimage, laplacian1image, 32);
	imwrite("mask1.bmp",laplacian1image);

	Mat laplacian2image = grayimage.clone();
	LaplacianMask2(grayimage, laplacian2image, 25);
	imwrite("mask2.bmp",laplacian2image);

	Mat minimumimage = grayimage.clone();
	MinimumVarianceLaplacian(grayimage, minimumimage, 20);
	imwrite("minimum.bmp",minimumimage);

	Mat laplacianofgaussianimage = grayimage.clone();
	LaplacianOfGaussian(grayimage, laplacianofgaussianimage, 4000);
	imwrite("laplacian of gaussian.bmp",laplacianofgaussianimage);

	Mat differenceoggaussianimage = grayimage.clone();
	DifferenceOfGaussian(grayimage, differenceoggaussianimage, 1);
	imwrite("difference of gussian.bmp", differenceoggaussianimage);

	

	waitKey(0);
	return 0;
}

void LaplacianMask1(Mat originimage, Mat resultimage, double threshold){
	for(int i=1; i<originimage.rows-1; i++){
		for(int j=1; j<originimage.cols-1; j++){
			int pixel = originimage.at<uchar>(i,j);
			int left = originimage.at<uchar>(i,j-1);
			int right = originimage.at<uchar>(i,j+1);
			int top = originimage.at<uchar>(i-1,j);
			int down = originimage.at<uchar>(i+1,j);

			int temp = pixel*(-4)+left+right+top+down;
			if(temp>=threshold){
				resultimage.at<uchar>(i,j) = 0;
			}else{
				resultimage.at<uchar>(i,j) = 255;
			}
		}
	}
}

void LaplacianMask2(Mat originimage, Mat resultimage, double threshold){
	for(int i=1; i<originimage.rows-1; i++){
		for(int j=1; j<originimage.cols-1; j++){
			double middle = originimage.at<uchar>(i,j);
			double left = originimage.at<uchar>(i,j-1);
			double right = originimage.at<uchar>(i,j+1);
			double top = originimage.at<uchar>(i-1,j);
			double down = originimage.at<uchar>(i+1,j);
			double lefttop = originimage.at<uchar>(i-1,j-1);
			double leftdown = originimage.at<uchar>(i+1,j-1);
			double righttop = originimage.at<uchar>(i-1,j+1);
			double rightdown = originimage.at<uchar>(i+1,j+1);

			double temp = (double)(middle*(-8)+left+right+top+down+lefttop+leftdown+righttop+rightdown)/(double)3;
			if(temp>=threshold){
				resultimage.at<uchar>(i,j) = 0;
			}else{
				resultimage.at<uchar>(i,j) = 255;
			}
		}
	}
}

void MinimumVarianceLaplacian(Mat originimage, Mat resultimage, double threshold){
	for(int i=1; i<originimage.rows-1; i++){
		for(int j=1; j<originimage.cols-1; j++){
			double middle = originimage.at<uchar>(i,j);
			double left = originimage.at<uchar>(i,j-1);
			double right = originimage.at<uchar>(i,j+1);
			double top = originimage.at<uchar>(i-1,j);
			double down = originimage.at<uchar>(i+1,j);
			double lefttop = originimage.at<uchar>(i-1,j-1);
			double leftdown = originimage.at<uchar>(i+1,j-1);
			double righttop = originimage.at<uchar>(i-1,j+1);
			double rightdown = originimage.at<uchar>(i+1,j+1);

			double temp = (double)((-4)*middle+(-1)*(left+right+top+down)+2*(leftdown+lefttop+rightdown+righttop))/(double)3;
			if(temp>=threshold){
				resultimage.at<uchar>(i,j) = 0;
			}else{
				resultimage.at<uchar>(i,j) = 255;
			}
		}
	}
}

void LaplacianOfGaussian(Mat originimage, Mat resultimage, double threshold){

	for(int i=0; i<resultimage.rows; i++){
		for(int j=0; j<resultimage.cols; j++){
			resultimage.at<uchar>(i,j) = 255;
		}
	}

	for(int i=5; i<originimage.rows-5; i++){
		for(int j=5; j<originimage.cols-5; j++){
			double middle = originimage.at<uchar>(i,j);
			double left = originimage.at<uchar>(i,j-1);
			double right = originimage.at<uchar>(i,j+1);
			double top = originimage.at<uchar>(i-1,j);
			double down = originimage.at<uchar>(i+1,j);
			double ll = originimage.at<uchar>(i,j-2);
			double lll = originimage.at<uchar>(i,j-3);
			double llll = originimage.at<uchar>(i,j-4);
			double lllll = originimage.at<uchar>(i,j-5);
			double rr = originimage.at<uchar>(i,j+2);
			double rrr = originimage.at<uchar>(i,j+3);
			double rrrr = originimage.at<uchar>(i,j+4);
			double rrrrr = originimage.at<uchar>(i,j+5);
			double tt = originimage.at<uchar>(i-2,j);
			double ttt = originimage.at<uchar>(i-3,j);
			double tttt = originimage.at<uchar>(i-4,j);
			double ttttt = originimage.at<uchar>(i-5,j);
			double dd = originimage.at<uchar>(i+2,j);
			double ddd = originimage.at<uchar>(i+3,j);
			double dddd = originimage.at<uchar>(i+4,j);
			double ddddd = originimage.at<uchar>(i+5,j);
			
			double lefttop = 52*originimage.at<uchar>(i-1,j-1)-14*originimage.at<uchar>(i-2,j-1)-14*originimage.at<uchar>(i-1,j-2)-22*originimage.at<uchar>(i-3,j-1)
				-22*originimage.at<uchar>(i-1,j-3)-24*originimage.at<uchar>(i-2,j-2)-8*originimage.at<uchar>(i-1,j-4)-8*originimage.at<uchar>(i-4,j-1)
				-15*originimage.at<uchar>(i-2,j-3)-15*originimage.at<uchar>(i-3,j-2)-originimage.at<uchar>(i-1,j-5)-originimage.at<uchar>(i-5,j-1)
				-4*originimage.at<uchar>(i-2,j-4)-4*originimage.at<uchar>(i-4,j-2)-7*originimage.at<uchar>(i-3,j-3)
				-originimage.at<uchar>(i-5,j-2)-originimage.at<uchar>(i-2,j-5)-2*originimage.at<uchar>(i-3,j-4)-2*originimage.at<uchar>(i-4,j-3);
			double leftdown = 52*originimage.at<uchar>(i+1,j-1)-14*originimage.at<uchar>(i+2,j-1)-14*originimage.at<uchar>(i+1,j-2)-22*originimage.at<uchar>(i+3,j-1)
				-22*originimage.at<uchar>(i+1,j-3)-24*originimage.at<uchar>(i+2,j-2)-8*originimage.at<uchar>(i+1,j-4)-8*originimage.at<uchar>(i+4,j-1)
				-15*originimage.at<uchar>(i+2,j-3)-15*originimage.at<uchar>(i+3,j-2)-originimage.at<uchar>(i+1,j-5)-originimage.at<uchar>(i+5,j-1)
				-4*originimage.at<uchar>(i+2,j-4)-4*originimage.at<uchar>(i+4,j-2)-7*originimage.at<uchar>(i+3,j-3)
				-originimage.at<uchar>(i+5,j-2)-originimage.at<uchar>(i+2,j-5)-2*originimage.at<uchar>(i+3,j-4)-2*originimage.at<uchar>(i+4,j-3);
			double righttop = 52*originimage.at<uchar>(i-1,j+1)-14*originimage.at<uchar>(i-2,j+1)-14*originimage.at<uchar>(i-1,j+2)-22*originimage.at<uchar>(i-3,j+1)
				-22*originimage.at<uchar>(i-1,j+3)-24*originimage.at<uchar>(i-2,j+2)-8*originimage.at<uchar>(i-1,j+4)-8*originimage.at<uchar>(i-4,j+1)
				-15*originimage.at<uchar>(i-2,j+3)-15*originimage.at<uchar>(i-3,j+2)-originimage.at<uchar>(i-1,j+5)-originimage.at<uchar>(i-5,j+1)
				-4*originimage.at<uchar>(i-2,j+4)-4*originimage.at<uchar>(i-4,j+2)-7*originimage.at<uchar>(i-3,j+3)
				-originimage.at<uchar>(i-5,j+2)-originimage.at<uchar>(i-2,j+5)-2*originimage.at<uchar>(i-3,j+4)-2*originimage.at<uchar>(i-4,j+3);
			double rightdown = 52*originimage.at<uchar>(i+1,j+1)-14*originimage.at<uchar>(i+2,j+1)-14*originimage.at<uchar>(i+1,j+2)-22*originimage.at<uchar>(i+3,j+1)
				-22*originimage.at<uchar>(i+1,j+3)-24*originimage.at<uchar>(i+2,j+2)-8*originimage.at<uchar>(i+1,j+4)-8*originimage.at<uchar>(i+4,j+1)
				-15*originimage.at<uchar>(i+2,j+3)-15*originimage.at<uchar>(i+3,j+2)-originimage.at<uchar>(i+1,j+5)-originimage.at<uchar>(i+5,j+1)
				-4*originimage.at<uchar>(i+2,j+4)-4*originimage.at<uchar>(i+4,j+2)-7*originimage.at<uchar>(i+3,j+3)
				-originimage.at<uchar>(i+5,j+2)-originimage.at<uchar>(i+2,j+5)-2*originimage.at<uchar>(i+3,j+4)-2*originimage.at<uchar>(i+4,j+3);

			double temp = leftdown+lefttop+rightdown+righttop+178*middle+103*right+103*left+103*top+103*down
				-ll-23*lll-9*llll-2*lllll-rr-23*rrr-9*rrrr-2*rrrrr-tt-23*ttt-9*tttt-2*ttttt-dd-23*ddd-9*dddd-2*ddddd;

			if(temp>=threshold){
				resultimage.at<uchar>(i,j) = 0;
			}else{
				resultimage.at<uchar>(i,j) = 255;
			}
		}
	}
}

void DifferenceOfGaussian(Mat originimage, Mat resultimage, double threshold){

	for(int i=0; i<resultimage.rows; i++){
		for(int j=0; j<resultimage.cols; j++){
			resultimage.at<uchar>(i,j) = 255;
		}
	}

	for(int i=5; i<originimage.rows-5; i++){
		for(int j=5; j<originimage.cols-5; j++){
			double kernel[11][11] = {
				{-1,-3,-4,-6,-7,-8,-7,-6,-4,-3,-1},
				{-3,-5,-8,-11,-13,-13,-13,-11,-8,-5,-3},
				{-4,-8,-12,-16,-17,-17,-17,-16,-12,-8,-4},
				{-6,-11,-16,-16,0,15,0,-16,-16,-11,-6},
				{-7,-13,-17,0,85,160,85,0,-17,-13,-7},
				{-8,-13,-17,15,160,283,160,15,-17,-13,-8},
				{-7,-13,-17,0,85,160,85,0,-17,-13,-7},
				{-6,-11,-16,-16,0,15,0,-16,-16,-11,-6},
				{-4,-8,-12,-16,-17,-17,-17,-16,-12,-8,-4},
				{-3,-5,-8,-11,-13,-13,-13,-11,-8,-5,-3},
				{-1,-3,-4,-6,-7,-8,-7,-6,-4,-3,-1}};

			double temp = 0;
			for(int k=i-5 ; k<i+6; k++){
				for(int l=j-5; l<j+6; l++){
					temp = temp + kernel[k-i+5][l-j+5]*(double)originimage.at<uchar>(k,l);
				}
			}
			if(temp>=threshold){
				resultimage.at<uchar>(i,j) = 255;
			}else{
				resultimage.at<uchar>(i,j) = 0;
			}
		}
	}
}