#include "mex.h"
#include <time.h>
#include <math.h>
#include <string.h>
#include <opencv2\opencv.hpp>
using namespace std;
using namespace cv;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int m, n;
	int a, b;
	double *sti_mask;
	double *ble_mask;
	m = (int)mxGetM(prhs[0]);
	n = (int)mxGetN(prhs[0])/3;
	sti_mask = mxGetPr(prhs[0]);
	ble_mask = mxGetPr(prhs[1]);

	//opencv读取视频
	IplImage *inframe1;
	IplImage *inframe2;
	inframe1 = cvLoadImage("WP_001.jpg");
	inframe2 = cvLoadImage("WP_002.jpg");
	//CvCapture *wp1 = cvCreateFileCapture("wp_1.mp4");
	//CvCapture *wp2 = cvCreateFileCapture("wp_2.mp4");
	//CvCapture *wp3 = cvCreateFileCapture("wp_3.mp4");
	cvNamedWindow("input", CV_WINDOW_AUTOSIZE);
	/*cvNamedWindow("output", 0);
	cvResizeWindow("output", 480, 360);*/
	/*int numframes = (int)cvGetCaptureProperty(wp1, CV_CAP_PROP_FRAME_COUNT);*/
	IplImage *outframe = cvCreateImage(cvSize(n, m), IPL_DEPTH_8U, 3);

	//查表
	for (int k = 0; k < 1;k++)
	{

		//运算时间
		double t1 = double(cvGetTickCount());
		double f = cvGetTickFrequency();
		//查表
		int in_step1 = inframe1->widthStep;
		int in_step2 = inframe2->widthStep;
		int out_step = outframe->widthStep;
		int channels = outframe->nChannels;
		for (int i = 0; i < m; i++)
		{
			for (int j = 0; j < n; j++)
				switch ((int)sti_mask[2*m*n + i*n + j])
			{
				case 1:
					//cvGrabFrame(wp1);
					//inframe1 = cvQueryFrame(wp1);
					if (!inframe1)  //没有它程序没办法跳出循环
					{
						break;
					}
					/*int channels = outframe->nChannels;*/
					outframe->imageData[i*out_step + j*channels + 0] = (inframe1->imageData[(int)sti_mask[i*n + j] * in_step1 + (int)sti_mask[m*n + i*n + j] * channels + 0]);
					outframe->imageData[i*out_step + j*channels + 1] = (inframe1->imageData[(int)sti_mask[i*n + j] * in_step1 + (int)sti_mask[m*n + i*n + j] * channels + 1]);
					outframe->imageData[i*out_step + j*channels + 2] = (inframe1->imageData[(int)sti_mask[i*n + j] * in_step1 + (int)sti_mask[m*n + i*n + j] * channels + 2]);
				case 2:
	/*				cvGrabFrame(wp2);
					inframe1 = cvQueryFrame(wp2);*/
					if (!inframe1)  //没有它程序没办法跳出循环
					{
						break;
					}
					/*int channels = outframe->nChannels;*/
					outframe->imageData[i*out_step + j*channels + 0] = inframe2->imageData[(int)sti_mask[i*n + j] * in_step2 + (int)sti_mask[m*n + i*n + j] * channels + 0];
					outframe->imageData[i*out_step + j*channels + 1] = inframe2->imageData[(int)sti_mask[i*n + j] * in_step2 + (int)sti_mask[m*n + i*n + j] * channels + 1];
					outframe->imageData[i*out_step + j*channels + 2] = inframe2->imageData[(int)sti_mask[i*n + j] * in_step2 + (int)sti_mask[m*n + i*n + j] * channels + 2];
				//case 3:
				//	cvGrabFrame(wp3);
				//	inframe1 = cvQueryFrame(wp3);
				//	if (!inframe1)  //没有它程序没办法跳出循环
				//	{
				//		break;
				//	}
				//	int in_step = inframe1->widthStep / sizeof(uchar);
				//	int out_step = outframe->widthStep / sizeof(uchar);
				//	int channels = outframe->nChannels;
				//	outframe->imageData[i*out_step + j*channels + 0] = inframe1->imageData[sti_mask[i*n + j] * in_step + sti_mask[m*n + i*n + j] * channels + 0];
				//	outframe->imageData[i*out_step + j*channels + 1] = inframe1->imageData[sti_mask[i*n + j] * in_step + sti_mask[m*n + i*n + j] * channels + 1];
				//	outframe->imageData[i*out_step + j*channels + 2] = inframe1->imageData[sti_mask[i*n + j] * in_step + sti_mask[m*n + i*n + j] * channels + 2];
				//case 12:
				//	cvGrabFrame(wp1);
				//	inframe1 = cvQueryFrame(wp1);
				//	int in_step1 = inframe1->widthStep / sizeof(uchar);
				//	cvGrabFrame(wp2);
				//	inframe2 = cvQueryFrame(wp2);
				//	if (!inframe1 || !inframe2)  //没有它程序没办法跳出循环
				//	{
				//		break;
				//	}
				//	int in_step2 = inframe2->widthStep / sizeof(uchar);
				//	int out_step = outframe->widthStep / sizeof(uchar);
				//	int channels = outframe->nChannels;
				//	outframe->imageData[i*out_step + j*channels + 0] = sti_mask[i*n + j] * inframe1->imageData[ble_mask[i*n + j] * in_step + ble_mask[m*n + i*n + j] * channels + 0] + sti_mask[m*n + i*n + j] * inframe2->imageData[ble_mask[2*m*n + i*n + j] * in_step + ble_mask[3*m*n + i*n + j] * channels + 0];
				//	outframe->imageData[i*out_step + j*channels + 0] = sti_mask[i*n + j] * inframe1->imageData[ble_mask[i*n + j] * in_step + ble_mask[m*n + i*n + j] * channels + 1] + sti_mask[m*n + i*n + j] * inframe2->imageData[ble_mask[2*m*n + i*n + j] * in_step + ble_mask[3*m*n + i*n + j] * channels + 1];
				//	outframe->imageData[i*out_step + j*channels + 0] = sti_mask[i*n + j] * inframe1->imageData[ble_mask[i*n + j] * in_step + ble_mask[m*n + i*n + j] * channels + 2] + sti_mask[m*n + i*n + j] * inframe2->imageData[ble_mask[2*m*n + i*n + j] * in_step + ble_mask[3*m*n + i*n + j] * channels + 2];
				default:
					break; 
			}
			cvShowImage("output", outframe);
			double t2 = double(cvGetTickCount());
			cout << "运行时间: " << (t2 - t1) / f / 1000 << "ms" << endl;
			cvWaitKey(10);
		}

	}

}