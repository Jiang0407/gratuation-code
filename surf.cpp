/**
 * @file SURF_Homography
 * @brief SURF detector + descriptor + FLANN Matcher + FindHomography
 * @author A. Huaman
 */

#include <stdio.h>
#include <iostream>
#include <string>
#include <opencv2/opencv.hpp> 
// #include <opencv2/core/core.hpp> 
// #include <opencv2/highgui/highgui.hpp>
// #include <opencv2/calib3d/calib3d.hpp>
// #include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/nonfree/nonfree.hpp>  
#include <opencv2/nonfree/features2d.hpp>
#include "mex.h"
using namespace cv; 

void mexFunction(int nout, mxArray *out[], int nin, const mxArray *in[])
{
    // Get the name of the image
    int len = mxGetNumberOfElements(in[0]);
    std::string imgName;
    imgName.resize( len + 1 );    
    mxGetString( in[0], &imgName[0], len + 1 );
  
    
	initModule_nonfree();//初始化模块，使用SIFT或SURF时用到  
	Ptr<FeatureDetector> detector = FeatureDetector::create( "SURF" );//创建SIFT特征检测器，可改成SURF/ORB
	Ptr<DescriptorExtractor> descriptor_extractor = DescriptorExtractor::create( "SIFT" );//创建特征向量生成器，可改成SURF/ORB
	Ptr<DescriptorMatcher> descriptor_matcher = DescriptorMatcher::create( "BruteForce" );//创建特征匹配器  
	 
    Mat image= cv::imread(imgName );
    
	//特征点检测   
	vector<KeyPoint> m_LeftKey;  
	detector->detect((Mat)image, m_LeftKey );//检测img1中的SIFT特征点，存储到m_LeftKey中   

	//根据特征点计算特征描述子矩阵，即特征向量矩阵  
	Mat descriptors1;  
	descriptor_extractor->compute((Mat)image, m_LeftKey, descriptors1 );  
    
	// convert the result to Matlab-supported format for returning  
	
	int sizei = m_LeftKey.size();   
	int sizej=4;
	KeyPoint tmp;
	out[0] = mxCreateDoubleMatrix(sizej,sizei, mxREAL);  
	double *keyMat; 
    keyMat = mxGetPr(out[0]);  
    for (int i = 0; i < sizei; i++)   
    {
		tmp=m_LeftKey[i];
		*(keyMat + i * sizej) =   (double)tmp.pt.x; 
		*(keyMat + i * sizej+1) = (double)tmp.pt.y; 
		*(keyMat + i * sizej+2) = (double)tmp.size; 
		*(keyMat + i * sizej+3) = (double)tmp.angle;
	} 		
		
    int rows = descriptors1.rows;  
    int cols = descriptors1.cols;  
    out[1] = mxCreateDoubleMatrix(cols, rows, mxREAL);  
    double *imgMat;  
    imgMat = mxGetPr(out[1]);  
    for (int i = 0; i < rows; i++)  
        for (int j = 0; j < cols; j++)  
           // *(imgMat + i + j * rows) = (double)descriptors1.at<uchar>(i, j);  
             *(imgMat + i * cols + j ) = (double)descriptors1.at<uchar>(i, j);  
	
	
}

