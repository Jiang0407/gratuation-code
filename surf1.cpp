/**
 * @file SURF_Homography
 * @brief SURF detector + descriptor + FLANN Matcher + FindHomography
 * @author A. Huaman
 */

#include <stdio.h>
#include <iostream>
#include <cv.h>
#include "opencv2/core/core.hpp"
#include <opencv2/opencv.hpp>  
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/calib3d/calib3d.hpp"
#include "opencv2/nonfree/features2d.hpp"
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/nonfree/nonfree.hpp>  

#include "mex.h"

using namespace cv;
using namespace std;

void mexFunction(int nout, mxArray *out[], int nin, const mxArray *in[])
{
    // Get the name of the image
//     int len = mxGetNumberOfElements(in[0]);
//     std::string imgName;
//     imgName.resize( len + 1 );    
//     mxGetString( in[0], &imgName[0], len + 1 );
    
	cout<<"测试SURF算法"<<endl;  
	initModule_nonfree();//初始化模块，使用SIFT或SURF时用到  
	Ptr<FeatureDetector> detector = FeatureDetector::create( "SURF" );//创建SIFT特征检测器，可改成SURF/ORB
	Ptr<DescriptorExtractor> descriptor_extractor = DescriptorExtractor::create( "SURF" );//创建特征向量生成器，可改成SURF/ORB
	Ptr<DescriptorMatcher> descriptor_matcher = DescriptorMatcher::create( "BruteForce" );//创建特征匹配器  


	//读入图像  
	Mat img1 = cvLoadImage("E:\\程序\\code\\images\\case26\\new_img\\411\\main_building\\WPL_100.bmp");  
	Mat img2 = cvLoadImage("E:\\程序\\code\\images\\case26\\new_img\\411\\grass\\WP_003.jpg"); 

	//特征点检测  

	vector<KeyPoint> m_LeftKey,m_RightKey;  
	detector->detect( img1, m_LeftKey );//检测img1中的SIFT特征点，存储到m_LeftKey中  
	detector->detect( img2, m_RightKey );  	


	KeyPoint tmp;
	double *keyMat=new double[7]; 
	//根据特征点计算特征描述子矩阵，即特征向量矩阵  
	Mat descriptors1,descriptors2;   
	descriptor_extractor->compute( img1, m_LeftKey, descriptors1 );  
	descriptor_extractor->compute( img2, m_RightKey, descriptors2 );  
    
	//画出特征点  
	Mat img_m_LeftKey,img_m_RightKey;  
	drawKeypoints(img1,m_LeftKey,img_m_LeftKey,Scalar::all(-1),0);  
	drawKeypoints(img2,m_RightKey,img_m_RightKey,Scalar::all(-1),0);  
	imshow("Src1",img_m_LeftKey);  
	imshow("Src2",img_m_RightKey);  

	//特征匹配  
	vector<DMatch> matches;//匹配结果  
	descriptor_matcher->match( descriptors1, descriptors2, matches );//匹配两个图像的特征矩阵  


	//计算匹配结果中距离的最大和最小值  
	//距离是指两个特征向量间的欧式距离，表明两个特征的差异，值越小表明两个特征点越接近  
	double max_dist = 0;  
	double min_dist = 100;  
	for(int i=0; i<matches.size(); i++)  
	{  
		double dist = matches[i].distance;  
		if(dist < min_dist) min_dist = dist;  
		if(dist > max_dist) max_dist = dist;  
	}  
	

	//筛选出较好的匹配点  
	vector<DMatch> goodMatches;  
	for(int i=0; i<matches.size(); i++)  
	{  
		if(matches[i].distance < 0.3 * max_dist)  
		{  
			goodMatches.push_back(matches[i]);  
		}  
	}  

	//画出匹配结果  
	Mat img_matches;  
	//红色连接的是匹配的特征点对，绿色是未匹配的特征点  
	drawMatches(img1,m_LeftKey,img2,m_RightKey,goodMatches,img_matches,  
		Scalar::all(-1)/*CV_RGB(255,0,0)*/,CV_RGB(0,255,0),Mat(),2);  

	imshow("MatchSURF",img_matches);  
	IplImage result=img_matches; 


	//RANSAC匹配过程
	vector<DMatch> m_Matches=goodMatches;
	// 分配空间
	int ptCount = (int)m_Matches.size();
	Mat p1(ptCount, 2, CV_32F);
	Mat p2(ptCount, 2, CV_32F);

	// 把Keypoint转换为Mat
	Point2f pt;
	for (int i=0; i<ptCount; i++)
	{
		pt = m_LeftKey[m_Matches[i].queryIdx].pt;
		p1.at<float>(i, 0) = pt.x;
		p1.at<float>(i, 1) = pt.y;

		pt = m_RightKey[m_Matches[i].trainIdx].pt;
		p2.at<float>(i, 0) = pt.x;
		p2.at<float>(i, 1) = pt.y;
	}

	// 用RANSAC方法计算F
	Mat m_Fundamental;
	vector<uchar> m_RANSACStatus;       // 这个变量用于存储RANSAC后每个点的状态
	findFundamentalMat(p1, p2, m_RANSACStatus, FM_RANSAC,0.5,0.01);

	// 计算野点个数
	int OutlinerCount = 0;
	for (int i=0; i<ptCount; i++)
	{
		if (m_RANSACStatus[i] == 0)    // 状态为0表示野点
		{
			OutlinerCount++;
		}
	}
	int InlinerCount = ptCount - OutlinerCount;   // 计算内点
	


	// 这三个变量用于保存内点和匹配关系
	vector<Point2f> m_LeftInlier;
	vector<Point2f> m_RightInlier;
	vector<DMatch> m_InlierMatches;

	m_InlierMatches.resize(InlinerCount);
	m_LeftInlier.resize(InlinerCount);
	m_RightInlier.resize(InlinerCount);
	InlinerCount=0;
	float inlier_minRx=img1.cols;        //用于存储内点中右图最小横坐标，以便后续融合
	for (int i=0; i<ptCount; i++)
	{
		if (m_RANSACStatus[i] != 0)
		{
			m_LeftInlier[InlinerCount].x = p1.at<float>(i, 0);
			m_LeftInlier[InlinerCount].y = p1.at<float>(i, 1);
			m_RightInlier[InlinerCount].x = p2.at<float>(i, 0);
			m_RightInlier[InlinerCount].y = p2.at<float>(i, 1);
			m_InlierMatches[InlinerCount].queryIdx = InlinerCount;
			m_InlierMatches[InlinerCount].trainIdx = InlinerCount;
			if(m_RightInlier[InlinerCount].x<inlier_minRx) inlier_minRx=m_RightInlier[InlinerCount].x;   //存储内点中右图最小横坐标
			InlinerCount++;
		}
	}
	// 把内点转换为drawMatches可以使用的格式
	vector<KeyPoint> key1(InlinerCount);
	vector<KeyPoint> key2(InlinerCount);
	KeyPoint::convert(m_LeftInlier, key1);
	KeyPoint::convert(m_RightInlier, key2);
	// 显示计算F过后的内点匹配
	Mat OutImage;
	drawMatches(img1, key1, img2, key2, m_InlierMatches, OutImage);
	cvNamedWindow( "Match features", 1);
	cvShowImage("Match features", &IplImage(OutImage));
	waitKey(0);

	cvDestroyAllWindows();

	//矩阵H用以存储RANSAC得到的单应矩阵
	Mat H = findHomography( m_LeftInlier, m_RightInlier, RANSAC );
    
    int rows = H.rows;  
    int cols = H.cols;  
    out[0] = mxCreateDoubleMatrix(cols, rows, mxREAL);  
    double *tmpH;  
    tmpH = mxGetPr(out[0]);  
    for (int i = 0; i < rows; i++)  
        for (int j = 0; j < cols; j++)  
           // *(imgMat + i + j * rows) = (double)descriptors1.at<uchar>(i, j);  
             *(tmpH + i * cols + j ) = (double)H.at<double>(i, j);  
}