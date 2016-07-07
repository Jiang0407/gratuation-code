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
    
	cout<<"����SURF�㷨"<<endl;  
	initModule_nonfree();//��ʼ��ģ�飬ʹ��SIFT��SURFʱ�õ�  
	Ptr<FeatureDetector> detector = FeatureDetector::create( "SURF" );//����SIFT������������ɸĳ�SURF/ORB
	Ptr<DescriptorExtractor> descriptor_extractor = DescriptorExtractor::create( "SURF" );//���������������������ɸĳ�SURF/ORB
	Ptr<DescriptorMatcher> descriptor_matcher = DescriptorMatcher::create( "BruteForce" );//��������ƥ����  


	//����ͼ��  
	Mat img1 = cvLoadImage("E:\\����\\code\\images\\case26\\new_img\\411\\main_building\\WPL_100.bmp");  
	Mat img2 = cvLoadImage("E:\\����\\code\\images\\case26\\new_img\\411\\grass\\WP_003.jpg"); 

	//��������  

	vector<KeyPoint> m_LeftKey,m_RightKey;  
	detector->detect( img1, m_LeftKey );//���img1�е�SIFT�����㣬�洢��m_LeftKey��  
	detector->detect( img2, m_RightKey );  	


	KeyPoint tmp;
	double *keyMat=new double[7]; 
	//����������������������Ӿ��󣬼�������������  
	Mat descriptors1,descriptors2;   
	descriptor_extractor->compute( img1, m_LeftKey, descriptors1 );  
	descriptor_extractor->compute( img2, m_RightKey, descriptors2 );  
    
	//����������  
	Mat img_m_LeftKey,img_m_RightKey;  
	drawKeypoints(img1,m_LeftKey,img_m_LeftKey,Scalar::all(-1),0);  
	drawKeypoints(img2,m_RightKey,img_m_RightKey,Scalar::all(-1),0);  
	imshow("Src1",img_m_LeftKey);  
	imshow("Src2",img_m_RightKey);  

	//����ƥ��  
	vector<DMatch> matches;//ƥ����  
	descriptor_matcher->match( descriptors1, descriptors2, matches );//ƥ������ͼ�����������  


	//����ƥ�����о����������Сֵ  
	//������ָ���������������ŷʽ���룬�������������Ĳ��죬ֵԽС��������������Խ�ӽ�  
	double max_dist = 0;  
	double min_dist = 100;  
	for(int i=0; i<matches.size(); i++)  
	{  
		double dist = matches[i].distance;  
		if(dist < min_dist) min_dist = dist;  
		if(dist > max_dist) max_dist = dist;  
	}  
	

	//ɸѡ���Ϻõ�ƥ���  
	vector<DMatch> goodMatches;  
	for(int i=0; i<matches.size(); i++)  
	{  
		if(matches[i].distance < 0.3 * max_dist)  
		{  
			goodMatches.push_back(matches[i]);  
		}  
	}  

	//����ƥ����  
	Mat img_matches;  
	//��ɫ���ӵ���ƥ���������ԣ���ɫ��δƥ���������  
	drawMatches(img1,m_LeftKey,img2,m_RightKey,goodMatches,img_matches,  
		Scalar::all(-1)/*CV_RGB(255,0,0)*/,CV_RGB(0,255,0),Mat(),2);  

	imshow("MatchSURF",img_matches);  
	IplImage result=img_matches; 


	//RANSACƥ�����
	vector<DMatch> m_Matches=goodMatches;
	// ����ռ�
	int ptCount = (int)m_Matches.size();
	Mat p1(ptCount, 2, CV_32F);
	Mat p2(ptCount, 2, CV_32F);

	// ��Keypointת��ΪMat
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

	// ��RANSAC��������F
	Mat m_Fundamental;
	vector<uchar> m_RANSACStatus;       // ����������ڴ洢RANSAC��ÿ�����״̬
	findFundamentalMat(p1, p2, m_RANSACStatus, FM_RANSAC,0.5,0.01);

	// ����Ұ�����
	int OutlinerCount = 0;
	for (int i=0; i<ptCount; i++)
	{
		if (m_RANSACStatus[i] == 0)    // ״̬Ϊ0��ʾҰ��
		{
			OutlinerCount++;
		}
	}
	int InlinerCount = ptCount - OutlinerCount;   // �����ڵ�
	


	// �������������ڱ����ڵ��ƥ���ϵ
	vector<Point2f> m_LeftInlier;
	vector<Point2f> m_RightInlier;
	vector<DMatch> m_InlierMatches;

	m_InlierMatches.resize(InlinerCount);
	m_LeftInlier.resize(InlinerCount);
	m_RightInlier.resize(InlinerCount);
	InlinerCount=0;
	float inlier_minRx=img1.cols;        //���ڴ洢�ڵ�����ͼ��С�����꣬�Ա�����ں�
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
			if(m_RightInlier[InlinerCount].x<inlier_minRx) inlier_minRx=m_RightInlier[InlinerCount].x;   //�洢�ڵ�����ͼ��С������
			InlinerCount++;
		}
	}
	// ���ڵ�ת��ΪdrawMatches����ʹ�õĸ�ʽ
	vector<KeyPoint> key1(InlinerCount);
	vector<KeyPoint> key2(InlinerCount);
	KeyPoint::convert(m_LeftInlier, key1);
	KeyPoint::convert(m_RightInlier, key2);
	// ��ʾ����F������ڵ�ƥ��
	Mat OutImage;
	drawMatches(img1, key1, img2, key2, m_InlierMatches, OutImage);
	cvNamedWindow( "Match features", 1);
	cvShowImage("Match features", &IplImage(OutImage));
	waitKey(0);

	cvDestroyAllWindows();

	//����H���Դ洢RANSAC�õ��ĵ�Ӧ����
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