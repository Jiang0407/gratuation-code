/*
 * A Demo to OpenCV Implementation of SURF
 * Further Information Refer to "SURF: Speed-Up Robust Feature"
 * Author: Liu Liu
 * liuliu.1987+opencv@gmail.com
 */
#include "opencv2/objdetect/objdetect.hpp"
#include "opencv2/features2d/features2d.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/calib3d/calib3d.hpp"
#include "opencv2/nonfree/nonfree.hpp"
#include "opencv2/imgproc/imgproc_c.h"
#include "opencv2/legacy/legacy.hpp"
#include "opencv2/legacy/compat.hpp"

#include <iostream>
#include <vector>
#include <stdio.h>
#include "mex.h"

using namespace std;

void mexFunction(int nout, mxArray *out[], int nin, const mxArray *in[])
{ 
    // Get the name of the image
    int nStringLen;
    nStringLen = mxGetNumberOfElements(in[0]);
    std::string szImageName;
    szImageName.resize( nStringLen + 1 );
    
    mxGetString( in[0], &szImageName[0], nStringLen + 1 );
      
    cv::initModule_nonfree();
    IplImage* object = cvLoadImage( szImageName.c_str(), CV_LOAD_IMAGE_GRAYSCALE ); 
    if( !object)
    { 
        exit(-1);
    }
//     else
//     {
//         cvNamedWindow("1");
//         cvShowImage("1",object);
//     }
    CvMemStorage* storage = cvCreateMemStorage(0);
//     IplImage* object_color = cvCreateImage(cvGetSize(object), 8, 3);
//     cvCvtColor( object, object_color, CV_GRAY2BGR );

    CvSeq* objectKeypoints = 0, *objectDescriptors = 0; 
    
    int i;
    CvSURFParams params = cvSURFParams(500, 1);

    double tt = (double)cvGetTickCount();
    cvExtractSURF( object, 0, &objectKeypoints, &objectDescriptors, storage, params );
//     printf("Object Descriptors: %d\n", objectDescriptors->total);

    //返回特征点矩阵
    int rows = objectKeypoints ->total;  
    int cols = 4;  
    CvSURFPoint *r = (CvSURFPoint*)cvGetSeqElem(objectKeypoints,i);
    out[0] = mxCreateDoubleMatrix(cols, rows, mxREAL);  
    double *k;  
    k = mxGetPr(out[0]);     
	CvSeqReader reader; 
	cvStartReadSeq(objectKeypoints, &reader, 0);
	for(int i = 0; i < objectKeypoints->total; i++)
	{
		  CvSURFPoint *r = (CvSURFPoint*)cvGetSeqElem(objectKeypoints,i);
          *(k + i*4 )= r->pt.x;
          *(k + i*4 + 1)= r->pt.y;
          *(k + i*4 + 2)= r->size;
          *(k + i*4 + 3)= r->dir; 
	}
    
    //返回特征描述矩阵
    rows = objectDescriptors->total;  
    cols = 128;  
    out[1] = mxCreateDoubleMatrix(cols, rows, mxREAL);  
    double *ds;  
    ds = mxGetPr(out[1]);      
	cvStartReadSeq(objectDescriptors, &reader, 0);
	for(int i = 0; i < objectDescriptors->total; i++)
	{
		  const float* descriptor = (const float*)reader.ptr; // descriptor 指向第i 特征描述符（float数组）的指针
          // 数组长度为:reader.seq->elem_size/sizeof(float) or objectDescriptors->elem_size/sizeof(float) 
		  for (int j=0;j<127;j++)
          {
              *(ds + i*128 + j)=descriptor[j];	
          }
		  CV_NEXT_SEQ_ELEM(reader.seq->elem_size, reader); //读取下一个特征点
	}
    
//     cvWaitKey(0); 
}

