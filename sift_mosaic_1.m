function Hg=sift_mosaic_1(img1,img2,num)
% save 0517-1.mat
% close all;
% clear all;
% clc;
% load('0517-1.mat')
addpath('modelspecific');
addpath('mexfiles');
addpath('multigs');
clear global;
global fitfn resfn degenfn psize numpar
fitfn = 'homography_fit';
resfn = 'homography_res';
degenfn = 'homography_degen';
psize   = 4;
numpar  = 9;

M     = 500;  % Number of hypotheses for RANSAC.
thr   = 0.1;  % RANSAC threshold.

%% make grayscale
im1 = im2single(img1) ;
im2 = im2single(img2) ;
if size(im1,3) > 1, im1g = rgb2gray(im1) ; else im1g = im1 ; end
if size(im2,3) > 1, im2g = rgb2gray(im2) ; else im2g = im2 ; end

%% SIFT matches
% fprintf('  Keypoint detection and matching...');tic;
% [ kp1,ds1 ] = vl_sift(im1g);
% [ kp2,ds2 ] = vl_sift(im2g);

row_img=size(img1,2);
vo_sift=floor((1-1/num)*row_img+1); %对于大图像取一半进行配准
[ kp1,ds1 ] = vl_sift( single( rgb2gray(img1(:,vo_sift:end,:) ) ),'PeakThresh', 0,'edgethresh',200);
kp1(1,:)=kp1(1,:)+vo_sift-1; %准确坐标
%kp2 第一维是列， 第二维是行
[ kp2,ds2 ] = vl_sift(single(rgb2gray(img2)),'PeakThresh', 0,'edgethresh',200);    

matches   = vl_ubcmatch(ds1,ds2); 

%% Normalise point distribution.
% fprintf('  Normalising point distribution...');tic;
data_orig = [ kp1(1:2,matches(1,:)) ; ones(1,size(matches,2)) ; kp2(1:2,matches(2,:)) ; ones(1,size(matches,2)) ];
[ dat_norm_img1,T1 ] = normalise2dpts(data_orig(1:3,:));
[ dat_norm_img2,T2 ] = normalise2dpts(data_orig(4:6,:));
data_norm = [ dat_norm_img1 ; dat_norm_img2 ];

% fprintf('done (%fs)\n',toc);
% fprintf('Outlier removal\n');tic;

%% Multi-GS
rng(0);
[ ~,res,~,~ ] = multigsSampling(100,data_norm,M,10);
con = sum(res<=thr);
[ ~, maxinx ] = max(con);
inliers = find(res(:,maxinx)<=thr);
outliers = find(res(:,maxinx)>thr);

if size(img1,1) == size(img2,1)
    % Show results of RANSAC.
%     fprintf('  Showing results of RANSAC...');tic;
    figure;
    imshow([img1 img2]);
    hold on;
    plot(data_orig(1,:),data_orig(2,:),'ro','LineWidth',2);
    plot(data_orig(4,:)+size(img1,2),data_orig(5,:),'ro','LineWidth',2);

    
     
%     for i=1:1:length(outliers)
%         plot(data_orig(1,outliers(i)),data_orig(2,outliers(i)),'ro','LineWidth',1);
%         plot(data_orig(4,outliers(i))+size(img1,2),data_orig(5,outliers(i)),'ro','LineWidth',1);
%         plot([data_orig(1,outliers(i)) data_orig(4,outliers(i))+size(img1,2)],[data_orig(2,outliers(i)) data_orig(5,inliers(i))],'r-','LineWidth',2);
%     end
   % 37 39 1 2 3
%     inliers=inliers([1:19 38 40:45 47 49 51 55 56 58 59 61:64 68 70:74 76:79 80 81 84:87 89 92:95 97 98 103 105:end]);
    for i=1:1:length(inliers)
        plot(data_orig(1,inliers(i)),data_orig(2,inliers(i)),'go','LineWidth',1);
        plot(data_orig(4,inliers(i))+size(img1,2),data_orig(5,inliers(i)),'go','LineWidth',1);
        plot([data_orig(1,inliers(i)) data_orig(4,inliers(i))+size(img1,2)],[data_orig(2,inliers(i)) data_orig(5,inliers(i))],'g-');
    end 
    title('Ransac''s results');
    fprintf('done (%fs)\n',toc);
end

%-----------------------
% Global homography (H).
%-----------------------
% fprintf('DLT (projective transform) on inliers\n');
% Refine homography using DLT on inliers.
% fprintf('> Refining homography (H) using DLT...');tic;
% % data_517 = data_orig(:,inliers);
% % data_orig = distortion_elimination(data_517);
% % [ dat_norm_img1,T1 ] = normalise2dpts(data_orig(1:3,:));
% % [ dat_norm_img2,T2 ] = normalise2dpts(data_orig(4:6,:));
% % data_norm = [ dat_norm_img1 ; dat_norm_img2 ];
[ h,A,D1,D2 ] = feval(fitfn,data_norm(:,inliers));
Hg = T2\(reshape(h,3,3)*T1);
% fprintf('done (%fs)\n',toc);