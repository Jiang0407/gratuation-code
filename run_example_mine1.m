function [warp_img,center00]=run_example_mine1(Hg,Hmdlt,invHmdlt,InmeshX,InmeshY,img1,img2,num,center)
% function warp_image=run_example_mine(Hg,Hmdlt,invHmdlt,InmeshX,InmeshY,img1,img2,num)
close all
run('vlfeat-0.9.14/toolbox/vl_setup');

% Example 1: temple
img_n = 2; % The number of input images
in_name = cell(2,1);
% in_name{1} = 'images/case26/WP_7.jpg';
% in_name{2} = 'images/case26/WP_8.jpg';
%------------
in_name{1} = img1;
in_name{2} = img2;
%------------
edge_list = [1,2]; % Each row represents an image pair to be aligned.
% In this example, there is only one pair (image 1 and image 2).
ref = 1; % the index of the reference image
tar = 2; % the index of the target image. Our warp is constructed from the homgraphy that maps from ref to tar
warp_type = 'ours'; % 'ours' for our warp. 'hom' for homography warp
zeroR_ON = 1; % whether we restrict the similarity warp in our warp to be no rotation (zeroR_ON=1) or not (zeroR_ON=0)
[warp_img,center00]=SPHP_stitching_mine(in_name, edge_list, ref, tar, warp_type, zeroR_ON, Hg, Hmdlt,invHmdlt,InmeshX, InmeshY,num,center); % run stitching
title('Example 1: temple. (Press any key to continue with the next example.)');

%------------------------
%Ö÷º¯Êý
%------------------------

    