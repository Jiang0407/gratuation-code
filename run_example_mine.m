function [warp_img,center00,sti_mask,ble_mask]=run_example_mine(img1,img2,num,center,s_mask,b_mask)
% function warp_image=run_example_mine(Hg,Hmdlt,invHmdlt,InmeshX,InmeshY,img1,img2,num)

 
in_name = cell(2,1);
% in_name{1} = 'images/case26/WP_7.jpg';
% in_name{2} = 'images/case26/WP_8.jpg';
%------------
in_name{1} = img1;
in_name{2} = img2; 

if num==1
   [warp_img,center00,sti_mask,ble_mask]=SPHP_stitching_mine(in_name, num,center); 
else
   [warp_img,center00,sti_mask,ble_mask]=SPHP_stitching_mine(in_name, num,center,s_mask,b_mask); % run stitching
end
% title('Example 1: temple. (Press any key to continue with the next example.)');

%------------------------
%Ö÷º¯Êý
%------------------------

    