function [warp_image,sti_mask,b_mask]=SPHP_stitching_mine(in_name, num,s_mask,b_mask)

ref=1;
tar=2;
edge_list = [1,2];
img_n = size(in_name, 1);
edge_n = size(edge_list, 1);

% check
for ei = 1 : edge_n
    i = edge_list(ei, 1);
    j = edge_list(ei, 2);
    if i == j, fprintf('Index pair error.\n');  end
end

% load
I = cell(img_n, 1);
for i = 1 : img_n
     I{i} = imread(in_name{i}); 
end
% preprocessing
% T{i} is the coor dinate transform of image i that change from the image
% coordinate (origin at the upper-left corner of the image, x towards right, y towards down) to
% the standard coordinate (origin at the center of the image, x towards
% right, y towards up) 

imgw = zeros(img_n, 1);
imgh = zeros(img_n, 1);
resized_in_name = cell(img_n, 1);
T = cell(img_n, 1);
for i = 1 : img_n   
    resized_in_name{i} = sprintf('%s_resized%s', in_name{i}(1:end-4), in_name{i}(end-3:end));
    imwrite(I{i}, resized_in_name{i});
    imgw(i) = size(I{i}, 2);
    imgh(i) = size(I{i}, 1);
    T{i} = [1  0 -(imgw(i)+1)/2; ...
            0 -1  (imgh(i)+1)/2; ...
            0  0  1];
end
for i = 1 : img_n
    I{i} = im2double(I{i});
end

%% compute homography of image pairs
H = cell(img_n, img_n); % H{i,j} is the transformation from I{j} to I{i} (using standard coordinates)
for i = 1 : img_n
    H{i, i} = eye(3);
end

for ei = 1 : edge_n
    i = edge_list(ei, 1);
    j = edge_list(ei, 2);
%     tmpH = surf1();%sift_mosaic(I{i}, I{j},num); % tmpH is the transform from I{i} to I{j} (using image coordinates)
    tmpH = sift_mosaic_1(I{i}, I{j},num); % tmpH is the transform from I{i} to I{j} (using image coordinates)
%    tmpH=surf128_mosaic(I{i}, I{j},in_name{1},in_name{2});
    H{j, i} = T{j} * tmpH * inv(T{i}); % change of coord.
    H{i, j} = inv(H{j, i});
end 
mdltpara = [];
% save('TT.mat','tmpH','I','H');
% normalize all homography s.t. h9 = 1
for hi = 1 : numel(H)
    if hi==2 || hi==3
        H{hi} = H{hi} / H{hi}(3, 3);
    end
end

%% let H0 be the homography from reference to target, extract its parameters for computing our warp
H0 = H{tar, ref};
A = H0(1:2, 1:2);
t = H0(1:2, 3);
h31 = H0(3, 1);
h32 = H0(3, 2);
B = A - t*[h31 h32];
c = sqrt(h31^2+h32^2);
theta = atan2(-h32, -h31); % for compute ub1 and ub2

%% compute ub1 and ub2
H12=H{1,2};
T1=T{1,1};
T2=T{2,1}; 

% close all;
% clear all;
% 
% load('1.mat')
[ub1,ub2] = calculateU1U2(img_n,H12, c, theta,imgw,imgh,T1,T2);
 
fprintf('ub1 = %f, ub2 = %f \n',ub1,ub2); 
%-----------
%将ub1和ub2赋值为算法改动后的值
%-----------
%-----------
% fprintf('init (u1,u2) = (%.0f,%.0f). offset = (%.0f, %.0f).\n', oriub1, oriub2, offset_table1(idx1, idx2), offset_table2(idx1, idx2));


%% Compute parameters/coefficients of our warp
c1para = compute_c1_warp_coeff(H0, t, c, theta, ub1, ub2); %ref-->tar

%% c1 warp for each image
% [c1out, c1omask, Td] = texture_mapping(resized_in_name, imgw, imgh, img_grid_size, T, warp_type, H, ref, tar, c1para, mdltpara, 0, 'white');
%---------------------------

%  close all
%  clear all
%  clc
%  load('5.mat')
if num==1
    [c1out,c1omask] = texture_mapping_mine(resized_in_name, imgw, imgh,T, H, ref, tar, c1para, mdltpara, num);
else
    [c1out,c1omask,ble_mask] = texture_mapping_mine(resized_in_name, imgw, imgh, T, H, ref, tar, c1para,mdltpara, num, s_mask,b_mask);
end
%----------------------------num
% for i = 1 : img_n
%     figure; imshow(c1out{i}); % the warped result of each image
% end

%---------------------
%SPHP效果融合-直接平均法
%---------------------
% img_result=imageblending(c1out{1},c1out{2});
% figure;imshow(img_result);
%---------------------
%融合
%---------------------
% [warp_image,~,~] = imgblending_mine(c1out,'global');
% close all;
% clear all;
% clc;
% load('4.mat')
endb=boundb(c1out{1});
if num==1
    [warp_image,sti_mask,b_mask] = imgblending_mine(c1out,'regional',500,endb,c1omask,num,I); %global为标志字符串，进行全局（global）or局部（original）融合。
else
    [warp_image,sti_mask,b_mask] = imgblending_mine(c1out,'regional',50,endb,c1omask,num,I,ble_mask); %global为标志字符串，进行全局（global）or局部（original）融合。
end
warp_image=uint8(warp_image);
% %%%%%%%%%%%%%%%%%%%
% %测试-直接平均法+线性过渡法
% %%%%%%%%%%%%%%%%%%%
% outcanvas=imageblending1(c1out{1},c1out{2},c1omask{1}(:,:,1),c1omask{2}(:,:,1));
% %%
% warp_image=uint8(warp_image);
% figure(3);imshow(warp_image);
%---------------------
% for i = 1 : img_n
%     eval(['delete ' resized_in_name{i}]);
% end
% 
% %% composite by linear blending
% for i = 1 : img_n
%     if i == 1
%         out = c1out{1};
%         out_mask = c1omask{1};
%         % center of out
%         [r, c] = find(c1omask{1}(:, :, 1));
%         out_center = [mean(r) mean(c)];
%     else % blend out and c1out{i}
%         % center of c1out{i}
%         [r, c] = find(c1omask{i}(:, :, 1));
%         out_i_center = [mean(r) mean(c)];
%         % compute weighted mask
%         vec = out_i_center - out_center; % vector from out_center to out_i_center
%         intsct_mask = c1omask{i}(:, :, 1) & out_mask(:, :, 1); % 1 channel
%         out_only_mask = out_mask(:, :, 1) - intsct_mask;
%         out_i_only_mask = c1omask{i}(:, :, 1) - intsct_mask;
%         
%         [r, c] = find(intsct_mask(:, :, 1));
%         idx = sub2ind(size(c1omask{i}(:, :, 1)), r, c);
%         out_wmask = zeros(size(c1omask{i}(:, :, 1)));
%         proj_val = (r - out_center(1))*vec(1) + (c- out_center(2))*vec(2); % inner product
%         out_wmask(idx) = (proj_val - (min(proj_val)+(1e-3))) / ...
%                          ((max(proj_val)-(1e-3)) - (min(proj_val)+(1e-3))); % weight map (of overlapped area) for c1out{i}, 1 channel
%         % blending
%         mask1 = out_mask(:, :, 1)&(out_wmask==0);
%         mask2 = out_wmask;
%         mask3 = c1omask{i}(:, :, 1)&(out_wmask==0);
%         mask1 = cat(3, mask1, mask1, mask1); mask2 = cat(3, mask2, mask2, mask2); mask3 = cat(3, mask3, mask3, mask3);
%         out = out.*(mask1+(1-mask2).*(mask2~=0)) + c1out{i}.*(mask2+mask3);
%         % update
%         out_mask = out_mask | c1omask{i};
%         out_center = out_i_center; % update out_center by assign center of c1out{i}
%     end
% end
% c1all = out;
% denom = zeros(size(c1out{1}));
% for i = 1 : img_n
%     denom = denom + c1omask{i};
% end
% bgcolor = 'white'; % background color
% if strcmp(bgcolor, 'black')
%     c1all(denom==0) = 1;
% else
%     c1all(denom==0) = 0;
% end
% warp_image=uint8(c1all*255);
% %%
% figure(4); imshow(warp_image); % the final result
% %%
% for i = 1 : img_n
%     imwrite(c1out{i}, sprintf('%s_out%d.png', warp_type, i));
% end
% imwrite(c1all, [warp_type '_outall.png']);