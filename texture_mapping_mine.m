function [c1out,c1omask,ble_mask]= texture_mapping_mine(in_name,imgw, imgh, T, H, ref, tar, c1para,mdltpara,num,s_mask,b_mask)
img_grid_size=10;
mapping_location = 'texture_mapping/';

img_n = numel(in_name);
c1out = cell(img_n, 1);
c1omask = cell(img_n, 1);
outmesh_X = cell(img_n, 1);     outmesh_Y = cell(img_n, 1);
curr_minx = zeros(img_n, 1);    curr_miny = zeros(img_n, 1);
curr_maxx = zeros(img_n, 1);    curr_maxy = zeros(img_n, 1);

for i = 1 : 2 %img_n
    meshW = ceil(imgw(i) / img_grid_size);
    meshH = ceil(imgh(i) / img_grid_size);
    X0list = linspace(1, imgw(i), meshW);
    Y0list = linspace(1, imgh(i), meshH);
    [inmesh_X{i}, inmesh_Y{i}] = meshgrid(X0list, Y0list);

    [outmesh_X{i}, outmesh_Y{i}] = mesh_warp(img_grid_size, imgw(i), imgh(i), T, i, ref, tar, H, c1para, mdltpara);
    
    curr_minx(i) = min(outmesh_X{i}(:)); curr_miny(i) = min(outmesh_Y{i}(:));
    curr_maxx(i) = max(outmesh_X{i}(:)); curr_maxy(i) = max(outmesh_Y{i}(:));
end % up to now, the coord. system is x towards right, y towards up
all_minx = min(curr_minx); all_miny = min(curr_miny);
all_maxx = max(curr_maxx); all_maxy = max(curr_maxy);
% fprintf('x from %.2f ~ %.2f\n', all_minx, all_maxx);
% fprintf('y from %.2f ~ %.2f\n', all_miny, all_maxy);

%% reorganize
% change of coord. s.t. (1) the upperleft of bounding box of composite corresponds to (0,0)
%                       (2) x towards right, y-axis towards down
% load('0506.mat')
T1 = [1, 0, -all_minx; 0, -1, all_maxy; 0, 0, 1];
for i = 1 :2 %img_n
    [outmesh_X{i}, outmesh_Y{i}] = apply_transform(outmesh_X{i}, outmesh_Y{i}, T1);
    % recompute min & max value
    curr_minx(i) = min(outmesh_X{i}(:)); curr_miny(i) = min(outmesh_Y{i}(:));
    curr_maxx(i) = max(outmesh_X{i}(:)); curr_maxy(i) = max(outmesh_Y{i}(:));
end % now the coord. system is x towards right, y towards down. bounding box aligns with both axis
all_minx = min(curr_minx); all_miny = min(curr_miny);
all_maxx = max(curr_maxx); all_maxy = max(curr_maxy);
bbarea = (all_maxx-all_minx)*(all_maxy-all_miny); % the area of bounding box
%if all_minx~=0 || all_miny~=0, fprintf('%.12f %.12f\n', all_minx, all_miny); pause; end

% scale down to avoid memory issue
%%
% numerical post-processing (set near zero to be zero)
for i = 1 : 2 %img_n
    outmesh_X{i}( outmesh_X{i}<=(1e-7) ) = 0;
    outmesh_Y{i}( outmesh_Y{i}<=(1e-7) ) = 0;
end

% adding borders by translating the composite
%-------------
% %取消5个像素点的边缘
% tr_x = 5;
% tr_y = 5;
% T3 = [1, 0, tr_x; 0, 1, tr_y; 0, 0, 1];
for i = 1 : 2 %img_n
    outmesh_X{i} = outmesh_X{i}+2; %令i,j从2开始.不是从0开始为了适应MATLAB；不从1开始为了融合方便
    outmesh_Y{i} = outmesh_Y{i}+2;
end
outimg_w = ceil(all_maxx - all_minx) + 2 ;
outimg_h = ceil(all_maxy - all_miny) + 2 ;
%-------------

%% opencv
for i = 1 : 2 %img_n
    if num==1
        [c1out{i},c1omask{i}]=interpolat_location(inmesh_X{i},outmesh_X{i},inmesh_Y{i},outmesh_Y{i},in_name{i},outimg_w,outimg_h,i);
    else
        
        if i==1
            [c1out{i},c1omask{i},ble_mask]=interpolat(inmesh_X{i},outmesh_X{i},inmesh_Y{i},outmesh_Y{i},in_name{i},outimg_w,outimg_h,s_mask,b_mask);
        else
            num=num+1;
            [c1out{i},c1omask{i}]=interpolat_location(inmesh_X{i},outmesh_X{i},inmesh_Y{i},outmesh_Y{i},in_name{i},outimg_w,outimg_h,num);
        end
    end
    figure,imshow(c1out{i});
    title(num2str(i));
%---------------------
end

%--------------------
%不能使用作者方法进行融合
%--------------------
T4 = [1, 0, 1; 0, 1, 1; 0, 0, 1]; % translation by 1 pixels because matlab coordinate starts from 1

% T = T4*T3*T2*T1;
    
    