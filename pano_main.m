close all;
clear all;
clc;

run('vlfeat-0.9.14/toolbox/vl_setup');
num_img=2;
img_name=cell(1,num_img);
% img_name{1}='images/case26/new_img/411/sta/zuo_01.bmp';
% img_name{2}='images/case26/new_img/411/sta/you_01.bmp';
% img_name{1}='images/case26/new_img/video2/aleft_01.bmp';
% img_name{2}='images/case26/new_img/video2/aright_01.bmp';
img_name{1}='images/case26/new_img/411/grass/WP_001.jpg';
img_name{2}='images/case26/new_img/411/grass/WP_002.jpg';
img_name{3}='images/case26/new_img/411/grass/WP_003.jpg';
img_name{4}='images/case26/new_img/411/grass/WP_004.jpg';
% img_name{5}='images/case26/new_img/411/grass/WP_005.jpg';
% img_name{6}='images/case26/new_img/411/grass/WP_006.jpg';
% img_name{7}='images/case26/new_img/411/grass/WP_007.jpg';
% img_name{8}='images/case26/new_img/411/grass/WP_008.jpg';
% img_name{9}='images/case26/new_img/411/grass/WP_009.jpg'; 
% img_name{5}='images/case26/new_img/411/grass/WP_005.jpg';
% img_name{1}='images/case26/new_img/411/machinery_building/WPL_001.jpg';
% img_name{2}='images/case26/new_img/411/machinery_building/WPL_002.jpg';
% img_name{1}='images/case26/new_img/411/grass/WP_001.jpg';
% img_name{2}='images/case26/new_img/411/grass/WP_002.jpg';
% img_name{5}='images/case26/new_img/411/grass/WP_005.jpg';
% img_name{1}='images/case26/new_img/411/main_building/WP_001.jpg';
% img_name{2}='images/case26/new_img/411/main_building/WP_002.jpg';
% img_name{3}='images/case26/new_img/411/main_building/WP_003.jpg';
% img_name{4}='images/case26/new_img/411/main_building/WP_004.jpg';
% img_name{5}='images/case26/new_img/411/main_building/WP_005.jpg';
% img_name{6}='images/case26/new_img/411/main_building/WP_006.jpg';
% img_name{7}='images/case26/new_img/411/main_building/WP_007.jpg';
img_name{1}='E:/Images/new/1.jpg';
img_name{2}='E:/Images/new/2.jpg'; 
img_name{3}='E:/Images/new/3.jpg';
img_name{4}='E:/Images/new/4.jpg'; 
img_name{5}='E:/Images/new/5.jpg';
img_name{6}='E:/Images/new/6.jpg'; 
img_name{1}='E:/Images/1#/1.jpg';
img_name{2}='E:/Images/1#/2.jpg'; 
img_name{3}='E:/Images/1#/3.jpg';
img_name{4}='E:/Images/1#/4.jpg'; 
img_name{5}='E:/Images/1#/5.jpg';
img_name{6}='E:/Images/1#/6.jpg'; 
img_name{1}='E:/Images/0525/3.jpg';
img_name{2}='E:/Images/0525/4.jpg';
img_name{3}='E:/Images/0525/3.jpg';
img_name{4}='E:/Images/0525/4.jpg';
img_name{5}='E:/Images/0525/5.jpg';
img_name{6}='E:/Images/0525/6.jpg'; 

img_name{1}='F:/桌面杂/0704/videoL/imgL1.bmp';
img_name{2}='F:/桌面杂/0704/videoR/imgR1.bmp';
img_name{1}='F:/桌面杂/0706-2/imageL/imgL150.bmp';
img_name{2}='F:/桌面杂/0706-2/imageR/imgR150.bmp';
img_name{1}='F:/桌面杂/0707/imageL/imgL100.bmp';
img_name{2}='F:/桌面杂/0707/imageR/imgR100.bmp';
% img_name{1}='F:/桌面杂/0706/l.jpg';
% img_name{2}='F:/桌面杂/0706/r.jpg';
% img_name{1}='C:/Users/Administrator/Desktop/0702/img771.bmp';
% img_name{2}='C:/Users/Administrator/Desktop/0702/ret772.bmp';
% img_name{3}='E:/Images/new/9.jpg';
% img_name{4}='E:/Images/new/4.jpg'; 
% img_name{1}='E:/Images/0516实验视频/91.jpg';
% img_name{2}='E:/Images/0516实验视频/92.jpg';
% img_name{3}='E:/程序/code/images/0516_1 (3).jpg';
% img_name{4}='E:/程序/code/images/0516_1 (4).jpg';

% 0516_4.jpg
% %-----------------------------------0
% img_name{1}='images/case26/new_img/411/SIX_SYSTEM/b1.jpg';
% img_name{2}='images/case26/new_img/411/SIX_SYSTEM/b2.jpg';
% img_name{3}='images/case26/new_img/411/SIX_SYSTEM/b3.jpg';
% img_name{4}='images/case26/new_img/411/SIX_SYSTEM/b4.jpg';
% img_name{5}='images/case26/new_img/411/SIX_SYSTEM/b5.jpg';
% img_name{6}='images/case26/new_img/411/SIX_SYSTEM/b6.jpg';
% 
% for i=1:num_img
%     tmp=imread(img_name{i});
%     [m n ~]=size(tmp);
%     if m*n>800*600
%         imwrite(imresize(tmp,0.3),img_name{i});
%     end
% end
in_name = cell(2,1); 
tic
for i=1:num_img-1
    if i==1
        in_name{1} = img_name{i};
        in_name{2} = img_name{i+1}; 
        [warp_img,sti_mask,ble_mask]=SPHP_stitching_mine(in_name,i);
        med_img=warp_img;
        %imwrite(med_img,'images/case26/new_img/411/main_building/WPL_100.bmp');
        figure,imshow(med_img)
        pause
    else
        med_name='images/case26/new_img/411/main_building/WPL_100.bmp';
        s_mask=sti_mask;
        b_mask=ble_mask;
         in_name{1} = med_name;
         in_name{2} = img_name{i+1}; 
        [warp_img,sti_mask,ble_mask]=SPHP_stitching_mine(in_name,i,s_mask,b_mask);
        imwrite(warp_img,'images/case26/new_img/411/main_building/WPL_100.bmp');
        %center0=center00; 
        med_img=warp_img;  
        figure,imshow(med_img)
        pause
    end
end
toc;
figure;imshow(warp_img); 
imwrite(warp_img,'main_b5.jpg')
save BBX.mat
% sti_mask=uint16(sti_mask);
%%
stiz=sti_mask(:,:,1);
stiy=sti_mask(:,:,2);
stix=sti_mask(:,:,3);
stii=stiz';
stij=stiy';
stil=stix';

% ble_mask=uint16(ble_mask);
ble11=ble_mask(:,:,1);
ble12=ble_mask(:,:,2);
ble21=ble_mask(:,:,3);
ble22=ble_mask(:,:,4);
ble1i=ble11';
ble1j=ble12';
ble2i=ble21';
ble2j=ble22';
save BBXX.mat
% % --------------------
% % % 验证制表法时间
% % % --------------------
% cd SPHP_code11;
% mex ../test2.cpp;
% load BBX.mat
% sti=sti_mask';
% ble=ble_mask';
% st=sti_mask(1:100,1);
% bl=ble_mask(1:100,1);
% createxml('st_table',sti_mask(1:100,1:100,1));
% createxml('bl_table',ble_mask(1:100,1:100,1));
%%
img=cell(num_img,1);
for i=1:num_img
    img{i}=imread( img_name{i} );
end
    
[m,n,~]=size(warp_img)
IMG=uint8(zeros(m,n,3));
% mex test2.cpp -ID:\E\OPENCV\opencv\build\include -LD:\E\OPENCV\opencv\build\x64\vc12\lib;
% test2(sti_mask,ble_mask);

fprintf( '查表所需时间： ');tic;
for i=1:m
    for j=1:n
        switch sti_mask(i,j,3)
            case 1
                IMG(i,j,:)=img{1}(sti_mask(i,j,1),sti_mask(i,j,2),:);
            case 12
                IMG(i,j,:)=sti_mask(i,j,1)*img{1}(ble_mask(i,j,1),ble_mask(i,j,2),:)+sti_mask(i,j,2)*img{2}(ble_mask(i,j,3),ble_mask(i,j,4),:);
            case 2
                IMG(i,j,:)=img{2}(sti_mask(i,j,1),sti_mask(i,j,2),:);
            case 23
                IMG(i,j,:)=sti_mask(i,j,1)*img{2}(ble_mask(i,j,1),ble_mask(i,j,2),:)+sti_mask(i,j,2)*img{3}(ble_mask(i,j,3),ble_mask(i,j,4),:);
            case 3
                IMG(i,j,:)=img{3}(sti_mask(i,j,1),sti_mask(i,j,2),:);
            case 34
                IMG(i,j,:)=sti_mask(i,j,1)*img{3}(ble_mask(i,j,1),ble_mask(i,j,2),:)+sti_mask(i,j,2)*img{4}(ble_mask(i,j,3),ble_mask(i,j,4),:);
            case 4
                IMG(i,j,:)=img{4}(sti_mask(i,j,1),sti_mask(i,j,2),:);
            case 45
                IMG(i,j,:)=sti_mask(i,j,1)*img{4}(ble_mask(i,j,1),ble_mask(i,j,2),:)+sti_mask(i,j,2)*img{5}(ble_mask(i,j,3),ble_mask(i,j,4),:); 
            case 5
                IMG(i,j,:)=img{5}(sti_mask(i,j,1),sti_mask(i,j,2),:);
            case 56
                IMG(i,j,:)=sti_mask(i,j,1)*img{5}(ble_mask(i,j,1),ble_mask(i,j,2),:)+sti_mask(i,j,2)*img{6}(ble_mask(i,j,3),ble_mask(i,j,4),:); 
            case 6
                IMG(i,j,:)=img{6}(sti_mask(i,j,1),sti_mask(i,j,2),:);
        end
    end
end
fprintf('%fs\n',toc);
figure;imshow(IMG);

%------------------
%向量化方法
%------------------
% for j=1:n
%     jsum=sum(sti_mask(:,j,1));
%     jf=find(jsum~=0,1,'first');
%     if size(jf)~=0
%         aj0(j)=1;
%     else
%         aj0(j)=0;
%     end
% end
% jfir=find(aj0==1,1,'first');
% jlas=find(aj0==1,1,'last');
% mstim=sti_mask(:,jfir:jlas,:);
% mblem=ble_mask(:,jfir:jlas,:);
% for i=1:m
%     i0=find(mstim(i,:,3)==0);
%     if size(i0,2)==0
%         ai0(i)=1;
%     else
%         ai0(i)=0;
%     end
% end
% ifir=find(ai0==1,1,'first');
% ilas=find(ai0==1,1,'last');
% fstim=mstim(ifir:ilas,:,:);
% fblem=mblem(ifir:ilas,:,:);
%     
% bound=zeros(2,num_img-1,1);
% [~,a]=find(fstim(:,:,3)==12);
% amax=max(a);
% amin=min(a);
% bound(1,1)=amin;
% bound(2,1)=amax;
% [~,b]=find(fstim(:,:,3)==23);
% bmax=max(b);
% bmin=min(b);
% bound(1,2)=bmin;
% bound(2,2)=bmax;
% [~,c]=find(fstim(:,:,3)==34);
% cmax=max(c);
% cmin=min(c);
% bound(1,3)=cmin;
% bound(2,3)=cmax;
% [~,d]=find(fstim(:,:,3)==45);
% dmax=max(d);
% dmin=min(d);
% bound(1,4)=dmin;
% bound(2,4)=dmax;
% aa=fstim(:,1:bound(1,1)-1,1);
% ab=fstim(:,1:bound(1,1)-1,2);
% IMG(:,1:bound(1,1),:)=img{1}(aa(:),ab(:),:);
        

    






