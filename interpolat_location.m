function [warp_img,warp_mask]=interpolat_location( inmeshX, outmeshX, inmeshY, outmeshY, in_name,outimg_w, outimg_h,num)
outmeshX=round(outmeshX);
outmeshY=round(outmeshY);
inmeshX=round(inmeshX);
inmeshY=round(inmeshY);

img1=imread(in_name); %不需要函数输入in_name
[m,n]=size(outmeshX);

%----------------------
warp_img=uint8( zeros(outimg_h,outimg_w,3) );
warp_mask=zeros(outimg_h,outimg_w,3);%融合用模板,第三维用来指示坐标来自哪个图像
%----------------------

outb_Lx=outmeshX(1,1);
outb_Uy=outmeshY(1,1);

for i=1:m-1
    for j=1:n-1
        outb_Lx(i,j)=outmeshX(i,j);
        outb_Uy(i,j)=outmeshY(i,j);
        
        if i<m-1 && j<n-1
            outspaceX=outmeshX(i,j+1)-outmeshX(i,j);
            outspaceY=outmeshY(i+1,j)-outmeshY(i,j);
        elseif i==m-1 && j<n-1
            outspaceX=outmeshX(i,j+1)-outmeshX(i,j);
            outspaceY=outmeshY(i+1,j)-outmeshY(i,j)+1; 
        elseif i<m-1 && j==n-1
            outspaceX=outmeshX(i,j+1)-outmeshX(i,j)+1;
            outspaceY=outmeshY(i+1,j)-outmeshY(i,j);
        else
        outspaceX=outmeshX(i,j+1)-outmeshX(i,j)+1;
        outspaceY=outmeshY(i+1,j)-outmeshY(i,j)+1;
        end     
        
        if outspaceX>0 && outspaceY>0   
            if i<m-1 && j<n-1
                [inpatch_x,inpatch_y]=meshgrid(inmeshX(i,j):inmeshX(i,j+1)-1,inmeshY(i,j):inmeshY(i+1,j)-1);
            elseif i==m-1 && j<n-1
                [inpatch_x,inpatch_y]=meshgrid(inmeshX(i,j):inmeshX(i,j+1),inmeshY(i,j):inmeshY(i+1,j)-1);
            elseif i<m-1 && j==n-1
                [inpatch_x,inpatch_y]=meshgrid(inmeshX(i,j):inmeshX(i,j+1)-1,inmeshY(i,j):inmeshY(i+1,j));
            else 
                [inpatch_x,inpatch_y]=meshgrid(inmeshX(i,j):inmeshX(i,j+1),inmeshY(i,j):inmeshY(i+1,j));
            end
            
            outpatch_x=round( imresize(inpatch_x,[outspaceY,outspaceX]) );
            outpatch_y=round( imresize(inpatch_y,[outspaceY,outspaceX]) );
            
            outb_Rx(i,j)=outb_Lx(i,j)+outspaceX-1;
            outb_Dy(i,j)=outb_Uy(i,j)+outspaceY-1; 
            
%---------------------------
%MATLAB第一维数组为i坐标，第二维为j坐标
%---------------------------
            warp_mask(outb_Uy(i,j):outb_Dy(i,j),outb_Lx(i,j):outb_Rx(i,j),1)=outpatch_y; %第一维为i坐标
            warp_mask(outb_Uy(i,j):outb_Dy(i,j),outb_Lx(i,j):outb_Rx(i,j),2)=outpatch_x; %第二维为j坐标
            warp_mask(outb_Uy(i,j):outb_Dy(i,j),outb_Lx(i,j):outb_Rx(i,j),3)=num*ones( (outb_Dy(i,j)-outb_Uy(i,j)+1),outb_Rx(i,j)-outb_Lx(i,j)+1); %第三位标记
        end
    end
end

%---------------------
%继续插值，消除位置空洞
%---------------------
for k=2:outimg_h %从2开始
    inXf=find(warp_mask(k,:,1)~=0,1);
    inXl=find(warp_mask(k,:,1)~=0,1,'last');
    for h=inXf+1:inXl-1
        if warp_mask(k,h,1)==0 && warp_mask(k-1,h,1)~=0
            warp_mask(k,h,:)=warp_mask(k,h-1,:);
        end
    end
end 
% % %--------------------------
% % %查找鱼眼图像中的像素点的位置
% % %--------------------------
% % i0 = 267; 
% % j0 = 362;
% % di = 232; 
% % dj = 235;
% % lambda=0.25;
% % for i=1:outimg_h
% %     for j=1:outimg_w
% %         if warp_mask(i,j,1)~=0
% %                 x = (warp_mask(i,j,1)-i0)/di;%(293,359)是原始坐标
% %                 y = (warp_mask(i,j,2)-j0)/dj;
% %                 r=x^2+y^2;
% %                 warp_mask(i,j,1)=ceil(i0+(warp_mask(i,j,1)-i0)*((1+4*lambda*r)^(1/2)-1)/(2*lambda*r));
% %                 warp_mask(i,j,2)=ceil(j0+(warp_mask(i,j,2)-j0)*((1+4*lambda*r)^(1/2)-1)/(2*lambda*r)); 
% %         else
% %             warp_img(i,j,:)=uint8(zeros(1,1,3));
% %         end
% %     end
% % end
%--------------------------
%图像查表显示图像
%--------------------------
% fprintf('查表所需时间： ');tic;
for i=1:outimg_h
    for j=1:outimg_w
        if warp_mask(i,j,1)~=0 && warp_mask(i,j,1)>0 &&warp_mask(i,j,2)>0
            a=warp_mask(i,j,1);
            b=warp_mask(i,j,2);
            warp_img(i,j,:)=img1(a,b,:);
        else
            warp_img(i,j,:)=uint8(zeros(1,1,3));
        end
    end
end
% fprintf('%fs\n',toc);
% figure;imshow(warp_img);
% title(num2str(num)); 

            
        