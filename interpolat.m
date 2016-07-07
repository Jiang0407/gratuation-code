function [warp_img,warp_mask,ble_mask]=interpolat( inmeshX, outmeshX, inmeshY, outmeshY,in_name,outimg_w, outimg_h, s_mask,b_mask)
outmeshX=round(outmeshX);
outmeshY=round(outmeshY);
inmeshX=round(inmeshX);
inmeshY=round(inmeshY);

[m,n]=size(outmeshX);
img1=imread(in_name);
warp_img=zeros(outimg_h,outimg_w,3);
warp_mask=zeros(outimg_h,outimg_w,3);%融合用模板
ble_mask=zeros(outimg_h,outimg_w,4);
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
        outb_Rx(i,j)=outmeshX(i,j)+outspaceX-1;
        outb_Dy(i,j)=outmeshY(i,j)+outspaceY-1; 
        if outspaceX>0 && outspaceY>0   
            if i<m-1 && j<n-1
                xa=inmeshX(i,j);
                xb=inmeshX(i,j+1)-1;
                ya=inmeshY(i,j);
                yb=inmeshY(i+1,j)-1;
            elseif i==m-1 && j<n-1
                xa=inmeshX(i,j);
                xb=inmeshX(i,j+1);
                ya=inmeshY(i,j);
                yb=inmeshY(i+1,j)-1;
            elseif i<m-1 && j==n-1
                xa=inmeshX(i,j);
                xb=inmeshX(i,j+1)-1;
                ya=inmeshY(i,j);
                yb=inmeshY(i+1,j);
            else 
                xa=inmeshX(i,j);
                xb=inmeshX(i,j+1);
                ya=inmeshY(i,j);
                yb=inmeshY(i+1,j);
            end
            patch_s=s_mask(ya:yb,xa:xb,:);
            ir_patchs=imresize(patch_s,[outspaceY,outspaceX],'nearest');
            warp_mask(outb_Uy(i,j):outb_Dy(i,j),outb_Lx(i,j):outb_Rx(i,j),:)=ir_patchs;
            
            patch_b1=b_mask(ya:yb,xa:xb,1);
            patch_b2=b_mask(ya:yb,xa:xb,2);
            patch_b3=b_mask(ya:yb,xa:xb,3);
            patch_b4=b_mask(ya:yb,xa:xb,4);
            ir_patchb1=imresize(patch_b1,[outspaceY,outspaceX],'nearest');
            ir_patchb2=imresize(patch_b2,[outspaceY,outspaceX],'nearest');
            ir_patchb3=imresize(patch_b3,[outspaceY,outspaceX],'nearest');
            ir_patchb4=imresize(patch_b4,[outspaceY,outspaceX],'nearest');
            ble_mask(outb_Uy(i,j):outb_Dy(i,j),outb_Lx(i,j):outb_Rx(i,j),1)=ir_patchb1;
            ble_mask(outb_Uy(i,j):outb_Dy(i,j),outb_Lx(i,j):outb_Rx(i,j),2)=ir_patchb2;
            ble_mask(outb_Uy(i,j):outb_Dy(i,j),outb_Lx(i,j):outb_Rx(i,j),3)=ir_patchb3;
            ble_mask(outb_Uy(i,j):outb_Dy(i,j),outb_Lx(i,j):outb_Rx(i,j),4)=ir_patchb4;
%------------------
%生成图像
%------------------
%             [i,j]
            patch_img=img1(ya:yb,xa:xb,:);
            ir_patchi=imresize(patch_img,[outspaceY,outspaceX],'nearest');
            warp_img(outb_Uy(i,j):outb_Dy(i,j),outb_Lx(i,j):outb_Rx(i,j),:)=ir_patchi;
        end
    end
end
%---------------------
%继续插值，消除空洞
%---------------------
for k=2:outimg_h %从2开始
    inXf=find(warp_mask(k,:,1)~=0,1);
    inXl=find(warp_mask(k,:,1)~=0,1,'last');
    for h=inXf+1:inXl-1
        if warp_mask(k,h,1)==0 && warp_mask(k-1,h,1)~=0
            warp_mask(k,h,:)=warp_mask(k,h-1,:); %表格空洞
            warp_img(k,h,:)=warp_img(k,h-1,:); %图像空洞
        end
    end
end
warp_img=uint8(warp_img);
% figure(1);imshow(warp_img);
