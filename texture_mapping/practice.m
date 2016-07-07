close all;
outmeshX=round(load('outmeshX.txt'));
outmeshY=round(load('outmeshY.txt'));
inmeshX=load('inmeshX.txt');
inmeshY=load('inmeshY.txt');
img1=imread('4.jpg');
grid_size=2;
[m,n]=size(outmeshX);
cw=max(outmeshX(:,end));
ch=max(outmeshY(end,:));
outb_Lx=outmeshX(1,1);
outb_Uy=outmeshY(1,1);
warp_img=zeros(ch,cw,3);
% warp_img=zeros(500,1000,3);
A=[];
% for i=1:m-1
%     for j=1:n-1
%         outbx=outmeshX(i,j);
%         outby=outmeshY(i,j);
%         outspaceX=outmeshX(i,j+1)-outmeshX(i,j);
%         outspaceY=outmeshY(i+1,j)-outmeshY(i,j);
%         inimgpatch=img1( (1+(i-1)*grid_size):(i*grid_size), (1+(j-1)*grid_size): (j*grid_size), :);
%         RE_inimgpatch=imresize(inimgpatch,[outspaceY,outspaceX]);
%         warp_img(outby:outby+outspaceY-1,outbx:outbx+outspaceX-1,:)=RE_inimgpatch;
%         A=[A;outspaceX outspaceY];
% %         if sum(sum(warp_img(outby:outby+outspaceY-1,outbx:outbx+outspaceX-1,1)))==0
% %             [i j]
% %             pause;
% %         end
%     end
% end
%---------------------------
%以上为原始gridsize=10时的算法
%---------------------------
% for i=1:m-1
%     for j=1:n-1
%         outbx=outmeshX(i,j);
%         outby=outmeshY(i,j);
%         outspaceX=outmeshX(i,j+1)-outmeshX(i,j);
%         outspaceY=outmeshY(i+1,j)-outmeshY(i,j);
%         if outspaceX>0 && outspaceY>0
%         inimgpatch=img1( (1+(i-1)*grid_size):(i*grid_size), (1+(j-1)*grid_size): (j*grid_size), :);
%         RE_inimgpatch=imresize(inimgpatch,[outspaceY,outspaceX]);
%         warp_img(outby:outby+outspaceY-1,outbx:outbx+outspaceX-1,:)=RE_inimgpatch;
%         end
% %         if sum(sum(warp_img(outby:outby+outspaceY-1,outbx:outbx+outspaceX-1,1)))==0
% %             [i j]
% %             pause;
% %         end
%     end
% end
%--------------------
%最终的方法
%--------------------
for i=1:m-1
    for j=1:n-1
        outb_Lx(i,j)=outmeshX(i,j);
        outb_Uy(i,j)=outmeshY(i,j);
        outspaceX=outmeshX(i,j+1)-outmeshX(i,j);
        outspaceY=outmeshY(i+1,j)-outmeshY(i,j);
        inimgpatch=img1( (1+(i-1)*grid_size):(i*grid_size), (1+(j-1)*grid_size): (j*grid_size), :);
        RE_inimgpatch=imresize(inimgpatch,[outspaceY,outspaceX]);
        outb_Rx(i,j)=outb_Lx(i,j)+outspaceX-1;
        outb_Dy(i,j)=outb_Uy(i,j)+outspaceY-1;
        warp_img(outb_Uy(i,j):outb_Dy(i,j),outb_Lx(i,j):outb_Rx(i,j),:)=RE_inimgpatch;
%         if i>1 && j>1
%             if outb_Rx(i,j-1)-outb_Lx(i-1,j)>0 && outb_Dy(i-1,j-1)-outb_Uy(i,j)>0
%                 warp_img(outb_Uy(i,j),outb_Lx(i,j)-1,:)=warp_img(outb_Uy(i,j),outb_Lx(i,j),:);
%             end
%             if outb_Rx(i-1,j-1)-outb_Lx(i,j)>0 && outb_Dy(i-1,j)-outb_Uy(i,j-1)>0
%                 warp_img(outb_Uy(i,j)-1,outb_Lx(i,j)+1,:)=warp_img(outb_Uy(i,j),outb_Lx(i,j),:);
%             end
%         end
%         if sum(sum(warp_img(outby:outby+outspaceY-1,outbx:outbx+outspaceX-1,1)))==0
%             [i j]
%             pause;
%         end
    end
end
%---------------------
%继续插值，消除空洞
%---------------------
for k=1:ch
    inXf=find(warp_img(k,:,1)~=0,1);
    inXl=find(warp_img(k,:,1)~=0,1,'last');
    for h=inXf:inXl
        if warp_img(k,h,1)==0
            warp_img(k,h,:)=warp_img(k,h-1,:);
        end
    end
end
%--------------------
warp_img=uint8(warp_img);
figure;imshow(warp_img);
        
        
        
R=interp2(inmeshX,inmeshY,double(img1(:,:,1)),outmesh_X,outmeshY,'nearest');
G=interp2(inmeshX,inmeshY,double(img1(:,:,2)),outmesh_X,outmeshY,'nearest');
B=interp2(inmeshX,inmeshY,double(img1(:,:,3)),outmesh_X,outmeshY,'nearest');
R=round(R);
G=round(G);
B=round(B);
outimg1=cat(3,R,G,B);
figure;imshow(uint8(outimg1));
