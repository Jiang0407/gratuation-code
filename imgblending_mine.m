function [stitched_img,reg_mask,b_mask]=imgblending_mine(IMG,way,blendsize,endb,c1omask,num,I,ble_mask)
img1=IMG{1};
img2=IMG{2};
[m,n,~]=size(img1);
stitched_img=zeros(m,n,3);
if ~( isa(img1,'uint8') && isa(img2,'uint8') )
    error('Both input images should be with Class Uint8');
end  
 %------------------------
 %以下用于局部融合
 %------------------------
 
%         endb=find(j1first~=0,1,'last');
        beginj=endb-blendsize+1;
        reg_mask=zeros(m,n,3);
%         b_mask=cell(1,2);
%         b_mask{1}=zeros(m,n,2); %用于存储融合区域坐标。第一维为i坐标，第二维为j坐标
%         b_mask{2}=zeros(m,n,2);
          b_mask=zeros(m,n,4);        
        img1(:,endb+1:end,:)=zeros(m,n-endb,3);
        img2(:,1:beginj-1,:)=zeros(m,beginj-1,3);
        figure,imshow(img1)
        figure,imshow(img2)
        pause
        I{1}=uint8(I{1}*255);
        I{2}=uint8(I{2}*255);

 %------------------
 %用于制表
 %------------------
        reg_mask(:,1:beginj-1,:)=c1omask{1}(:,1:beginj-1,:);
        reg_mask(:,endb+1:end,:)=c1omask{2}(:,endb+1:end,:);
        if nargin==8
            b_mask(:,1:beginj-1,:)=ble_mask(:,1:beginj-1,:);
            b_mask(:,endb+1:end,:)=ble_mask(:,endb+1:end,:);
        end
%-------------------        
%         c1omask{1}(:,endb+1:end,:)=zeros(m,n-endb,3);
%         c1omask{2}(:,1:beginj-1,:)=zeros(m,beginj-1,3);
%         blend_mask1=c1omask{1}(:,beginj:endb,:);
%         blend_mask2=c1omask{2}(:,beginj:endb,:);
        
        bw1=im2bw(img1,0);
        bw2=im2bw(img2,0);
        for i=1:m
            A1=sum(bw1(i,:));
            A2=sum(bw2(i,:));
            if A1~=0     
                i1first(i)=find(bw1(i,:),1,'first');
                i1last(i)=find(bw1(i,:),1,'last');
            else
                i1first(i)=0;
                i1last(i)=0;
            end
            if A2~=0
                i2first(i)=find(bw2(i,:),1,'first');
                i2last(i)=find(bw2(i,:),1,'last');
            else
                i2first(i)=0;
                i2last(i)=0;
            end
        end

        for j=1:n
            B1=sum(bw1(:,j));
            B2=sum(bw2(:,j));
            if B1~=0        
                j1first(j)=find(bw1(:,j),1,'first');
                j1last(j)=find(bw1(:,j),1,'last');
            else
                j1first(j)=0;
                j1last(j)=0;
            end
            if B2~=0
                j2first(j)=find(bw2(:,j),1,'first');
                j2last(j)=find(bw2(:,j),1,'last');
            else
                j2first(j)=0;
                j2last(j)=0;  
            end
        end        
%%      
        fprintf('实时融合所需时间： ');tic;
        for i=1:m
            for j=1:beginj-1
                stitched_img(i,j,:)=img1(i,j,:);
            end
            for j=endb+1:n
                stitched_img(i,j,:)=img2(i,j,:);
            end
        end
%         stitched_img(:,1:beginj-1,:)=img1(:,1:beginj-1,:);
%         stitched_img(:,endb+1:end,:)=img2(:,endb+1:end,:);
        for i=1:m
            for j=beginj:endb
                if bw1(i,j)==1 && bw2(i,j)==1
                    di1=min( abs(i-j1first(j)), abs(i-j1last(j)) );
                    dj1=min( abs(j-i1first(i)), abs(j-i1last(i)) );
                    d1=min(di1,dj1);
                    di2=min( abs(i-j2first(j)), abs(i-j2last(j)) );
                    dj2=min( abs(j-i2first(i)), abs(j-i2last(i)) );
                    d2=min(di2,dj2);
                    if d1==0 && d2==0
                        w1=1;%范围由0-1变为0-100的整数.5/16改
                    else
                        w1=d1/(d1+d2);%范围由0-1变为0-100的整数.5/16改
                    end
                    w2=1-w1;
%                     w1=(floor(w1*10))/10;
%                     w2=(floor(w2*10))/10;
                    stitched_img(i,j,:)=(w1*img1(i,j,:)+w2*img2(i,j,:));%像素值回到0-255.5/16改。
 %------------------
 %查表用
 %------------------
                    reg_mask(i,j,3)=num*10+num+1;
                    reg_mask(i,j,1)=w1;
                    reg_mask(i,j,2)=w2;
%                     b_mask{1}(i,j,:)=c1omask{1}(i,j,1:2);
%                     b_mask{2}(i,j,:)=c1omask{2}(i,j,1:2);
                    b_mask(i,j,1:2)=c1omask{1}(i,j,1:2);
                    b_mask(i,j,3:4)=c1omask{2}(i,j,1:2);
                elseif bw1(i,j)==1
                    stitched_img(i,j,:)=img1(i,j,:);
                    reg_mask(i,j,:)=c1omask{1}(i,j,:); %查表用
                else
                    stitched_img(i,j,:)=img2(i,j,:);
                    reg_mask(i,j,:)=c1omask{2}(i,j,:); %查表用
                end
            end
        end
        fprintf('%fs\n',toc);
%         figure;imshow(uint8(stitched_img));
%         title('融合图像')
 
 %--------------------
 %验证查表结果
 %只有两幅图像时使用
 %--------------------
%  s_img=zeros(m,n,3);
%  fprintf('查表融合所需时间： ');tic
%  for i=1:m
%      for j=1:n
%          if reg_mask(i,j,3)~=0
%              a=reg_mask(i,j,1);
%              b=reg_mask(i,j,2);
%              c=reg_mask(i,j,3);
%              switch c
%                  case num
%                      s_img(i,j,:)=I{1}(a,b,:);
%                  case num+1
%                      s_img(i,j,:)=I{2}(a,b,:);
%                  case num*10+num+1
%                      p1=b_mask{1}(i,j,:);
%                      p2=b_mask{2}(i,j,:);
%                      s_img(i,j,:)=a*I{1}(p1(1),p1(2),:)+...
%                          b*I{2}(p2(1),p2(2),:); %改到这4.27
%              end
%          end
%      end
%  end
%  fprintf('%fs\n',toc);
%   figure(10);imshow(uint8(s_img)); 
end

