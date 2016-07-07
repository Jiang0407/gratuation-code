function output_canvas = imageblending1(warped_img1,warped_img2,w1,w2)
    
%     w1 = imfill(im2bw(uint8(warped_img1), 0),'holes');
%     w2 = imfill(im2bw(uint8(warped_img2), 0),'holes');
    w1=im2bw(w1,0);
    w2=im2bw(w2,0);

%     w2 = mat2gray(w2,[0.1,1]);
% 
    warped_img1 = double(warped_img1);
    warped_img2 = double(warped_img2);
%%
%直接平均法
    output_canvas(:,:,1) = ((warped_img1(:,:,1).*w1)+(warped_img2(:,:,1).*w2))./(w1+w2);
    output_canvas(:,:,2) = ((warped_img1(:,:,2).*w1)+(warped_img2(:,:,2).*w2))./(w1+w2);
    output_canvas(:,:,3) = ((warped_img1(:,:,3).*w1)+(warped_img2(:,:,3).*w2))./(w1+w2);
    output_canvas = uint8(output_canvas);
figure;imshow(output_canvas);
title('改进的SPHP算法配准效果图');
%%
%线性融合法
[m,n]=size(w1);
output_canvas=zeros(m,n,3);
for i=1:m
    for j=1:n
        w=w1(i,j)+w2(i,j);
        if w==1
            if w1(i,j)==1
                output_canvas(i,j,:)=warped_img1(i,j,:);
            else
                output_canvas(i,j,:)=warped_img2(i,j,:);
            end
        elseif w==2
            j1last=find(w1(i,:)==1,1,'last');
            j2first=find(w2(i,:)==1,1,'first');
            wht1=(j1last-j)/(j1last-j2first);
            wht2=1-wht1;
            output_canvas(i,j,:)=wht1*warped_img1(i,j,:)+wht2*warped_img2(i,j,:);
        else continue;
        end
    end
end
output_canvas=uint8(output_canvas);
figure;imshow(output_canvas);

            
        
    
    