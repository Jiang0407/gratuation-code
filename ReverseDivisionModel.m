clear all;
clc; 
img_name{1}='E:/Images/0516ʵ����Ƶ/91.jpg';
img_name{2}='E:/Images/0516ʵ����Ƶ/92.jpg';
img_name{3}='E:/Images/0516ʵ����Ƶ/93.jpg';
img_name{4}='E:/Images/0516ʵ����Ƶ/94.jpg';
img = imread(img_name{2});
% imwrite(img(1:400,:,:),img_name{1});
%figure,imshow(img);
[m,n,cc] = size(img);%�õ�rgbͼƬ����������,����ߴ�,�����rgb��ccֵ3
mm = m;
nn = n;
C = zeros(mm,nn,cc);%��ԭͼ��ͬ����С��û������ֵ��0����

BS=200;
hBS=BS/2;
i0 = 237;%ͼ����������(293,359)
j0 = 362;
di = 232;%ÿ�����صĸ���
dj = 235;
lambda=0.32;

fprintf('�ٽ���ֵ��Ҫʱ�䣺'); tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=-hBS+1:mm+hBS
    for j=-hBS+1:nn+hBS
        x = (i-i0)/di;%(293,359)��ԭʼ����
        y = (j-j0)/dj;
        r=x^2+y^2;
        a=ceil(i0+(i-i0)*((1+4*lambda*r)^(1/2)-1)/(2*lambda*r));
        b=ceil(j0+(j-j0)*((1+4*lambda*r)^(1/2)-1)/(2*lambda*r)); 
        if a>0 && b>0 && a <= m && b<=n
            C(i+hBS,j+hBS,:)=img(a,b,:);
        end
    end
end
fprintf(' %fs\n',toc);
C = uint8(C);%��CתΪuint8��
figure,imshow(C);title('�ٽ���ֵ');
imwrite(C,'�ٽ���ֵ.jpg');
C=[];
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fprintf('˫���Բ�ֵ��Ҫʱ�䣺'); tic;
% for i=-hBS+1:mm+hBS
%     for j=-hBS+1:nn+hBS
%         x = (i-i0)/di;%(293,359)��ԭʼ����
%         y = (j-j0)/dj;
%         r=x^2+y^2;
%         a=(i0+(i-i0)*((1+4*lambda*r)^(1/2)-1)/(2*lambda*r));
%         b=(j0+(j-j0)*((1+4*lambda*r)^(1/2)-1)/(2*lambda*r)); 
%         a0=floor(a);a1=a0+1;sa0=a-a0;sa1=1-sa0;
%         b0=floor(b);b1=b0+1;sb0=b-b0;sb1=1-sb0;
%         if a0>0 && b0>0 && a1 <= m && b1<=n
%             tmp1=sa1*img(a0,b0,:)+sa0*img(a1,b0,:);
%             tmp2=sa1*img(a0,b1,:)+sa0*img(a1,b1,:);
%             C(i+hBS,j+hBS,:)=sb1*tmp1+sb0*tmp2;
%         end
%     end
% end
% fprintf(' %fs\n',toc);
% C = uint8(C);%��CתΪuint8��
% figure,imshow(C);title('˫���Բ�ֵ');
% imwrite(C,'˫���Բ�ֵ.jpg');