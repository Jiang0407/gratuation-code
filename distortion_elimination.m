function data_out= distortion_elimination(data)

%�������������Ҫ�޸Ĳ������������������Ҫ���棬��Ϊ���ڻ���ʹ��
%��һ��ͼ���������
i0 = 267; 
j0 = 362;
di = 232; 
dj = 235;
lambda=0.25;
 
[~,n,~]=size(data);
data_out = ones(size(data));

for k=1:n
    i=data(2,k); %��
    j=data(1,k); %��
    x = (i-i0)/di;%(293,359)��ԭʼ����
    y = (j-j0)/dj;
    r=x^2+y^2;
    data_out(2,k) = (i0 + (i-i0) * (1-lambda*r)^(1/2));
    data_out(1,k) = (j0 + (j-j0) * (1-lambda*r)^(1/2));
end

%�ڶ���ͼ���������
i0 = 237; 
j0 = 362;
di = 232; 
dj = 235;
lambda=0.32; 
for k=1:n
    i=data(5,k); %��
    j=data(4,k); %��
    x = (i-i0)/di;%(293,359)��ԭʼ����
    y = (j-j0)/dj;
    r=x^2+y^2;
    data_out(5,k) = ceil(i0 + (i-i0) * (1-lambda*r)^(1/2));
    data_out(4,k) = ceil(j0 + (j-j0) * (1-lambda*r)^(1/2));
end

for k=1:n
    if(data_out(1,k)>600 || data_out(1,k)<0 || data_out(4,k)>600 || data_out(4,k)<0 || data_out(2,k)>800 || data_out(2,k)<0 || data_out(5,k)>800 || data_out(5,k)<0)
        data_out(3,k)=0;
    end
end
data_out1=[];
for k=1:n
    if(data_out(3,k) ~= 0)
        data_out1(:,k)=data_out(:,k);
    end
end
data_out=data_out1;