function [H,Index_final]=RANSAC(RE_COO)
%RANSAC算法用于角点匹配求取单应性矩阵H。
%输入re_coo为2*n*2维矩阵，其中re_coo(:,:,1)存放标准图像特征点，re_coo(:,:,2)
%存放经过初步匹配过的，与标准图像对应的待配准图像特征点。re_coo(1,:,:)存放横坐标
%re_coo(2,:,:)存放纵坐标。（注意matlab中ij坐标系和xy坐标系）
%中间量T_pp为经过RANSAC算法所提取的匹配角点T_pp(1,:,:)内为x坐标
%输出H为单应性矩阵，T_PP为与输入RE_COO相同数据格式的坐标
%当RE_COO与T_pp格式不同时，T_PP每个二维数据表为T_pp对应数据表的转置，相同则不变。
%H矩阵求解过程中是按照xy坐标系求取的，不是ij坐标系。
% [m,n,v]=size(re_coo);
% C=zeros(4,n*(n-1)*(n-2)*(n-3)/(4*3*2*1));   %C中每一列存放4个特征点列坐标值。大小为4*Cn4
% i0=1;j0=2;k0=3;l0=4;
% q=1;
% for i=i0:n-3
%     for j=j0:n-2
%         for k=k0:n-1
%             for l=l0:n
%                 C(:,q)=[i;j;k;l];
%                 q=q+1;
%             end
%             l0=l0+1;
%         end
%         k0=k0+1;
%         l0=k0+1;
%     end
%     j0=j0+1;
%     k0=j0+1;
%     l0=k0+1;
% end
% [mc,nc]=size(C);
% NT=0;
% for i_m=1:nc
%     D=re_coo(:,(C(:,i_m))',:);
%     x11=D(1,1,1);
%     x21=D(1,2,1);
%     x31=D(1,3,1);
%     x41=D(1,4,1);
%     y11=D(2,1,1);
%     y21=D(2,2,1);
%     y31=D(2,3,1);
%     y41=D(2,4,1);
%     
%     x12=D(1,1,2);
%     x22=D(1,2,2);
%     x32=D(1,3,2);
%     x42=D(1,4,2);
%     
%     y12=D(2,1,2);
%     y22=D(2,2,2);
%     y32=D(2,3,2);
%     y42=D(2,4,2);
%     
%     A1=[x12,0,x22,0,x32,0,x42,0;...
%         y12,0,y22,0,y32,0,y42,0;...
%         1,0,1,0,1,0,1,0;...
%         0,x12,0,x22,0,x32,0,x42;...
%         0,y12,0,y22,0,y32,0,y42;...
%         0,1,0,1,0,1,0,1];
%     A2=-[x11*x12,y11*x12,x21*x22,y21*x22,x31*x32,y31*x32,x41*x42,y41*x42;...
%         x11*y12,y11*y12,x21*y22,y21*y22,x31*y32,y31*y32,x41*y42,y41*y42];
%     A=[A1',A2'];
%     B=[x11,y11,x21,y21,x31,y31,x41,y41]';
%     G=(A'*A)\(A'*B);
%     G=G';
%     H_inal=[G(1:3);G(4:6);G(7:8),1];
%     re_hg_coo2=[re_coo(:,:,2);ones(1,n)];%变为齐次坐标
%     re_hg_coo1=[re_coo(:,:,1);ones(1,n)];
%     %求出来的H与xoy坐标系相对应
%     T=H_inal*re_hg_coo2;                 %transform后的齐次坐标
%     T_fenmu=[T(3,:);T(3,:);T(3,:)];
%     T_minus=T./T_fenmu-re_hg_coo1;
%     T_dist=T_minus(1,:).^2+T_minus(2,:).^2;
%     T_index=find(T_dist<8);         %设置阈值
%     if ~isempty(T_index)
%         nt=length(T_index);
%         if nt>NT
%            NT=nt;
%            Index_final=T_index;
%            H=H_inal;
%         end
%     end        
% %     else T_sum(i)=0;
% %     end 
% end
% if NT<4
%     error('无足够多的匹配点');
% else
%     T_pp=re_coo(:,Index_final,:);
% end
%以上为RANSAC算法实现正确匹配角点的提取
%%
%RANSAC_2
w=0.3;
NT=0;
k_diedai=0;
i=1;
i_loop=1;
% loopmin=0.05*( ceil( log(0.01)/( log(1-w^4) ) ) ); % --> 为了增加配准精度，增加默认的最少迭代次数
% RE_COO=RE_COO(:,[2 1],:);
[m,n,v]=size(RE_COO);
%判定输入量RE_COO的数据维数信息，当RE_COO为n*d(n为数据个数)时，对其没个二维数据
%进行转置，如果RE_COO为d*n时，不进行转置。14/10/28改。
if m>n
    re_coo(:,:,1)=RE_COO(:,:,1)';
    re_coo(:,:,2)=RE_COO(:,:,2)';
    [m,n,v]=size(re_coo);
    roll_sign=1;                  %roll_sign为标志，保证输出数据形式和输入一样
else
    re_coo=RE_COO;
    roll_sign=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k_loop=ceil( log(0.01)/( log(1-w^4) ) );
while i<=k_loop 
    k=ceil( log(0.02)/( log(1-w^4) ) ); %抽样次数的计算复制到下面。14/10/29改。
    while i_loop==1;                 %Index的值应该不同，用while循环来满足这一点
        Index=round(rand(1,4)*(n-1))+1;
        if Index(1)==Index(2) || Index(1)==Index(3) || Index(1)==Index(4)...
                 || Index(2)==Index(3) || Index(2)==Index(4)...
                 || Index(3)==Index(4)
             continue;
        else break;
        end
    end
    
    D=re_coo(:,Index,:);
    x11=D(1,1,1);
    x21=D(1,2,1);
    x31=D(1,3,1);
    x41=D(1,4,1);
    y11=D(2,1,1);
    y21=D(2,2,1);
    y31=D(2,3,1);
    y41=D(2,4,1);
    
    x12=D(1,1,2);
    x22=D(1,2,2);
    x32=D(1,3,2);
    x42=D(1,4,2);
    
    y12=D(2,1,2);
    y22=D(2,2,2);
    y32=D(2,3,2);
    y42=D(2,4,2);
    
    A1=[x12,0,x22,0,x32,0,x42,0;...
        y12,0,y22,0,y32,0,y42,0;...
        1,0,1,0,1,0,1,0;...
        0,x12,0,x22,0,x32,0,x42;...
        0,y12,0,y22,0,y32,0,y42;...
        0,1,0,1,0,1,0,1];
    A2=-[x11*x12,y11*x12,x21*x22,y21*x22,x31*x32,y31*x32,x41*x42,y41*x42;...
        x11*y12,y11*y12,x21*y22,y21*y22,x31*y32,y31*y32,x41*y42,y41*y42];
    A=[A1',A2'];
%当Index所指向的点中包含坐标相同但是方向不同的特征点的时候，rank(A)<8,求不出
%变换矩阵中的8个参数，因此对A的秩进行判断，小于8重新进行取值。 14/10/28改。
    if rank(A)<8    
        continue;
    else i=i+1;
    end
    B=[x11,y11,x21,y21,x31,y31,x41,y41]';
    G=(A'*A)\(A'*B);
    G=G';
    H_inal=[G(1:3);G(4:6);G(7:8),1];
    re_hg_coo2=[re_coo(:,:,2);ones(1,n)];%变为齐次坐标
    re_hg_coo1=[re_coo(:,:,1);ones(1,n)];
    %求出来的H与xoy坐标系相对应
    T=H_inal*re_hg_coo2;                 %transform后的齐次坐标
    T_fenmu=[T(3,:);T(3,:);T(3,:)];
    T_minus=T./T_fenmu-re_hg_coo1;
    T_dist=T_minus(1,:).^2+T_minus(2,:).^2;
    T_index=find(T_dist<0.1);         %设置阈值
    nt=length(T_index);
    if nt>NT
        NT=nt;
        Index_final=T_index;
    end
    w_1=nt/n;
    k_diedai=k_diedai+1;
    if w_1>w  
        w=w_1;
        k=ceil( log(0.02)/( log(1-w^4) ) ); %此处也进行所需迭代次数的计算。14/10/29改。
    end
    if k_diedai>=k 
        break;
    end
end
T_pp=re_coo(:,Index_final,:);
%%
%测试用程序
% figure;imshow([X1 X2]);hold on;
% n_x1=size(X1,2);
% plot(T_pp(1,:,1),T_pp(2,:,1),'r+',T_pp(1,:,2)+n_x1,T_pp(2,:,2),'y+');
%测试用程序终了
%%
n_index_final=size(Index_final,2);
AA=[];
BB=[];
for kk=1:n_index_final
    x_1_1=T_pp(1,kk,1);
    y_1_1=T_pp(2,kk,1);   
    x_1_2=T_pp(1,kk,2); 
    y_1_2=T_pp(2,kk,2);
    AA=[AA;x_1_2,y_1_2,1,0,0,0,-x_1_1*x_1_2,-x_1_1*y_1_2;
        0,0,0,x_1_2,y_1_2,1,-x_1_2*y_1_1,-y_1_1*y_1_2];
    BB=[BB;x_1_1;y_1_1];
end
h=(AA'*AA)\(AA'*BB);
h=h';
H=[h(1:3);h(4:6);[h(7:8),1]];
% H=H(:,[2 1 3]);  %由xy坐标系改为ij坐标系。14/10/28改。
if roll_sign==1
    T_PP(:,:,1)=T_pp(:,:,1)';
    T_PP(:,:,2)=T_pp(:,:,2)';
else
    T_PP=T_pp;
end



    
        
    
        
    
    
    