function [H,Index_final]=RANSAC(RE_COO)
%RANSAC�㷨���ڽǵ�ƥ����ȡ��Ӧ�Ծ���H��
%����re_cooΪ2*n*2ά��������re_coo(:,:,1)��ű�׼ͼ�������㣬re_coo(:,:,2)
%��ž�������ƥ����ģ����׼ͼ���Ӧ�Ĵ���׼ͼ�������㡣re_coo(1,:,:)��ź�����
%re_coo(2,:,:)��������ꡣ��ע��matlab��ij����ϵ��xy����ϵ��
%�м���T_ppΪ����RANSAC�㷨����ȡ��ƥ��ǵ�T_pp(1,:,:)��Ϊx����
%���HΪ��Ӧ�Ծ���T_PPΪ������RE_COO��ͬ���ݸ�ʽ������
%��RE_COO��T_pp��ʽ��ͬʱ��T_PPÿ����ά���ݱ�ΪT_pp��Ӧ���ݱ��ת�ã���ͬ�򲻱䡣
%H�������������ǰ���xy����ϵ��ȡ�ģ�����ij����ϵ��
% [m,n,v]=size(re_coo);
% C=zeros(4,n*(n-1)*(n-2)*(n-3)/(4*3*2*1));   %C��ÿһ�д��4��������������ֵ����СΪ4*Cn4
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
%     re_hg_coo2=[re_coo(:,:,2);ones(1,n)];%��Ϊ�������
%     re_hg_coo1=[re_coo(:,:,1);ones(1,n)];
%     %�������H��xoy����ϵ���Ӧ
%     T=H_inal*re_hg_coo2;                 %transform����������
%     T_fenmu=[T(3,:);T(3,:);T(3,:)];
%     T_minus=T./T_fenmu-re_hg_coo1;
%     T_dist=T_minus(1,:).^2+T_minus(2,:).^2;
%     T_index=find(T_dist<8);         %������ֵ
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
%     error('���㹻���ƥ���');
% else
%     T_pp=re_coo(:,Index_final,:);
% end
%����ΪRANSAC�㷨ʵ����ȷƥ��ǵ����ȡ
%%
%RANSAC_2
w=0.3;
NT=0;
k_diedai=0;
i=1;
i_loop=1;
% loopmin=0.05*( ceil( log(0.01)/( log(1-w^4) ) ) ); % --> Ϊ��������׼���ȣ�����Ĭ�ϵ����ٵ�������
% RE_COO=RE_COO(:,[2 1],:);
[m,n,v]=size(RE_COO);
%�ж�������RE_COO������ά����Ϣ����RE_COOΪn*d(nΪ���ݸ���)ʱ������û����ά����
%����ת�ã����RE_COOΪd*nʱ��������ת�á�14/10/28�ġ�
if m>n
    re_coo(:,:,1)=RE_COO(:,:,1)';
    re_coo(:,:,2)=RE_COO(:,:,2)';
    [m,n,v]=size(re_coo);
    roll_sign=1;                  %roll_signΪ��־����֤���������ʽ������һ��
else
    re_coo=RE_COO;
    roll_sign=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k_loop=ceil( log(0.01)/( log(1-w^4) ) );
while i<=k_loop 
    k=ceil( log(0.02)/( log(1-w^4) ) ); %���������ļ��㸴�Ƶ����档14/10/29�ġ�
    while i_loop==1;                 %Index��ֵӦ�ò�ͬ����whileѭ����������һ��
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
%��Index��ָ��ĵ��а���������ͬ���Ƿ���ͬ���������ʱ��rank(A)<8,�󲻳�
%�任�����е�8����������˶�A���Ƚ����жϣ�С��8���½���ȡֵ�� 14/10/28�ġ�
    if rank(A)<8    
        continue;
    else i=i+1;
    end
    B=[x11,y11,x21,y21,x31,y31,x41,y41]';
    G=(A'*A)\(A'*B);
    G=G';
    H_inal=[G(1:3);G(4:6);G(7:8),1];
    re_hg_coo2=[re_coo(:,:,2);ones(1,n)];%��Ϊ�������
    re_hg_coo1=[re_coo(:,:,1);ones(1,n)];
    %�������H��xoy����ϵ���Ӧ
    T=H_inal*re_hg_coo2;                 %transform����������
    T_fenmu=[T(3,:);T(3,:);T(3,:)];
    T_minus=T./T_fenmu-re_hg_coo1;
    T_dist=T_minus(1,:).^2+T_minus(2,:).^2;
    T_index=find(T_dist<0.1);         %������ֵ
    nt=length(T_index);
    if nt>NT
        NT=nt;
        Index_final=T_index;
    end
    w_1=nt/n;
    k_diedai=k_diedai+1;
    if w_1>w  
        w=w_1;
        k=ceil( log(0.02)/( log(1-w^4) ) ); %�˴�Ҳ����������������ļ��㡣14/10/29�ġ�
    end
    if k_diedai>=k 
        break;
    end
end
T_pp=re_coo(:,Index_final,:);
%%
%�����ó���
% figure;imshow([X1 X2]);hold on;
% n_x1=size(X1,2);
% plot(T_pp(1,:,1),T_pp(2,:,1),'r+',T_pp(1,:,2)+n_x1,T_pp(2,:,2),'y+');
%�����ó�������
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
% H=H(:,[2 1 3]);  %��xy����ϵ��Ϊij����ϵ��14/10/28�ġ�
if roll_sign==1
    T_PP(:,:,1)=T_pp(:,:,1)';
    T_PP(:,:,2)=T_pp(:,:,2)';
else
    T_PP=T_pp;
end



    
        
    
        
    
    
    