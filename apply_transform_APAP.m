function [X0,Y0]=apply_transform_APAP(X, Y, Hmdlt)
[m,n]=size(X);
X0=zeros(m,n);
Y0=zeros(m,n);
for i=1:m
    for j=1:n
        Hap=Hmdlt{i,j};
        X0(i,j)=( Hap(1,1)*X(i,j)+Hap(1,2)*Y(i,j)+Hap(1,3) )/...
                ( Hap(3,1)*X(i,j)+Hap(3,2)*Y(i,j)+Hap(3,3) );
        Y0(i,j)=( Hap(2,1)*X(i,j)+Hap(2,2)*Y(i,j)+Hap(2,3) )/...
                ( Hap(3,1)*X(i,j)+Hap(3,2)*Y(i,j)+Hap(3,3) );     
    end
end
A=0;
        