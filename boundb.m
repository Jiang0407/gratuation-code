function endb=boundb(img)
[~,n,~]=size(img);
bw=im2bw(img,0);
jfirst=zeros(1,n);
for j=1:n
    B=sum(bw(:,j));
    if B~=0        
        jfirst(j)=find(bw(:,j),1,'first');
    end
end
endb=find(jfirst~=0,1,'last');