Ax=InmeshX{1}(95:100,98:100);
Bx=InmeshY{1}(95:100,98:100);
for i=1:6
    for j=1:3
        H=Hmdlt{i+94,j+97};
        x11=Ax(i,j);
        y11=Bx(i,j);
        Bwx(i,j)=(H(1, 1)*x11 + H(1, 2)*y11 + H(1, 3)) ./ ...
                (H(3, 1)*x11 + H(3, 2)*y11 + H(3, 3));
        Bwy(i,j) = (H(2, 1)*x11 + H(2, 2)*y11 + H(2, 3)) ./ ...
                (H(3, 1)*x11 + H(3, 2)*y11 + H(3, 3));
    end
end

for i=1:100
    for j=90:100
        H=re_Hmdlt{i,j};
        AAA{i,j-89}=H-re_Hmdlt{1,100};
    end
end