function Number = numsetup(nz,nx)
%setup numbering system
Number = zeros(nz,nx);
num = 1;
for j=1:nx    
for i=1:nz
Number(i,j) = num;
num = num+1;
end
end
end