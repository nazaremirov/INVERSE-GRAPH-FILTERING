function H=PolFilter_H(hc,MM)

[n,nn]=size(MM);
[m,l]=size(hc);

q=zeros(n,n);

for i=1:m
    CC=eye(n,n);
    for j=1:l-1
        CC=CC*sparse(MM(1:n,(j-1)*n+1:j*n))^hc(i,j);
    end
    q=q+hc(i,l)*CC;
end

H=q;