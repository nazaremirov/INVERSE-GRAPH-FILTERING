function G=ChebFilter_G(ch,MM,K,a)

[mm,l]=size(ch);

NC=[];
tk=1;
for i=1:mm
    if sum(ch(i,1:l-1))<=K
        NC(tk,:)=ch(i,:);
        tk=tk+1;
    end
end

[m,l]=size(NC);

[n,nn]=size(MM);


T=[];

for i=1:nn/n
    S{i}=(2/(a(2*i)-a(2*i-1)))*MM(1:n,(i-1)*n+1:i*n)-((a(2*i-1)+a(2*i))/(a(2*i)-a(2*i-1)))*sparse(eye(n,n));
    T=[T;ChebShift(S{i},K)];
end



q=zeros(n,n);

for i=1:m
    Q=eye(n,n);
    for j=1:l-1
        Q=Q*sparse(T((j-1)*n+1:j*n,NC(i,j)*n+1:(NC(i,j)+1)*n));
    end
    q=q+(NC(i,l))*Q;
end
G=q;