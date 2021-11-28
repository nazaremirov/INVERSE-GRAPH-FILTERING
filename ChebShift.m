function T=ChebShift(S,k)

[n,n]=size(S);
T(1:n,1:n)=eye(n,n);

if k>=1
T(1:n,n+1:2*n)=S;
end

if k>=2
    
for i=3:k+1
    T(1:n,(i-1)*n+1:i*n)=2*S*T(1:n,(i-2)*n+1:(i-1)*n)-T(1:n,(i-3)*n+1:(i-2)*n);
end

end
