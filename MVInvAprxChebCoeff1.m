function C=MVInvAprxChebCoeff1(hc,K,a1,a2)
syms x;
[m,l]=size(hc);
p=hc(m);
for i=1:m-1
    p=p+hc(i)*x^(m-i);
end
h=matlabFunction(p);
f=@(x,i) cos((i-1).*x)./h((a1+a2)/2+(a2-a1)/2*cos(x));
for i=1:K+1
    C(i)=2/pi*integral(@(x)f(x,i),0,pi);
end 
C(1)=C(1)/2;