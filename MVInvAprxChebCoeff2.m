function C=MVInvAprxChebCoeff2(hc,K,a1,a2,b1,b2)
[m,l]=size(hc);
%X=sym('X',l,1);
syms x y;
pp=0;
for i=1:m
    pp=pp+hc(i,l)*prod([x,y].^hc(i,1:l-1)');
end

h=matlabFunction(pp,'Vars',[x,y]);

f=@(x,y,k1,k2) cos(k1.*x).*cos(k2.*y)./h((a1+a2)/2+(a2-a1)/2*cos(x),(b1+b2)/2+(b2-b1)/2*cos(y));

o=0;
for i=1:K+1
    for j=1:K+1
            if i+j-2<=K
                [mm,pp]=size(find([i-1,j-1]==0));            
                c=2^(l-1-pp)/pi^(l-1)*integral2(@(x,y)f(x,y,i-1,j-1),0,pi,0,pi);
                o=o+1;
                C(o,1:3)=[i-1,j-1,c];
            end
    end
end
