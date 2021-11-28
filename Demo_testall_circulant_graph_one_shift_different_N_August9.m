% Main Data we need from this simulation are
% p : number of iteration
% table{NN} for each NN: is an average relative error E(m) for m=1,2,...,p for each NodeNumber(i)
% NumberNode : is a values for N, ex: NumberNode=[50,100,..., 2000]
% Table_Iter :
% TABLE_Time :
% But if you can you need to save complete workspace
% for plotting use different matlab code "Circulant_Numeric_Conv_calculation_plot_Aug31.m" that I attached (it gives all three plots that we may use)

% please uncomment the line 1015 "%%%save('Revision_CirculantGraph_Different_N_trial_1000_until_N2000')" at the end of the code to save the complete workspace

clear;clc;format short
A_L=[];
B_K=[];
Conv_rate_icpa0={};
Conv_rate_icpa1={};
Conv_rate_icpa2={};
Conv_rate_icpa3={};
Conv_rate_icpa4={};
Conv_rate_icpa5={};
Conv_rate_iopa0={};
Conv_rate_iopa1={};
Conv_rate_iopa2={};
Conv_rate_iopa3={};
Conv_rate_iopa4={};
Conv_rate_iopa5={};
Conv_rate_gradient={};
Conv_rate_arma={};
table={};

epsilon=0.001;

% Change this parameters (take largr numbers)
NumberNode=[50,100,200,400,600,800,1000,1200,1400,1600,1800,2000];
%NumberNode=[250,500,1000,2000,4000,8000,16000,32000,64000,128000];

for NN=1:length(NumberNode)
    
S1=[];S2=[];S3=[];S=[];L_C=[];H=[];c=[];Shifted_L=[];T=[];a_L=[];b_K=[];

N=NumberNode(NN);

s=[1,2,5]';  % Generating set

[l,l1]=size(s);

Z=[];

for u=1:l
TT=zeros(N,N);
for i=1:N
    for j=1:N
        if mod(i-j,N)==s(u)
           TT(i,j)=1;
        end
    end
end
Z(1:N,1+N*(u-1):N*u)=TT+TT';
end

S1=Z(1:N,1:N);  % Shift matrix C_1
S2=Z(1:N,N+1:2*N);  % Shift matrix C_2
S3=Z(1:N,2*N+1:3*N);  % Shift matrix C_5

S=S1+S2+S3; 

%--------------Normalized Laplacian--------------------------

L_C=eye(N,N)-diag(sum(S').^(-1/2))*S*diag(sum(S').^(-1/2));

arma_conv_rate_theor(NN)=norm(L_C,2)*4/9;

L_C=sparse(L_C);



%% -----------------------Parameters--------------------------

K=5; %maximum approximation power (until K)

L=5; % Order of an approx. polynomial

p=30; %number of iteration

h=[-1,-3/4,27/4]'; % filter coefficients

a=[0,2]; % spectrum bound


%% ---------- Filter H --------------

H=h(3)*eye(N,N)+h(2)*L_C+h(1)*L_C^2;

H=sparse(H);

%% -----------------Inverse for Centralized---------------

G_cent=inv(H);

%% -------------Chebyshev Coefficient and Filter G_K----------
  
c=MVInvAprxChebCoeff1(h,K,0,2);


Shifted_L=(2/(a(2)-a(1)))*L_C-((a(1)+a(2))/(a(2)-a(1)))*sparse(eye(N,N));

T=ChebShift(Shifted_L,K);
  
q=c(1)*T(1:N,1:N);
G_K{1}=q;

b_K(1)=max(abs(eig(eye(N,N)-H*G_K{1})));

for k=2:K+1
    q=q+(c(k))*T(1:N,N*(k-1)+1:N*k);
    G_K{k}=q; % G{k} is a cheb. pol. approx. of order k-1
    b_K(k)=max(abs(eig(eye(N,N)-H*G_K{k})));
end


%% ---------------------- Optimal polynomial approx. G_L ------------------

eig_L=eig(L_C); 

eig_matrix= fliplr(vander(eig_L));

E = h(3)*ones(size(eig_L))+h(2).*eig_L+h(1).*eig_L.*eig_L;

matrix_A_temp= diag(E)*eig_matrix; 

G_L={};

for jj=1:L+1

matrix_A= matrix_A_temp(:, 1:jj); 

A_eq=[-ones(N, 1)  matrix_A;
     -ones(N, 1)  -matrix_A]; 
 
b_eq=[ones(N, 1)
      -ones(N, 1)];
  
b_eq=b_eq';
  
coe_f=[1 zeros(1, jj)];

s_op=linprog(coe_f,A_eq,b_eq); % s_op contains the largest eigenvalue of G and 
                               % the coefficient for GL 
GL=[];
                               
GL=s_op(2).*eye(size(L_C));

for ll=3:length(s_op)
    GL=GL+ s_op(ll).*(L_C^(ll-2));
end

G_L{jj}=GL;

a_L(jj)=max(abs(eig(eye(N,N)-H*G_L{jj})));
end

gamma=2/(min(abs(eig(H)))+max(abs(eig(H))));



%% ---------- Start the simulation for 1000 trial -----------------
trial=1000;

sum_error_icpa0=zeros(1,p);
sum_error_icpa1=zeros(1,p);
sum_error_icpa2=zeros(1,p);
sum_error_icpa3=zeros(1,p);
sum_error_icpa4=zeros(1,p);
sum_error_icpa5=zeros(1,p);

sum_error_iopa0=zeros(1,p);
sum_error_iopa1=zeros(1,p);
sum_error_iopa2=zeros(1,p);
sum_error_iopa3=zeros(1,p);
sum_error_iopa4=zeros(1,p);
sum_error_iopa5=zeros(1,p);

sum_error_arma=zeros(1,p);
sum_error_gradient=zeros(1,p);

sum_X_icpa0=zeros(N,p);
sum_X_icpa1=zeros(N,p);
sum_X_icpa2=zeros(N,p);
sum_X_icpa3=zeros(N,p);
sum_X_icpa4=zeros(N,p);
sum_X_icpa5=zeros(N,p);

sum_X_iopa0=zeros(N,p);
sum_X_iopa1=zeros(N,p);
sum_X_iopa2=zeros(N,p);
sum_X_iopa3=zeros(N,p);
sum_X_iopa4=zeros(N,p);
sum_X_iopa5=zeros(N,p);

sum_X_arma=zeros(N,p);
sum_X_gradient=zeros(N,p);

for ii=1:trial

signal_o=2*rand(N, 1)-1;

denom_signal=norm(signal_o,2);

observ= H* signal_o;     
          
dsp3=['trial=', num2str(ii)]; disp(dsp3)

%% ---------------Centralized --------------------------

tic
x_cent=G_cent*observ;
time_cent{NN}{ii}=toc;

%% --------------------ARMA------------------------------

ob_ini_arma1=observ;
ob_ini_arma2=observ;

x_initial1=zeros(N,1);
x_initial2=zeros(N,1);

X_arma=[];
Error_arma=[] ;

final_iter_arma{NN}{ii}=0;

time_arma{ii}=[];

for m=1:p
    
   tic
    
   x_initial1= (4/9) * L_C * x_initial1 + observ;
%    
   x_initial2= (-1/3) * L_C * x_initial2+ observ;
   
   x_arma=(16/189).*x_initial1+(4/63).*x_initial2; 
   
   time_arma{ii}(m)=toc;
   
   X_arma=[X_arma x_arma];
   
   Error_arma=[Error_arma (norm(signal_o-x_arma,2)/denom_signal)];
   
   if (norm(signal_o-x_arma,2)/denom_signal)>epsilon
       %Error_arma=[Error_arma (norm(signal_o-x_arma,2)/denom_signal)];
       final_iter_arma{NN}{ii}=final_iter_arma{NN}{ii}+1;
   else
       %Error_arma=[Error_arma 0];
       continue
   end
      
end



%% ---------------------ICPA-----------------------------

bb0=observ;

x_icpa0=sparse(zeros(N,1));
z0=sparse(zeros(N,1));

X_icpa0=[];
Error_icpa0=[];
final_iter_icpa0{NN}{ii}=0;
time_icpa0{ii}=[];

for m=1:p
    tic

    z0=sparse(G_K{1}*bb0);
    x_icpa0=sparse(x_icpa0+z0);
    bb0=sparse(bb0-H*z0);
    
    time_icpa0{ii}(m)=toc;
    
    X_icpa0=[X_icpa0 x_icpa0];
    
    Error_icpa0=[Error_icpa0 (norm(signal_o-x_icpa0,2)/denom_signal)];
    
   if (norm(signal_o-x_icpa0,2)/denom_signal)>epsilon
       
       %Error_icpa0=[Error_icpa0 (norm(signal_o-x_icpa0,2)/denom_signal)];
       final_iter_icpa0{NN}{ii}=final_iter_icpa0{NN}{ii}+1;
   else
       %Error_icpa0=[Error_icpa0 0];
       continue;
   end
    

end


bb1=observ;

x_icpa1=sparse(zeros(N,1));
z1=sparse(zeros(N,1));

X_icpa1=[];
Error_icpa1=[];

final_iter_icpa1{NN}{ii}=0;
time_icpa1{ii}=[];

for m=1:p
    tic
    
    z1=sparse(G_K{2}*bb1);
    x_icpa1=sparse(x_icpa1+z1);
    bb1=sparse(bb1-H*z1);
    
    time_icpa1{ii}(m)=toc;
    
    X_icpa1=[X_icpa1 x_icpa1];
    
    Error_icpa1=[Error_icpa1 (norm(signal_o-x_icpa1,2)/denom_signal)];

   if (norm(signal_o-x_icpa1,2)/denom_signal)>epsilon
       
       %Error_icpa1=[Error_icpa1 (norm(signal_o-x_icpa1,2)/denom_signal)];
       final_iter_icpa1{NN}{ii}=final_iter_icpa1{NN}{ii}+1;
   else
       %Error_icpa1=[Error_icpa1 0];
       continue;
   end
end

bb2=observ;

x_icpa2=sparse(zeros(N,1));
z2=sparse(zeros(N,1));

X_icpa2=[];
Error_icpa2=[];

final_iter_icpa2{NN}{ii}=0;
time_icpa2{ii}=[];

for m=1:p
    tic
    
    z2=sparse(G_K{3}*bb2);
    x_icpa2=sparse(x_icpa2+z2);
    bb2=sparse(bb2-H*z2);
    
    time_icpa2{ii}(m)=toc;
    
    X_icpa2=[X_icpa2 x_icpa2];
    
    Error_icpa2=[Error_icpa2 (norm(signal_o-x_icpa2,2)/denom_signal)];

   if (norm(signal_o-x_icpa2,2)/denom_signal)>epsilon
       
       %Error_icpa2=[Error_icpa2 (norm(signal_o-x_icpa2,2)/denom_signal)];
       final_iter_icpa2{NN}{ii}=final_iter_icpa2{NN}{ii}+1;
   else
       %Error_icpa2=[Error_icpa2 0];
       continue;
   end
end

bb3=observ;

x_icpa3=sparse(zeros(N,1));
z3=sparse(zeros(N,1));

X_icpa3=[];
Error_icpa3=[];

final_iter_icpa3{NN}{ii}=0;
time_icpa3{ii}=[];

for m=1:p
    tic
    
    z3=sparse(G_K{4}*bb3);
    x_icpa3=sparse(x_icpa3+z3);
    bb3=sparse(bb3-H*z3);
    
    time_icpa3{ii}(m)=toc;
    
    X_icpa3=[X_icpa3 x_icpa3];
    
    Error_icpa3=[Error_icpa3 (norm(signal_o-x_icpa3,2)/denom_signal)];

   if (norm(signal_o-x_icpa3,2)/denom_signal)>epsilon
       
       %Error_icpa3=[Error_icpa3 (norm(signal_o-x_icpa3,2)/denom_signal)];
       final_iter_icpa3{NN}{ii}=final_iter_icpa3{NN}{ii}+1;
   else
       %Error_icpa3=[Error_icpa3 0];
       continue;
   end
end

bb4=observ;

x_icpa4=sparse(zeros(N,1));
z4=sparse(zeros(N,1));

X_icpa4=[];
Error_icpa4=[];

final_iter_icpa4{NN}{ii}=0;
time_icpa4{ii}=[];

for m=1:p
    tic
    
    z4=sparse(G_K{5}*bb4);
    x_icpa4=sparse(x_icpa4+z4);
    bb4=sparse(bb4-H*z4);
    
    time_icpa4{ii}(m)=toc;
    
    X_icpa4=[X_icpa4 x_icpa4];
    
    Error_icpa4=[Error_icpa4 (norm(signal_o-x_icpa4,2)/denom_signal)];

   if (norm(signal_o-x_icpa4,2)/denom_signal)>epsilon
       
       %Error_icpa4=[Error_icpa4 (norm(signal_o-x_icpa4,2)/denom_signal)];
       final_iter_icpa4{NN}{ii}=final_iter_icpa4{NN}{ii}+1;
   else
       %Error_icpa4=[Error_icpa4 0];
       continue;
   end
end

bb5=observ;

x_icpa5=sparse(zeros(N,1));
z5=sparse(zeros(N,1));

X_icpa5=[];
Error_icpa5=[];

final_iter_icpa5{NN}{ii}=0;
time_icpa5{ii}=[];

for m=1:p
    tic
    
    z5=sparse(G_K{6}*bb5);
    x_icpa5=sparse(x_icpa5+z5);
    bb5=sparse(bb5-H*z5);
    
    time_icpa5{ii}(m)=toc;
    
    X_icpa5=[X_icpa5 x_icpa5];
    
    Error_icpa5=[Error_icpa5 (norm(signal_o-x_icpa5,2)/denom_signal)];

   if (norm(signal_o-x_icpa5,2)/denom_signal)>epsilon
       
       %Error_icpa5=[Error_icpa5 (norm(signal_o-x_icpa5,2)/denom_signal)];
       final_iter_icpa5{NN}{ii}=final_iter_icpa5{NN}{ii}+1;
   else
       %Error_icpa5=[Error_icpa5 0];
       continue;
   end
end
%% --------------------IOPA--------------------------------

x_iopa0=zeros(N,1);
zz0=zeros(N,1);

bbb0=observ;

X_iopa0=[];
Error_iopa0=[];

final_iter_iopa0{NN}{ii}=0;
time_iopa0{ii}=[];

for m=1:p
    tic
    
    zz0=sparse(G_L{1}*bbb0);
    x_iopa0=sparse(x_iopa0+zz0);
    bbb0=sparse(bbb0-H*zz0);
    
    time_iopa0{ii}(m)=toc;
    
    X_iopa0=[X_iopa0 x_iopa0];
    
    Error_iopa0=[Error_iopa0 (norm(signal_o-x_iopa0,2)/denom_signal)];

   if (norm(signal_o-x_iopa0,2)/denom_signal)>epsilon
       
       %Error_iopa0=[Error_iopa0 (norm(signal_o-x_iopa0,2)/denom_signal)];
       final_iter_iopa0{NN}{ii}=final_iter_iopa0{NN}{ii}+1;
   else
       %Error_iopa0=[Error_iopa0 0];
       continue
   end    
end


x_iopa1=zeros(N,1);
zz1=zeros(N,1);

bbb1=observ;
X_iopa1=[];
Error_iopa1=[];

final_iter_iopa1{NN}{ii}=0;
time_iopa1{ii}=[];

for m=1:p
    tic
    
    zz1=sparse(G_L{2}*bbb1);
    x_iopa1=sparse(x_iopa1+zz1);
    bbb1=sparse(bbb1-H*zz1);
    
    time_iopa1{ii}(m)=toc;
    
    X_iopa1=[X_iopa1 x_iopa1];
    
    Error_iopa1=[Error_iopa1 (norm(signal_o-x_iopa1,2)/denom_signal)];

   if (norm(signal_o-x_iopa1,2)/denom_signal)>epsilon
       
       %Error_iopa1=[Error_iopa1 (norm(signal_o-x_iopa1,2)/denom_signal)];
       final_iter_iopa1{NN}{ii}=final_iter_iopa1{NN}{ii}+1;
   else
       %Error_iopa1=[Error_iopa1 0];
       continue;
   end      
end


x_iopa2=zeros(N,1);
zz2=zeros(N,1);

bbb2=observ;
X_iopa2=[];
Error_iopa2=[];
time_iopa2{ii}=[];

final_iter_iopa2{NN}{ii}=0;

for m=1:p
    tic
    
    zz2=sparse(G_L{3}*bbb2);
    x_iopa2=sparse(x_iopa2+zz2);
    bbb2=sparse(bbb2-H*zz2);
    
    time_iopa2{ii}(m)=toc;
    
    X_iopa2=[X_iopa2 x_iopa2];
    
    Error_iopa2=[Error_iopa2 (norm(signal_o-x_iopa2,2)/denom_signal)];

   if (norm(signal_o-x_iopa2,2)/denom_signal)>epsilon
       
       %Error_iopa2=[Error_iopa2 (norm(signal_o-x_iopa2,2)/denom_signal)];
       final_iter_iopa2{NN}{ii}=final_iter_iopa2{NN}{ii}+1;
   else
       %Error_iopa2=[Error_iopa2 0];
       continue;
   end    
end

x_iopa3=zeros(N,1);
zz3=zeros(N,1);

bbb3=observ;
X_iopa3=[];
Error_iopa3=[];

final_iter_iopa3{NN}{ii}=0;
time_iopa3{ii}=[];

for m=1:p
    
    tic
    
    zz3=sparse(G_L{4}*bbb3);
    x_iopa3=sparse(x_iopa3+zz3);
    bbb3=sparse(bbb3-H*zz3);
    
    time_iopa3{ii}(m)=toc;
    
    X_iopa3=[X_iopa3 x_iopa3];
    
    Error_iopa3=[Error_iopa3 (norm(signal_o-x_iopa3,2)/denom_signal)];

   if (norm(signal_o-x_iopa3,2)/denom_signal)>epsilon
       
       %Error_iopa3=[Error_iopa3 (norm(signal_o-x_iopa3,2)/denom_signal)];
       final_iter_iopa3{NN}{ii}=final_iter_iopa3{NN}{ii}+1;
   else
       %Error_iopa3=[Error_iopa3 0];
       continue;
   end     
end


x_iopa4=zeros(N,1);
zz4=zeros(N,1);

bbb4=observ;
X_iopa4=[];
Error_iopa4=[];

final_iter_iopa4{NN}{ii}=0;
time_iopa4{ii}=[];

for m=1:p
    tic
    
    zz4=sparse(G_L{5}*bbb4);
    x_iopa4=sparse(x_iopa4+zz4);
    bbb4=sparse(bbb4-H*zz4);
    
    time_iopa4{ii}(m)=toc;
    
    X_iopa4=[X_iopa4 x_iopa4];
    
    Error_iopa4=[Error_iopa4 (norm(signal_o-x_iopa4,2)/denom_signal)];

   if (norm(signal_o-x_iopa4,2)/denom_signal)>epsilon
       
       %Error_iopa4=[Error_iopa4 (norm(signal_o-x_iopa4,2)/denom_signal)];
       final_iter_iopa4{NN}{ii}=final_iter_iopa4{NN}{ii}+1;
   else
       %Error_iopa4=[Error_iopa4 0];
       continue;
   end     
end


x_iopa5=zeros(N,1);
zz5=zeros(N,1);

bbb5=observ;
X_iopa5=[];
Error_iopa5=[];

final_iter_iopa5{NN}{ii}=0;
time_iopa5{ii}=[];

for m=1:p
    tic
    
    zz5=sparse(G_L{6}*bbb5);
    x_iopa5=sparse(x_iopa5+zz5);
    bbb5=sparse(bbb5-H*zz5);
    
    time_iopa5{ii}(m)=toc;
    
    X_iopa5=[X_iopa5 x_iopa5];
    
    Error_iopa5=[Error_iopa5 (norm(signal_o-x_iopa5,2)/denom_signal)];

   if (norm(signal_o-x_iopa5,2)/denom_signal)>epsilon
       
       %Error_iopa5=[Error_iopa5 (norm(signal_o-x_iopa5,2)/denom_signal)];
       final_iter_iopa5{NN}{ii}=final_iter_iopa5{NN}{ii}+1;
   else
       %Error_iopa5=[Error_iopa5 0];
       continue;
   end     
end

%% -----------------Gradient Descent--------------

y_gradient=zeros(N,1);
X_gradient=[];
Error_gradient=[];

final_iter_gradient{NN}{ii}=0;
time_gradient{ii}=[];

for m=1:p
    tic
    
    y_gradient=y_gradient-gamma*(H*y_gradient-observ);
    
    time_gradient{ii}(m)=toc;
    
    X_gradient=[X_gradient y_gradient];
    
    Error_gradient=[Error_gradient (norm(signal_o-y_gradient,2)/denom_signal)];
    
   if (norm(signal_o-y_gradient,2)/denom_signal)>epsilon
       
       %Error_gradient=[Error_gradient (norm(signal_o-y_gradient,2)/denom_signal)];
       final_iter_gradient{NN}{ii}=final_iter_gradient{NN}{ii}+1;
   else
       %Error_gradient=[Error_gradient 0];
       continue;
   end
end
%% ------------ Centralized  Methdod --------------

sum_error_icpa0=sum_error_icpa0+Error_icpa0;
sum_error_icpa1=sum_error_icpa1+Error_icpa1;
sum_error_icpa2=sum_error_icpa2+Error_icpa2;
sum_error_icpa3=sum_error_icpa3+Error_icpa3;
sum_error_icpa4=sum_error_icpa4+Error_icpa4;
sum_error_icpa5=sum_error_icpa5+Error_icpa5;

sum_error_iopa0=sum_error_iopa0+Error_iopa0;
sum_error_iopa1=sum_error_iopa1+Error_iopa1;
sum_error_iopa2=sum_error_iopa2+Error_iopa2;
sum_error_iopa3=sum_error_iopa3+Error_iopa3;
sum_error_iopa4=sum_error_iopa4+Error_iopa4;
sum_error_iopa5=sum_error_iopa5+Error_iopa5;

sum_error_arma=sum_error_arma+Error_arma;
sum_error_gradient=sum_error_gradient+Error_gradient;


sum_X_icpa0=sum_X_icpa0+X_icpa0;
sum_X_icpa1=sum_X_icpa1+X_icpa1;
sum_X_icpa2=sum_X_icpa2+X_icpa2;
sum_X_icpa3=sum_X_icpa3+X_icpa3;
sum_X_icpa4=sum_X_icpa4+X_icpa4;
sum_X_icpa5=sum_X_icpa5+X_icpa5;

sum_X_iopa0=sum_X_iopa0+X_iopa0;
sum_X_iopa1=sum_X_iopa1+X_iopa1;
sum_X_iopa2=sum_X_iopa2+X_iopa2;
sum_X_iopa3=sum_X_iopa3+X_iopa3;
sum_X_iopa4=sum_X_iopa4+X_iopa4;
sum_X_iopa5=sum_X_iopa5+X_iopa5;

sum_X_arma=sum_X_arma+X_arma;
sum_X_gradient=sum_X_gradient+X_gradient;

end

avg_error_icpa0=sum_error_icpa0/trial;
avg_error_icpa1=sum_error_icpa1/trial;
avg_error_icpa2=sum_error_icpa2/trial;
avg_error_icpa3=sum_error_icpa3/trial;
avg_error_icpa4=sum_error_icpa4/trial;
avg_error_icpa5=sum_error_icpa5/trial;

avg_error_iopa0=sum_error_iopa0/trial;
avg_error_iopa1=sum_error_iopa1/trial;
avg_error_iopa2=sum_error_iopa2/trial;
avg_error_iopa3=sum_error_iopa3/trial;
avg_error_iopa4=sum_error_iopa4/trial;
avg_error_iopa5=sum_error_iopa5/trial;

avg_error_arma=sum_error_arma/trial;
avg_error_gradient=sum_error_gradient/trial;

table{NN}=[avg_error_icpa0',avg_error_arma',avg_error_icpa1',avg_error_gradient',avg_error_iopa0',avg_error_icpa2',avg_error_iopa1',avg_error_icpa3',avg_error_icpa4',avg_error_iopa2',avg_error_icpa5',avg_error_iopa3',avg_error_iopa4',avg_error_iopa5'];


avg_X_icpa0=sum_X_icpa0/trial;
avg_X_icpa1=sum_X_icpa1/trial;
avg_X_icpa2=sum_X_icpa2/trial;
avg_X_icpa3=sum_X_icpa3/trial;
avg_X_icpa4=sum_X_icpa4/trial;
avg_X_icpa5=sum_X_icpa5/trial;

avg_X_iopa0=sum_X_iopa0/trial;
avg_X_iopa1=sum_X_iopa1/trial;
avg_X_iopa2=sum_X_iopa2/trial;
avg_X_iopa3=sum_X_iopa3/trial;
avg_X_iopa4=sum_X_iopa4/trial;
avg_X_iopa5=sum_X_iopa5/trial;

avg_X_arma=sum_X_arma/trial;
avg_X_gradient=sum_X_gradient/trial;


A_L=[A_L;a_L];
B_K=[B_K;b_K];

%Conv_rate_icpa0{NN}=avg_error_icpa0(2:final_iter_icpa0+1)./avg_error_icpa0(1:final_iter_icpa0);
Conv_rate_icpa1{NN}=avg_error_icpa1(2:max([final_iter_icpa1{NN}{1:trial}])+1)./avg_error_icpa1(1:max([final_iter_icpa1{NN}{1:trial}]));
Conv_rate_icpa2{NN}=avg_error_icpa2(2:max([final_iter_icpa2{NN}{1:trial}])+1)./avg_error_icpa2(1:max([final_iter_icpa2{NN}{1:trial}]));
Conv_rate_icpa3{NN}=avg_error_icpa3(2:max([final_iter_icpa3{NN}{1:trial}])+1)./avg_error_icpa3(1:max([final_iter_icpa3{NN}{1:trial}]));
Conv_rate_icpa4{NN}=avg_error_icpa4(2:max([final_iter_icpa4{NN}{1:trial}])+1)./avg_error_icpa4(1:max([final_iter_icpa4{NN}{1:trial}]));
Conv_rate_icpa5{NN}=avg_error_icpa5(2:max([final_iter_icpa5{NN}{1:trial}])+1)./avg_error_icpa5(1:max([final_iter_icpa5{NN}{1:trial}]));

Conv_rate_iopa0{NN}=avg_error_iopa0(2:max([final_iter_iopa0{NN}{1:trial}])+1)./avg_error_iopa0(1:max([final_iter_iopa0{NN}{1:trial}]));
Conv_rate_iopa1{NN}=avg_error_iopa1(2:max([final_iter_iopa1{NN}{1:trial}])+1)./avg_error_iopa1(1:max([final_iter_iopa1{NN}{1:trial}]));
Conv_rate_iopa2{NN}=avg_error_iopa2(2:max([final_iter_iopa2{NN}{1:trial}])+1)./avg_error_iopa2(1:max([final_iter_iopa2{NN}{1:trial}]));
Conv_rate_iopa3{NN}=avg_error_iopa3(2:max([final_iter_iopa3{NN}{1:trial}])+1)./avg_error_iopa3(1:max([final_iter_iopa3{NN}{1:trial}]));
Conv_rate_iopa4{NN}=avg_error_iopa4(2:max([final_iter_iopa4{NN}{1:trial}])+1)./avg_error_iopa4(1:max([final_iter_iopa4{NN}{1:trial}]));
Conv_rate_iopa5{NN}=avg_error_iopa5(2:max([final_iter_iopa5{NN}{1:trial}])+1)./avg_error_iopa5(1:max([final_iter_iopa5{NN}{1:trial}]));

Conv_rate_gradient{NN}=avg_error_gradient(2:max([final_iter_gradient{NN}{1:trial}])+1)./avg_error_gradient(1:max([final_iter_gradient{NN}{1:trial}]));

Conv_rate_arma{NN}=avg_error_arma(2:max([final_iter_arma{NN}{1:trial}])+1)./avg_error_arma(1:max([final_iter_arma{NN}{1:trial}]));

%% -----------------------------------

Av_Conv_rate_icpa1(NN)=sum(avg_error_icpa1(2:max([final_iter_icpa1{NN}{1:trial}])+1)./avg_error_icpa1(1:max([final_iter_icpa1{NN}{1:trial}])))/max([final_iter_icpa1{NN}{1:trial}]);
Av_Conv_rate_icpa2(NN)=sum(avg_error_icpa2(2:max([final_iter_icpa2{NN}{1:trial}])+1)./avg_error_icpa2(1:max([final_iter_icpa2{NN}{1:trial}])))/max([final_iter_icpa2{NN}{1:trial}]);
Av_Conv_rate_icpa3(NN)=sum(avg_error_icpa3(2:max([final_iter_icpa3{NN}{1:trial}])+1)./avg_error_icpa3(1:max([final_iter_icpa3{NN}{1:trial}])))/max([final_iter_icpa3{NN}{1:trial}]);
Av_Conv_rate_icpa4(NN)=sum(avg_error_icpa4(2:max([final_iter_icpa4{NN}{1:trial}])+1)./avg_error_icpa4(1:max([final_iter_icpa4{NN}{1:trial}])))/max([final_iter_icpa4{NN}{1:trial}]);
Av_Conv_rate_icpa5(NN)=sum(avg_error_icpa5(2:max([final_iter_icpa5{NN}{1:trial}])+1)./avg_error_icpa5(1:max([final_iter_icpa5{NN}{1:trial}])))/max([final_iter_icpa5{NN}{1:trial}]);

Av_Conv_rate_iopa0(NN)=sum(avg_error_iopa0(2:max([final_iter_iopa0{NN}{1:trial}])+1)./avg_error_iopa0(1:max([final_iter_iopa0{NN}{1:trial}])))/max([final_iter_iopa0{NN}{1:trial}]);
Av_Conv_rate_iopa1(NN)=sum(avg_error_iopa1(2:max([final_iter_iopa1{NN}{1:trial}])+1)./avg_error_iopa1(1:max([final_iter_iopa1{NN}{1:trial}])))/max([final_iter_iopa1{NN}{1:trial}]);
Av_Conv_rate_iopa2(NN)=sum(avg_error_iopa2(2:max([final_iter_iopa2{NN}{1:trial}])+1)./avg_error_iopa2(1:max([final_iter_iopa2{NN}{1:trial}])))/max([final_iter_iopa2{NN}{1:trial}]);
Av_Conv_rate_iopa3(NN)=sum(avg_error_iopa3(2:max([final_iter_iopa3{NN}{1:trial}])+1)./avg_error_iopa3(1:max([final_iter_iopa3{NN}{1:trial}])))/max([final_iter_iopa3{NN}{1:trial}]);
Av_Conv_rate_iopa4(NN)=sum(avg_error_iopa4(2:max([final_iter_iopa4{NN}{1:trial}])+1)./avg_error_iopa4(1:max([final_iter_iopa4{NN}{1:trial}])))/max([final_iter_iopa4{NN}{1:trial}]);
Av_Conv_rate_iopa5(NN)=sum(avg_error_iopa5(2:max([final_iter_iopa5{NN}{1:trial}])+1)./avg_error_iopa5(1:max([final_iter_iopa5{NN}{1:trial}])))/max([final_iter_iopa5{NN}{1:trial}]);

Av_Conv_rate_gradient(NN)=sum(avg_error_gradient(2:max([final_iter_gradient{NN}{1:trial}])+1)./avg_error_gradient(1:max([final_iter_gradient{NN}{1:trial}])))/max([final_iter_gradient{NN}{1:trial}]);

Av_Conv_rate_arma(NN)=sum(avg_error_arma(2:max([final_iter_arma{NN}{1:trial}])+1)./avg_error_arma(1:max([final_iter_arma{NN}{1:trial}])))/max([final_iter_arma{NN}{1:trial}]);



%% --------------------------------------
Time_icpa0(NN)=sum([time_icpa0{1:trial}])/(p*trial);
Time_icpa1(NN)=sum([time_icpa1{1:trial}])/(p*trial);
Time_icpa2(NN)=sum([time_icpa2{1:trial}])/(p*trial);
Time_icpa3(NN)=sum([time_icpa3{1:trial}])/(p*trial);
Time_icpa4(NN)=sum([time_icpa4{1:trial}])/(p*trial);
Time_icpa5(NN)=sum([time_icpa5{1:trial}])/(p*trial);

Time_iopa0(NN)=sum([time_iopa0{1:trial}])/(p*trial);
Time_iopa1(NN)=sum([time_iopa1{1:trial}])/(p*trial);
Time_iopa2(NN)=sum([time_iopa2{1:trial}])/(p*trial);
Time_iopa3(NN)=sum([time_iopa3{1:trial}])/(p*trial);
Time_iopa4(NN)=sum([time_iopa4{1:trial}])/(p*trial);
Time_iopa5(NN)=sum([time_iopa5{1:trial}])/(p*trial);

Time_gradient(NN)=sum([time_gradient{1:trial}])/(p*trial);

Time_arma(NN)=sum([time_arma{1:trial}])/(p*trial);

Time_cent(NN)=sum([time_cent{NN}{1:trial}])/(trial);


%% ---------------------------------------
Av_iter_icpa0(NN)=sum([final_iter_icpa0{NN}{1:trial}]+1)/(trial);
Av_iter_icpa1(NN)=sum([final_iter_icpa1{NN}{1:trial}]+1)/(trial);
Av_iter_icpa2(NN)=sum([final_iter_icpa2{NN}{1:trial}]+1)/(trial);
Av_iter_icpa3(NN)=sum([final_iter_icpa3{NN}{1:trial}]+1)/(trial);
Av_iter_icpa4(NN)=sum([final_iter_icpa4{NN}{1:trial}]+1)/(trial);
Av_iter_icpa5(NN)=sum([final_iter_icpa5{NN}{1:trial}]+1)/(trial);

Av_iter_iopa0(NN)=sum([final_iter_iopa0{NN}{1:trial}]+1)/(trial);
Av_iter_iopa1(NN)=sum([final_iter_iopa1{NN}{1:trial}]+1)/(trial);
Av_iter_iopa2(NN)=sum([final_iter_iopa2{NN}{1:trial}]+1)/(trial);
Av_iter_iopa3(NN)=sum([final_iter_iopa3{NN}{1:trial}]+1)/(trial);
Av_iter_iopa4(NN)=sum([final_iter_iopa4{NN}{1:trial}]+1)/(trial);
Av_iter_iopa5(NN)=sum([final_iter_iopa5{NN}{1:trial}]+1)/(trial);

Av_iter_gradient(NN)=sum([final_iter_gradient{NN}{1:trial}]+1)/(trial);

Av_iter_arma(NN)=sum([final_iter_arma{NN}{1:trial}]+1)/(trial);

end

% figure(1)
% plot(1:p, log10(avg_error_icpa0),'-*','LineWidth', 2)
% hold on
% plot(1:p, log10(avg_error_arma),'-x','LineWidth', 2)
% hold on;
% plot(1:p, log10(avg_error_icpa1),'-*','LineWidth', 2)
% hold on
% plot(1:p, log10(avg_error_gradient),'-o','LineWidth', 2)
% hold on
% plot(1:p, log10(avg_error_iopa0),'-x', 'LineWidth', 2)
% hold on
% plot(1:p, log10(avg_error_icpa2),'-*', 'LineWidth', 2)
% hold on
% plot(1:p, log10(avg_error_iopa1),'-x', 'LineWidth', 2)
% hold on
% plot(1:p, log10(avg_error_icpa3),'-*', 'LineWidth', 2)
% hold on
% plot(1:p, log10(avg_error_icpa4),'-*', 'LineWidth', 2)
% hold on
% plot(1:p, log10(avg_error_iopa2),'-x', 'LineWidth', 2)
% hold on
% plot(1:p, log10(avg_error_icpa5),'-*','LineWidth', 2)
% hold on
% plot(1:p, log10(avg_error_iopa3),'-x', 'LineWidth', 2)
% hold on
% plot(1:p, log10(avg_error_iopa4),'-x','LineWidth', 2)
% hold on 
% plot(1:p, log10(avg_error_iopa5),'-x', 'LineWidth', 2)
% hold on
% 
% 
% legend('ICPA0','ARMA','ICPA1','GRAD','IOPA0', 'ICPA2','IOPA1','ICPA3','ICPA4',...
% 'IOPA2', 'ICPA5','IOPA3','IOPA4','IOPA5','Location','best')
% ax = gca;
% ax.FontSize = 15;
% hold off

%method_orders={'ARMA','ICPA1','GRAD','IOPA0', 'ICPA2','IOPA1','ICPA3','ICPA4','IOPA2', 'ICPA5','IOPA3','IOPA4','IOPA5'};

%table=[avg_error_icpa0',avg_error_arma',avg_error_icpa1',avg_error_gradient',avg_error_iopa0',avg_error_icpa2',avg_error_iopa1',avg_error_icpa3',avg_error_icpa4',avg_error_iopa2',avg_error_icpa5',avg_error_iopa3',avg_error_iopa4',avg_error_iopa5'];

% Error table for each N

table;

format short

%Average conv rate (Simulation)
TABLE_Conv=[Av_Conv_rate_arma',Av_Conv_rate_icpa1',Av_Conv_rate_gradient',Av_Conv_rate_iopa0',Av_Conv_rate_icpa2',Av_Conv_rate_iopa1',Av_Conv_rate_icpa3',Av_Conv_rate_icpa4',Av_Conv_rate_iopa2',Av_Conv_rate_icpa5',Av_Conv_rate_iopa3',Av_Conv_rate_iopa4',Av_Conv_rate_iopa5']

%Cconv rate (Theoretical)
TABLE_Theor_conv_rate=[arma_conv_rate_theor',B_K(:,2),A_L(:,1),A_L(:,1),B_K(:,3),A_L(:,2),B_K(:,4),B_K(:,5),A_L(:,3),B_K(:,6),A_L(:,4),A_L(:,5),A_L(:,6)]

% Average number of iteration needed to reach 10^(-3)
Table_Iter=[Av_iter_arma',Av_iter_icpa1',Av_iter_gradient',Av_iter_iopa0',Av_iter_icpa2',Av_iter_iopa1',Av_iter_icpa3',Av_iter_icpa4',Av_iter_iopa2',Av_iter_icpa5',Av_iter_iopa3',Av_iter_iopa4',Av_iter_iopa5']

format shortE

% Average iteration time for each methods
TABLE_Time=[Time_arma',Time_icpa1',Time_gradient',Time_iopa0',Time_icpa2',Time_iopa1',Time_icpa3',Time_icpa4',Time_iopa2',Time_icpa5',Time_iopa3',Time_iopa4',Time_iopa5']

%name of the methods appear on table "TABLE_Total_time_new" (in order)
MethodsOrders={'1:ARMA','2:ICPA1','3:GD0','4:IOPA0','5:ICPA2','6:IOPA1','7:ICPA3','8:ICPA4','9:IOPA2','10:ICPA5','11:IOPA3','12:IOPA4','13:IOPA5'};

% Total time needed to reach 10^(-3)
TABLE_Total_time_new=Table_Iter.*TABLE_Time;

% Total time plot
%Please save this figure 
figure(1)
plot(log10(NumberNode), log10(Time_cent), 'LineWidth', 2)
hold on
plot(log10(NumberNode), log10(TABLE_Total_time_new(:,2)),'-x', 'LineWidth', 2)
hold on
plot(log10(NumberNode), log10(TABLE_Total_time_new(:,5)),'-o', 'LineWidth', 2)
hold on
plot(log10(NumberNode), log10(TABLE_Total_time_new(:,6)),'-.x','LineWidth', 2)
hold on;
plot(log10(NumberNode), log10(TABLE_Total_time_new(:,7)),'--o', 'LineWidth', 2)
hold on
plot(log10(NumberNode), log10(TABLE_Total_time_new(:,8)),'-.x', 'LineWidth', 2)
hold on
plot(log10(NumberNode), log10(TABLE_Total_time_new(:,9)),'-o', 'LineWidth', 2)
hold on
plot(log10(NumberNode), log10(TABLE_Total_time_new(:,10)),'-.x', 'LineWidth', 2)
hold on
plot(log10(NumberNode), log10(TABLE_Total_time_new(:,11)),'-o','LineWidth', 2)
hold on
plot(log10(NumberNode), log10(TABLE_Total_time_new(:,12)),'--x', 'LineWidth', 2)
hold on
plot(log10(NumberNode), log10(TABLE_Total_time_new(:,13)),'--o','LineWidth', 2)
hold on
plot(log10(NumberNode), log10(TABLE_Total_time_new(:,1)),'-.x','LineWidth', 2)
hold on
plot(log10(NumberNode), log10(TABLE_Total_time_new(:,3)),'-o', 'LineWidth', 2)
hold on

xlabel('$$log(N)$$','Interpreter','Latex')
ylabel('$$log(T)$$','Interpreter','Latex')

legend('Centralized','ICPA1','ICPA2', 'IOPA1','ICPA3','ICPA4','IOPA2', 'ICPA5','IOPA3','IOPA4','IOPA5','ARMA','GD0/IOPA0','Location','northwest')
ax = gca;
ax.FontSize = 15;
hold off

% Table for total time needed to reach 10^(-3) (with Centralized method at firts column)
TableOrders={'1:Centralized','2:ARMA','3:GD0','4:ICPA1','5:ICPA2','6:IOPA1','7:ICPA3','8:ICPA4','9:IOPA2','10:ICPA5','11:IOPA3','12:IOPA4','13:IOPA5'};
% Please save this data
TABLE_TOT_TIME=[Time_cent',TABLE_Total_time_new(:,1),TABLE_Total_time_new(:,3),TABLE_Total_time_new(:,2),TABLE_Total_time_new(:,5:13)]; 

%TABLE_Conv=[Av_Conv_rate_arma',Av_Conv_rate_icpa1',Av_Conv_rate_gradient',Av_Conv_rate_iopa0',Av_Conv_rate_icpa2',Av_Conv_rate_iopa1',Av_Conv_rate_icpa3',Av_Conv_rate_icpa4',Av_Conv_rate_iopa2',Av_Conv_rate_icpa5',Av_Conv_rate_iopa3',Av_Conv_rate_iopa4',Av_Conv_rate_iopa5'];

%TABLE=[TABLE_Time(:,1),TABLE_Time(:,3),TABLE_Time(:,2),TABLE_Time(:,5:13)]; 
 
% figure(2)
% plot(NumberNode, TABLE_Conv(:,1),'-o','LineWidth', 2)
% hold on
% plot(NumberNode, TABLE_Conv(:,2),'-o', 'LineWidth', 2)
% hold on
% plot(NumberNode, TABLE_Conv(:,3),'-o', 'LineWidth', 2)
% hold on
% plot(NumberNode, TABLE_Conv(:,5),'-o', 'LineWidth', 2)
% hold on
% plot(NumberNode, TABLE_Conv(:,6),'-o','LineWidth', 2)
% hold on;
% plot(NumberNode, TABLE_Conv(:,7),'-x', 'LineWidth', 2)
% hold on
% plot(NumberNode, TABLE_Conv(:,8),'-x', 'LineWidth', 2)
% hold on
% plot(NumberNode, TABLE_Conv(:,9),'-x', 'LineWidth', 2)
% hold on
% plot(NumberNode, TABLE_Conv(:,10),'-x', 'LineWidth', 2)
% hold on
% plot(NumberNode, TABLE_Conv(:,11),'-x','LineWidth', 2)
% hold on
% plot(NumberNode, TABLE_Conv(:,12),'-x', 'LineWidth', 2)
% hold on
% plot(NumberNode, TABLE_Conv(:,13),'-x','LineWidth', 2)
% hold on 
% xticks(NumberNode)
% 
% 
% legend('ARMA','ICPA1','GD0/IOPA0', 'ICPA2','IOPA1','ICPA3','ICPA4','IOPA2', 'ICPA5','IOPA3','IOPA4','IOPA5','Location','northeast')
% ax = gca;
% ax.FontSize = 15;
% hold off   


%%%save('Revision_CirculantGraph_Different_N_trial_1000_until_N2000')
