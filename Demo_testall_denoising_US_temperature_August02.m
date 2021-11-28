
% This simulation compares the denoising result of ICPA1, IOPA1 and GD for US_Hourly_2010_August_1st
% It requires external MVInvAprxChebCoeff2 and PolFilter_H, where
% MVInvAprxChebCoeff2 gives the coefficients of chebyshev pol.approx and
% PolFilter_H gives polynomial filter H

% For different noise level change the value of eta
% Before running the simulation change trial value to 1000
% save the data and figure at the end

% Fig1 and Fig2 are 1st and 12th snapshot of a signal

clear;clc;format short

%% ---------------- Data -----------------
load US_Hourly_2010_August_1st
%load JointEigenvalues_for_TVGD

X=signal;

Norm_X=norm(X,'fro');

[N,M]=size(X); 

%% -------------------------- Graph signal plot -----------------------

lim1=min(min(X(:,1)),min(X(:,M/2)));
lim2=max(max(X(:,1)),max(X(:,M/2)));

figure(1) 
pp=graph(A);
Plot=plot(pp,'NodeCdata',X(:,1));
Plot.XData=x;
Plot.YData=y;
Plot.MarkerSize=7;
Plot.LineWidth=1;
Plot.NodeFontSize=3;
xlabel('Temperature data $${\bf w}_1$$','Interpreter','Latex')
ax = gca;
ax.FontSize = 20;
colormap jet
colorbar
caxis([lim1,lim2])
xlim([min(x),max(x)])
ylim([min(y),max(y)])
set(gca,'Visible','off')

figure(2)
pp=graph(A);
Plot=plot(pp,'NodeCdata',X(:,M/12));
Plot.XData=x;
Plot.YData=y;
Plot.MarkerSize=7;
Plot.LineWidth=1;
Plot.NodeFontSize=3;
xlabel('Temperature data $${\bf w}_{12}$$','Interpreter','Latex')
ax = gca;
ax.FontSize = 20;
colormap jet
colorbar
caxis([lim1,lim2])
xlim([min(x),max(x)])
ylim([min(y),max(y)])
set(gca,'Visible','off')

%% -----------------------Parameters--------------------------

K=1; %maximum approximation power ICPA

L=1; %maximum approximation power for IOPA

p=6; %number of iteration

eta=5; %noise level 

trial=1000; % number of trials


%% ----------------- Shifts S1 and S2--------------------------------------

Adjacency_matrix=A; % Adjacency matrix of a graph W

[N,M]=size(X);

T=zeros(M,M);
T(1:M-1,2:M)=eye(M-1,M-1);
T(1,M)=1;
C_g=(T+T');
D=diag(sum(C_g));
L_C=D-C_g;
L_C=diag(diag(D).^(-1/2))*L_C*diag(diag(D).^(-1/2)); % normalized laplacian matrix of a cycle graph C


D_W=diag(sum(A));
L_W=D_W-A;
L_W=diag(diag(D_W).^(-1/2))*L_W*diag(diag(D_W).^(-1/2)); % normalized laplacian matrix of a graph W

S1=kron_sparse(eye(M,M),L_W); %Shift1 IxL_W
S2=kron_sparse(L_C,eye(N,N)); %Shift2 L_CxI

%% -------------------------alpha, beta ----------------------------------


% alpha=M*N*(eta^2/3)/(X(:)'*S1*X(:)+M*N*(eta^2/3));
% beta=M*N*(eta^2/3)/(X(:)'*S2*X(:)+2*(M-1)*N*(eta^2/3));

 alpha=M*N*(eta^2/3)/(X(:)'*S1*X(:)+M*N*(eta^2/3));
 beta=M*N*(eta^2/3)/(X(:)'*S2*X(:)+M*N*(eta^2/3));


AB=[alpha,0;                %(alpha,0)
    0,beta;                 %(0,beta)
    alpha,beta];            %(alpha,beta)

%% -------Start of the simulation for different cases (alpha,0),(0,beta) and (alpha,beta)

avg_snr_icpa={};
avg_snr_iopa={};
avg_snr_gradient={};
avg_snr_centralized={};
avg_snr_input={};

avg_X_icpa={};
avg_X_iopa={};
avg_X_gradient={};
avg_X_centralized={};
avg_noise={};

for qq=1:3 % for (alpha,0), (0,beta), and (alpha,beta)

alpha=[];beta=[];

alpha=AB(qq,1); beta=AB(qq,2);

dsp2=['Testing for alpha=', num2str(alpha), ', beta=', num2str(beta)];
disp(dsp2)

%% ----------------------------Filter H ----------------
h=[];

% h(x,y)=\sum_i [h(i,3)*(x^(h(i,1))*y^(h(i,2))]
% first column is a power of x, second column is a power of y, and 3rd
% column is a corresponding coefficient

h=[0,0,1;
   1,0,alpha;
   0,1,beta]; % polynomial coefficient


[r_h,c_h]=size(h);

H=eye(N*M)+alpha*S1+beta*S2; %Filter H

%% ----------------Filter G_K--------------------------------------------------

C=[];
C=MVInvAprxChebCoeff2(h,K,0,2,0,2); % Chebyshev coefficients until power K for two shifts

G_K=[];
G_K=C(1,3)*eye(N*M,N*M)+C(2,3)*(S2-eye(N*M))+C(3,3)*(S1-eye(N*M)); % G_K is a Kth order Chebyshev polynomial approximation of an inverse filter (K=1).

%G_K=ChebFilter_G(C,[S1,S2],K,[0,2,0,2]); % inverse filter G_K, where G_K{k} is k-1 power approximation 


%% -------------Filter G_L---------------------------------------------------

[Q1,QQ1]=schur(full(L_W));
[Q2,QQ2]=schur(full(L_C));

QQ=kron_sparse(Q2,Q1); % Common unitary matrix

E1=diag(QQ'*S1*QQ);  % eigenvalue of a shift 1
E2=diag(QQ'*S2*QQ); % eigenvalue of a shift 2 with corresponding order

EE=[E1,E2]; % Joint eigenvalues 

H_eig=zeros(N*M,1);

for j=1:r_h
    H_eig=H_eig+h(j,c_h)*prod((EE.^h(j,1:c_h-1))')'; % simple way to calculate h(lambda_i)
end

lmn=min(abs(H_eig));  % abs min eigenvalue of H
lmx=max(abs(H_eig));  % abs max eigenvalue of H


ch=[];

    
t=1;
for i=1:L+1
    for j=1:L+1
        if i+j-2<=L
           ch(t,:)=[j-1,i-1]; % defines powers of shifts such that i+j<L for x^i*y^j
           t=t+1;
        end
    end
end


[r_ch,c_ch]=size(ch);

M1=zeros(N*M,r_ch);

        
for j=1:r_ch
    M1(:,j)=prod((EE.^ch(j,:))')';
end

A1=diag(H_eig);

AA=A1*M1; %------diag(h(\lambda))*(eigenvalue power ordering)

[sa,sb]=size(AA);

% setup for linear programming

A_eq=[-ones(N*M,1) AA;-ones(N*M,1) -AA];
b_eq=[ones(1,N*M) -ones(1,N*M)] ;
coef_eq=[1 zeros(1,sb)];

% Linear programming to solve optimal pol. approx.
s_op=[];
s_op=linprog(coef_eq,A_eq,b_eq); % s_op contains the largest eigenvalue of G and 
                               % the coefficient for GL 
                               

ch(:,c_ch+1)=s_op(2:length(s_op)); % coefficient of polynomial approx.   

G_L=[];

G_L=PolFilter_H(ch,[S1,S2]); %optimal polynomial approx. filter G_L


gamma=2*(lmn+ lmx)^(-1); % optimal parameter for Gradient descent


sum_snr_icpa=zeros(1,p);
sum_snr_iopa=zeros(1,p);
sum_snr_gradient=zeros(1,p);
sum_snr_centralized=0;
sum_snr_input=0;

sum_X_icpa=zeros(N,p*M);
sum_X_iopa=zeros(N,p*M);
sum_X_gradient=zeros(N,p*M);
sum_X_centralized=zeros(N,M);
sum_noise=zeros(N,M);

%% ---------- Start the simulation for 1000 trial -----------------

for ii=1:trial

dsp3=['trial=', num2str(ii)];
disp(dsp3)

%% -------------Noise----------------------------
noise=unifrnd(-eta,eta,N,M);
BN=X+noise; %noisy data
b_noise=BN(:);

%% ---------------------ICPA-----------------------------

x_icpa=sparse(zeros(N*M,1));
z=sparse(zeros(N*M,1));
bb=b_noise;

X_icpa=[];
SNR_icpa=[];

for m=1:p
    z=sparse(G_K*bb);
    x_icpa=sparse(x_icpa+z);
    bb=sparse(bb-H*z);
    
    X_icpa=[X_icpa reshape(x_icpa,N,M)];

    SNR_icpa=[SNR_icpa -20*log10(norm(X-reshape(x_icpa,N,M),'fro')/Norm_X)];
end

%% --------------------IOPA--------------------------------

x_iopa=zeros(N*M,1);
zz=zeros(N*M,1);

bbb=b_noise;
X_iopa=[];
SNR_iopa=[];
for m=1:p
    zz=sparse(G_L*bbb);
    x_iopa=sparse(x_iopa+zz);
    bbb=sparse(bbb-H*zz);
    
    X_iopa=[X_iopa reshape(x_iopa,N,M)];
    
    SNR_iopa=[SNR_iopa -20*log10(norm(X-reshape(x_iopa,N,M),'fro')/Norm_X)];
    
end

%% -----------------Gradient Descent--------------

SNR_GD=[];
x_gradient=zeros(N*M,1);
X_gradient=[];
SNR_gradient=[];

for m=1:p
    x_gradient=x_gradient-gamma*(H*x_gradient-b_noise);
    
    X_gradient=[X_gradient reshape(x_gradient,N,M)];

    SNR_gradient=[SNR_gradient -20*log10(norm(X-reshape(x_gradient,N,M),'fro')/Norm_X)];
end
%% ------------ Centralized  Methdod --------------

X_centralized = reshape(H\b_noise, N,M);

SNR_centralized=-20*log10(norm(X-X_centralized,'fro')/Norm_X);

SNR_input=-20*log10(norm(X-BN,'fro')/Norm_X);


sum_snr_icpa=sum_snr_icpa+SNR_icpa;
sum_snr_iopa=sum_snr_iopa+SNR_iopa;
sum_snr_gradient=sum_snr_gradient+SNR_gradient;
sum_snr_centralized=sum_snr_centralized+SNR_centralized;
sum_snr_input=sum_snr_input+SNR_input;


sum_X_icpa=sum_X_icpa+X_icpa;
sum_X_iopa=sum_X_iopa+X_iopa;
sum_X_gradient=sum_X_gradient+X_gradient;

sum_X_centralized=sum_X_centralized+X_centralized;

sum_noise=sum_noise+noise;

end


avg_snr_icpa{qq}=sum_snr_icpa/trial;
avg_snr_iopa{qq}=sum_snr_iopa/trial;
avg_snr_gradient{qq}=sum_snr_gradient/trial;
avg_snr_centralized{qq}=sum_snr_centralized/trial;
avg_snr_input{qq}=sum_snr_input/trial;

avg_X_icpa{qq}=sum_X_icpa/trial;
avg_X_iopa{qq}=sum_X_iopa/trial;
avg_X_gradient{qq}=sum_X_gradient/trial;
avg_X_centralized{qq}=sum_X_centralized/trial;
avg_noise{qq}=sum_noise/trial;

param_alpha_beta{qq}=AB;

end


%% ----------------------Table -----------------------------------

Methods_order={'(alpha,0)';'(0,beta)';'(alpha,beta)'};

table_input=[avg_snr_input{1};avg_snr_input{2};avg_snr_input{3}];

% -----------each of this tables are 3x(p+1) where 3 is for
% (alpha,0),(0,beta) and (alpha,beta) and 1:p is for iteration and p+1 is
% a result for infinity (centralized)

table_iopa=[[avg_snr_iopa{1}; avg_snr_iopa{2}; avg_snr_iopa{3}],[avg_snr_centralized{1}; avg_snr_centralized{2}; avg_snr_centralized{3}]];

table_icpa=[[avg_snr_icpa{1}; avg_snr_icpa{2}; avg_snr_icpa{3}],[avg_snr_centralized{1}; avg_snr_centralized{2}; avg_snr_centralized{3}]];
    
table_gradient=[[avg_snr_gradient{1}; avg_snr_gradient{2}; avg_snr_gradient{3}],[avg_snr_centralized{1}; avg_snr_centralized{2}; avg_snr_centralized{3}]];

%table is 3x1
    
table_centralized=[avg_snr_centralized{1}; avg_snr_centralized{2}; avg_snr_centralized{3}];

%save('xxxxxRevision_Denoised_US_hourly_Data_eta_5','Methods_order','table_input','table_iopa','table_icpa','table_gradient','avg_noise','avg_X_icpa','avg_X_iopa','avg_X_gradient','avg_X_centralized','param_alpha_beta','X','A')

