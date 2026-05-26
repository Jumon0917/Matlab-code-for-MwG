%% 8-DOF mass-spring model
clear;close all;clc

%% model information
M = eye(8);
%Stiffness matrix
diagnal0=[2400,4000,5600,7200,8800,10400,12000,13600];
diagnal1=[-1600,-2400,-3200,-4000,-4800,-5600,-6400];
K=diag(diagnal0)+diag(diagnal1,1)+diag(diagnal1,-1);

Nd = 8; Nk = 9;
subK = zeros(Nd,Nd,Nk);
% floor stiffness
kk = 800*(1:9);
subK(1,1,1) = kk(1);
subK(Nd,Nd,Nk) = kk(Nk);
for j = 2:Nk-1
    subK(j-1:j,j-1:j,j) = [1 -1;-1 1]*kk(j);
end
K0 = sum(subK,3); % norminal stiffness
kappa = linspace(-0.2,0.2,9)';%
% alpha(12) = -0.1;
dK = bsxfun(@times,subK,reshape(kappa,[1 1 Nk]));  % deviation of stiffness
Kt = K0+sum(dK,3); % true stiffness

% eigenvalue/eigenvector
[V,Dt] = eig(Kt,M);
[W,k] = sort(diag(Dt));
Omega = sqrt(W); % analytical circular frequencies
fa = Omega/2/pi;
V = V(:,k);

f0 = sqrt(eig(K0,M))/2/pi; % nominal model

%% measurement
Lm = [1:8]; % measured dofs
No = length(Lm);
Lu = setdiff(1:Nd,Lm);
Sa = eye(Nd);
S = Sa(Lm,:);
Su = Sa(Lu,:);

% virtual measured modal information
Ns = 500; % # samples
py = 600;
T = Ns;
lambda = zeros(Nd,Ns);
phi_m = zeros(Nd,Nd,Ns);
rng(1e8);
alpha = mvnrnd(repmat(kappa',Ns,1),25e-4*eye(Nk))';
rng(1e9); n1 = randn(Ns,1);
rng(1.5e9); n2 = randn(Ns,1);
for s = 1:Ns
    dKs = bsxfun(@times,subK,reshape(alpha(:,s),[1 1 Nk]));
    Ks = K0+sum(dKs,3); % true stiffness
    [Vs,Ds] = eig(Ks,M);
    [Ws,ks] = sort(diag(Ds));
    Vs = Vs(:,ks);
    lambda(:,s) = Ws*(1+1e-2*n1(s));
    phi_m(:,:,s) = Vs*(1+1e-2*n2(s));
end

Kupt = [];

% Suppose 1-3 & 5 modes are identified
modesid =1:8; Nm = length(modesid);  % # modes
for t = 1:Ns
    phi_m0(:,:,t) = mynormc(S*phi_m(:,modesid,t));
end
phi_m = phi_m0;
lambda = lambda(modesid,:);
phi_m = phi_m(:,:,1:Ns);
lambda = lambda(:,1:Ns);
testl = sqrt(lambda)/2/pi;
row_std = std(testl, 0, 2); 
frow = mean(testl,2);
%% Gibbs sampling
% Known parameters
Theta.par.K0 = K0;
Theta.par.subK = subK;
Theta.par.M = M;
Theta.par.S = S;
Theta.par.Su = Su;

% Measurements
Theta.data.lambda = lambda;
Theta.data.phi_m = phi_m;

% Initial setting
Theta.prior.P.e0 = Nk;
Theta.prior.P.E0 = eye(Nk)/(Nk*1e-2*eye(Nk));
Theta.prior.kappa.kappa0 = 0*ones(Nk,1);
Theta.prior.kappa.p0 = 1;
Theta.prior.alpha.u0 = ones(Nm,1)*1e-5;
Theta.prior.alpha.h0 = ones(Nm,1)*1e-1;
Theta.prior.beta.r0 = ones(Nm,1)*1e-5;
Theta.prior.beta.s0 = ones(Nm,1)*1e-1;
OPT.Nsamp = py; % # samples

flag = false; % Step size needs to be adjusted
[a4,kappat3,kap2] = MwGsample(Theta,flag,Nm,Nk);


kl = size(a4,2)-600;
[rhat] = R4(a4,kl);
rhatcount = 601:50:kl+500;
hi = figure;
plot(rhatcount,rhat);
[rhatx,rhaty] = ginput(1);

rhatx = fix(rhatx);
kapt0 = zeros(Nk*4,T*kl);
for i = 1:4
    kapt0((i-1)*Nk+1:i*Nk,:) = kappat3(:,:,i);
end
kapt = zeros(4*Nk,T,kl);
for j = 1:kl
    kapt(:,:,j) = kapt0(:,T*(j-1)+1:T*j);
end
mean1 = mean(kapt(:,:,rhatx-600+1:end),3);% The mean of kappat after convergence

kappat6 = kapt(:,:,rhatx-600+1:end);
kappat7 = zeros(Nk,T,4*size(kappat6,3));
kappat8 = zeros(Nk,T*4*size(kappat6,3));

for j = 1:4
    kappat7(:,:,size(kappat6,3)*(j-1)+1:size(kappat6,3)*j) = kappat6(Nk*(j-1)+1:Nk*j,:,:);
end
for j = 1:4*size(kappat6,3)
    kappat8(:,T*(j-1)+1:T*j) = kappat7(:,:,j);% All samples of kappat  after convergence
end
kap3 = mean(kap2(:,rhatx-600+1:end,:),3);% The mean of kappa after convergence


%% Distribution of each kappat

kappat = reshape(kappat8,Nk,[]);
for i = 1:Nk
    std0(i) = std(kappat(i,:));
end
figure
for i = 1:Nk
    subplot(3,Nk/3,i);
    f1 = normpdf(kap3(i)-0.3:0.005:kap3(i)+0.3,kap3(i),std0(i));
    histogram(kappat(i,:),20,'Normalization','pdf'); hold on;
    plot(kap3(i)-0.3:0.005:kap3(i)+0.3,f1,'r','Linewidth',2);
    xlim([kap3(i)-3*std0(i) kap3(i)+3*std0(i)]); ylim([0 1.1*max(f1)]);
    hold off; xlabel(['\kappa_',num2str(i)]); ylabel('PDF');
end


[rhat] = R4(a4,3400);
rhatcount = 601:50:3400+500;
figure;
plot(rhatcount, rhat, 'LineWidth', 1.5);
hold on;
yline(1.05, 'r--', 'LineWidth', 1.5, 'Color', 'r');
ylim([0.98 1.1]);xlim([600 3851]);
xlabel('Iteration', 'FontSize', 14);
ylabel('Gelman-Rubin index', 'FontSize', 14);
set(gca, 'FontSize', 14); 
xticks([600 1000 1500 2000 2500 3000 3500 3900]);  
xticklabels({'600','1000','1500','2000','2500','3000','3500','3900'}); 

