% This is the main program for the eight-story shear model. 
% It presents the complete workflow under the condition of no missing data. 
% Users can modify and adapt it to other scenarios as needed.

clear;close all;clc

%% FEM
Nk = 24; % Unknown parameters
fem = FEM_4CSM_v2;
dofs = [25, 43, 49, 67, 73, 91, 97, 115, 121, 139, 145, 163, 169, 187, ...
    26, 44, 50, 68, 74, 92, 98, 116, 122, 140, 146, 164, 170, 188, ...
    193, 211, 194, 212];% Channel index

[K0,M0] = fem.matrix ;
[~,phifem0] = fem.L_Slv_Mod(1:20);
phifem = phifem0(dofs,:);
Nd = length(dofs);
[K0, M, Tirs] = simple_irs(K0, M0, dofs);
subK = makesubKsm(fem,Nk,Tirs,dofs);

%% Data
load shearmodel.mat
phishear([1:14 29 30],:,:) = -phishear([1:14 29 30],:,:);

Lm = 1:32; % measured dofs
bs = 12;% measured modes
fshear = fshear(1:bs,:);
lambda = (fshear*2*pi).^2;
[~,T] = size(lambda);
for i = 1:T
    phi_m(:,:,i) = mynormc(phishear(Lm,1:bs,i));
end
Lu = setdiff(1:Nd,Lm);
Sa = eye(Nd);
S = Sa(Lm,:);
Su = Sa(Lu,:);
Nm = size(fshear,1);

%% MwG sampling
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

flag = true; % Step size needs to be adjusted
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
    subplot(4,Nk/4,i);
    f1 = normpdf(kap3(i)-0.3:0.005:kap3(i)+0.3,kap3(i),std0(i));
    histogram(kappat(i,:),30,'Normalization','pdf'); hold on;
    %plot(kap3(i)-0.3:0.005:kap3(i)+0.3,f1,'r','Linewidth',2);
    xlim([kap3(i)-3*std0(i) kap3(i)+3*std0(i)]); ylim([0 1.1*max(f1)]);
    hold off; xlabel(['\kappa_',num2str(i)]); ylabel('PDF');
end

figure;
j = 1;
max_pdf_vals = zeros(1, ceil(Nk/3));
for i = 2:3:Nk
    [counts, edges] = histcounts(kappat(i,:), 30, 'Normalization', 'pdf');
    max_pdf_vals(j) = max(counts);
    j = j + 1;
end
max_y = max(max_pdf_vals) * 1.1;
j = 1;
for i = 2:3:Nk
    subplot(2,4,j);
    h = histogram(kappat(i,:), 30, 'Normalization', 'pdf', 'FaceAlpha', 1);
    current_max = max(h.Values);
    hold on;
    ylim([0 current_max * 1.1]); 
    line([mean(kap3(i,:)) mean(kap3(i,:))], ylim, ...
        'Color', 'r', 'LineWidth', 2, 'LineStyle', '-');
    xlim([mean(kap3(i,:))-3*std0(i) mean(kap3(i,:))+3*std0(i)]);
    ylabel('PDF', 'FontSize', 14);
    xlabel(['\it{\kappa_t}_{',num2str(i),'}']);
    hold off;
    set(gca, 'FontSize', 14);
    j = j + 1;
end

%% Evaluation of Update Effectiveness

akappat = zeros(Nk,T);
for i = 1:Nk
    abc = [i Nk+i 2*Nk+i 3*Nk+i];
    akappat(i,:) = mean(mean1(abc,:),1);
end
fem = FEM_4CSM_v2;
dofs = [25 43 49 67 73 91 97 115 121 139 145 163 169 187 26 44 50 68 74 92 98 116 122 140 146 164 170 188 193 211 194 212];
fem0 = FEMUpt_4CSM_v5(fem, akappat(:,1));
[feq,phifem0] = fem0.L_Slv_Mod(1:20);
load shearmodel.mat
phishear([1:14 29 30],:,:) = -phishear([1:14 29 30],:,:);

phiall = mynormc(phifem0(dofs,:));

mode_shapes1 = phishear(:,:,1);mode_shapes2 = phiall;
mac_value = [];
[~, num_modes1] = size(mode_shapes1);
[~, num_modes2] = size(mode_shapes2);

best_match = zeros(num_modes1, 1);
max_MAC = zeros(num_modes1, 1);
frequency_diff = zeros(num_modes1,1);
for i = 1:num_modes1
    max_mac_for_mode = 0;
    best_mode_index = 0;
    for j = 1:num_modes2
        phi1 = mode_shapes1(:, i);
        phi2 = mode_shapes2(:, j);
        mac_value(i,j) = (abs(phi1' * phi2)^2) / (phi1' * phi1) / (phi2' * phi2);
        if mac_value(i,j) > max_mac_for_mode
            max_mac_for_mode = mac_value(i,j);
            best_mode_index = j;
        end
    end
    best_match(i) = best_mode_index;
    max_MAC(i) = max_mac_for_mode;
    f1 = fshear(i,1);
    f2(i) = feq(best_mode_index);
    frequency_diff(i) = (f1 - f2(i)) / f1 * 100;
end

phi = mynormc(phishear(:,:,1));
phi1 = phi(1:32, 9);
phi2 = phi(1:32, 10);

A = [phi1, phi2];
phi3 = mode_shapes2(:,9);phi4 = mode_shapes2(:,10);
B = [phi3, phi4];

angle_rad = subspace(A,B);
cos_theta1 = cos(angle_rad);

phi1 = phi(1:32, 4);
phi2 = phi(1:32, 5);

A = [phi1, phi2];
phi3 = mode_shapes2(:,4);phi4 = mode_shapes2(:,5);
B = [phi3, phi4];

angle_rad = subspace(A,B);
cos_theta2 = cos(angle_rad);
macfinal = [max_MAC(1:3);cos_theta2;max_MAC(6:8);cos_theta1;max_MAC(11:12)];
feqfinal = [frequency_diff(1:3);mean(frequency_diff(4:5));frequency_diff(6:8);mean(frequency_diff(9:10));frequency_diff(11:12)];


%% Parameter evolution over time
figure;
eo = 1:T;
test1 = mean(kap3,2);
colors = lines(Nk); 

for i = 1:Nk
    subplot(4, Nk/4, i);
    std2 = std(squeeze(kappat7(i,:,:)),0,2);
    errorbar(eo, akappat(i,:), std2, ...
        'Color', colors(i,:), 'LineStyle', 'none', ...
        'LineWidth', 0.8, 'CapSize', 5);hold on   
    plot(eo, akappat(i,:), '--o', 'Color', colors(i,:), ...
        'MarkerFaceColor', colors(i,:), 'LineWidth', 1);hold on

    plot([1 T], [test1(i) test1(i)], '-', ...
        'Color', colors(i,:)*0.5, 'LineWidth', 1.5);

    grid on;
    xlabel('\itT'); ylim([test1(i)-std2(end)-0.1 test1(i)+std2(end)+0.1]);xlim([0 23]);
    ylabel(sprintf('\\kappa_{%d} ', i));set(gca, 'FontSize', 14);

    hold off;
end


