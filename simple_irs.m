function [K0, M0, T_irs] = simple_irs(K, M, dofs)
% 最小修正的IRS

Nd = length(dofs);
all_dofs = 1:size(K,1);
slave_dofs = setdiff(all_dofs, dofs);
Nn = length(slave_dofs);

% 分块
K_mm = K(dofs, dofs);
K_ms = K(dofs, slave_dofs);
K_sm = K_ms';
K_ss = K(slave_dofs, slave_dofs);

M_mm = M(dofs, dofs);
M_ms = M(dofs, slave_dofs);
M_sm = M_ms';
M_ss = M(slave_dofs, slave_dofs);

% 静态缩聚
Phi = -K_ss \ K_sm;  % 关键：这是从自由度的变换 [Nn × Nd]
T_static = [eye(Nd); Phi];

% 完整矩阵
K_full = [K_mm K_ms; K_sm K_ss];
M_full = [M_mm M_ms; M_sm M_ss];

% 静态缩聚矩阵
M_static = T_static' * M_full * T_static;
K_static = T_static' * K_full * T_static;


% 计算 M_combined = M_sm + M_ss * Phi
M_combined = M_sm + M_ss * Phi;  % [Nn × Nd]

% 计算 Theta
Theta = M_static \ K_static;

% 计算修正项
S = -K_ss \ (M_combined * Theta);  % [Nn × Nd]

% IRS变换
T_irs = T_static + [zeros(Nd, Nd); S];

% 缩聚矩阵
K0 = T_irs' * K_full * T_irs;
M0 = T_irs' * M_full * T_irs;

% 2. 振型计算和验证
% 完整模型
[V_full, D_full] = eigs(K, M, max(20), 'smallestabs' );
[~, idx_full] = sort(diag(D_full));
freq_full = sqrt(diag(D_full(idx_full, idx_full))) / (2*pi);

% 缩聚模型
[V_red, D_red] = eig(K0, M0);
[~, idx_red] = sort(diag(D_red));
freq_red = sqrt(diag(D_red(idx_red, idx_red))) / (2*pi);

% 频率误差
N_compare = min(15, length(freq_red));
freq_error = abs(freq_full(1:N_compare) - freq_red(1:N_compare)) ./ freq_full(1:N_compare);

% 振型重建
reorder_idx = [dofs, slave_dofs];
Phi_reconstructed = zeros(size(K,1), N_compare);

for i = 1:N_compare
    phi_reordered = T_irs * V_red(:, idx_red(i));
    Phi_reconstructed(reorder_idx, i) = phi_reordered;
end

% MAC计算
MAC = zeros(N_compare);
for i = 1:N_compare
    for j = 1:N_compare
        phi_i = V_full(dofs, idx_full(i));
        phi_j = Phi_reconstructed(dofs, j);
        MAC(i,j) = abs(phi_i' * phi_j)^2 / ((phi_i' * phi_i) * (phi_j' * phi_j));
    end
end
MAC_diag = zeros(N_compare, 1);

for i = 1:N_compare
    % 原始模型第i阶模态的主自由度部分
    phi_original = V_full(dofs, idx_full(i));
    phi_reduced = V_red(:, idx_red(i));

   scale = norm(phi_original) / norm(phi_reduced);
    phi_reduced_scaled = phi_reduced * scale;

    % 计算MAC值（只比较主自由度部分）
    mac_value = abs(phi_original' * phi_reduced_scaled)^2 / ...
        ((phi_original' * phi_original) * (phi_reduced_scaled' * phi_reduced_scaled));

    MAC_diag(i) = mac_value;
end


end