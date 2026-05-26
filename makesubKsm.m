function subK = makesubKsm(fem,n,T,dofs)
[K0,~,] = fem.matrix ;
Nd = length(dofs);
all_dofs = 1:size(K0,1);
slave_dofs = setdiff(all_dofs, dofs);
Nn = length(slave_dofs);

% 分块

for i = 1:n
    theta = zeros(n,1);theta(i) = 1;
    fem0 = FEMUpt_4CSM_v5(fem,theta);
    [K1,~,] = fem0.matrix;
    subK0 = K1-K0;
    K_mm = subK0(dofs, dofs);
    K_ms = subK0(dofs, slave_dofs);
    K_sm = K_ms';
    K_ss = subK0(slave_dofs, slave_dofs);
    K_full = [K_mm K_ms; K_sm K_ss];
    subK(:,:,i) = T' * K_full * T;
end
end
