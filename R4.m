% function rhat = R4(a4, kl)
% m = 4;
% n = 100;
% rhat = [];
% 
% for k = 601:50:kl+500
%     x = a4(1:m, k:k+99);
% 
%     % 链内均值
%     xbar_j = mean(x, 2);
%     xbar_all = mean(xbar_j);
% 
%     % B/n: 链间方差
%     B_over_n = sum((xbar_j - xbar_all).^2) / (m - 1);
% 
%     % W: 链内方差（无偏）
%     s_j2 = var(x, 0, 2); % 每个链的方差（除以 n-1）
%     W = mean(s_j2);
% 
%     % σ²₊
%     sigma2_plus = (n - 1)/n * W + B_over_n;
% 
%     % V_hat
%     V_hat = sigma2_plus + B_over_n / m;
% 
%     % R_hat（未校正）
%     R_hat_uncorrected = V_hat / W;
% 
%     % 取平方根
%     rhat_val = sqrt(R_hat_uncorrected);
% 
%     rhat = [rhat rhat_val];
% end
% end
function[rhat] = R4(a4, kl)
%%
m = 4; n = 100; rhat = [];

for k = 601:50:kl+500
    x = a4(1:m, k:k+99);

    % 计算链内均值
    xba = mean(x, 2);

    % 计算链间方差 B/n
    B_dividen = var(xba, 0); 
    B = B_dividen * n;

    % 计算链内方差 W
    W = mean(var(x, 0, 2)); % 每个链的方差，然后取平均

    % 计算总体方差估计 (公式1.1)
    sigma2_plus = ((n-1)/n)*W + B_dividen;

    % 计算 pooled posterior variance estimate (包含抽样变异性)
    V_hat = sigma2_plus + B_dividen/m;

    % 计算自由度估计 (基于矩方法)
    % 首先需要估计var(V_hat)
    % 计算s的方差
    s = var(x, 0, 2); % 每个链的方差
    var_s = var(s);

    % 计算协方差项
    covariance1 = cov(s, xba.^2);
    cov_s_xba2 = covariance1(1, 2);

    covariance2 = cov(s, xba);
    cov_s_xba = covariance2(1, 2);

    xmean = mean(xba);

    % 估计var(V_hat) (公式推导见论文)
    var_Vhat = ((n-1)/n)^2/m * var_s + ...
               ((m+1)/(m*n))^2 * 2/(m-1) * B^2 + ...
               2*(m+1)*(n-1)/(m*n^2) * n/m * (cov_s_xba2 - 2*xmean*cov_s_xba);

    % 估计自由度
    df = 2 * V_hat^2 / var_Vhat;

    % 计算校正后的Rhat (公式1.2和1.3)
    R_hat_uncorrected = V_hat / W;
    R_hat_corrected = ((df + 3) / (df + 1)) * R_hat_uncorrected;

    % 取平方根得到最终的scale reduction factor
    Rhat = sqrt(R_hat_corrected);

    rhat = [rhat Rhat];
end
end