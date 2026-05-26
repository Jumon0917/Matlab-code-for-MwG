function [a4,kappat2,kap2] = MwGsample(Theta,flag,Nm,Nk)
%% Gibbs sampling
% Known parameters
K0 = Theta.par.K0;              % norminal stiffness
subK = Theta.par.subK;          % substructure stiffness
M = Theta.par.M;                % mass
S = Theta.par.S;                % measured selection matrix
Su = Theta.par.Su;              % unmeasured selection matrix
lambda = Theta.data.lambda;     % eigenvalues Nm by T
phi_m = Theta.data.phi_m;       % mode shapes No by Nm by T

% Initial setting
Theta.prior.P.e0 = Nk;
Theta.prior.P.E0 = eye(Nk)/(Nk*1e-2*eye(Nk));
Theta.prior.kappa.kappa0 = 0*ones(Nk,1);
Theta.prior.kappa.p0 = 1;
Theta.prior.alpha.u0 = ones(Nm,1)*1e-5;
Theta.prior.alpha.h0 = ones(Nm,1)*1e-1;
Theta.prior.beta.r0 = ones(Nm,1)*1e-5;
Theta.prior.beta.s0 = ones(Nm,1)*1e-1;
OPT.Nsamp = 600; % # samples

% normalization
nm = norm(M);
lambda0 = median(eig(K0/nm));

Mp = M/nm;
K0p = K0/(nm*lambda0);
subKp = subK/(nm*lambda0);
lambdap = lambda/lambda0;

% size
[Nd,~,Nk] = size(subK);
No = size(S,1);
Nu = size(Su,1);
[Nm,T] = size(lambda);

% parameters corresponds to measured dofs
Km0p = K0p*S';
subKmp = zeros(Nd,No,Nk);
for i = 1:Nk
    subKmp(:,:,i) = subKp(:,:,i)*S';
end
Mmp = Mp*S';

Mup = [];
Ku0p = [];
subKup = [];
Kupt = [];

% parameters corresponds to unmeasured dofs
if ~isempty(Su)
    Ku0p = K0p*Su';
    subKup = zeros(Nd,Nu,Nk);
    for i = 1:Nk
        subKup(:,:,i) = subKp(:,:,i)*Su';
    end
    Mup = Mp*Su';
end


%% MwG sampler
% Initial values
kappa0 = Theta.prior.kappa.kappa0;
p0 = Theta.prior.kappa.p0;
kappa = kappa0;
kappat = repmat(kappa,1,T);
e0 = Theta.prior.P.e0;
E0 = Theta.prior.P.E0; iE0 = inv(E0);
P = e0*E0; sqrP = chol(P);
alpha = Theta.prior.alpha.u0./Theta.prior.alpha.h0;   % Nm by 1
beta = Theta.prior.beta.r0./Theta.prior.beta.s0;      % Nm by 1

dKmp = bsxfun(@times,subKmp,reshape(kappa,[1 1 Nk]));  % deviation of stiffness
Kmp = Km0p+sum(dKmp,3); % stiffness
dKp = bsxfun(@times,subKp,reshape(kappa,[1 1 Nk]));  % deviation of stiffness
Kp = K0p+sum(dKp,3); % normalized stiffness
Kmpt = repmat(Kmp,1,1,T);
Kpt = repmat(Kp,1,1,T);

if ~isempty(Su)
    dKup = bsxfun(@times,subKup,reshape(kappa,[1 1 Nk]));  % deviation of stiffness
    Kup = Ku0p+sum(dKup,3); % stiffness
    Kupt = repmat(Kup,1,1,T);
end


D = zeros(T,1);   % used in acceptance-rejection
for t = 1:T
    for k = 1:Nm
        K_M_cand = Kp-lambdap(k,t)*Mp;
        [~, S0, ~] = svd(K_M_cand);
        sigma_min = S0(end, end);
        f_lambda =  - 0.5*alpha(k)*sigma_min^2;
        D(t) = D(t)+f_lambda;
    end
end

iter = 1; accept = 0; loglik = -Inf;
Theta.loglik(:,iter) = -loglik;
Theta.posterior.kappa = [];
Theta.posterior.kappa(:,iter) = kappa;
Theta.posterior.kappat(:,:,iter) = kappat;
Theta.posterior.P(:,:,iter) = P;
Theta.posterior.alpha(:,iter) = alpha;
Theta.posterior.beta(:,iter) = beta;


while iter <= OPT.Nsamp
    % p(phi_u|phi_m,lambdap,kappat,beta)
    if ~isempty(Su)
        phi = zeros(Nd,Nm,T);         % samples
        for k = 1:Nm
            for t = 1:T
                Mm_Km = (lambdap(k,t)*Mmp-Kmpt(:,:,t))*phi_m(:,k,t);
                Ku_Mu = Kupt(:,:,t)-lambdap(k,t)*Mup;
                R = qr([Ku_Mu Mm_Km],0);
                R11 = triu(R(1:Nu,1:Nu)); R12 = R(1:Nu,end);
                phi_u = R11\(R12+randn(Nu,1)/sqrt(beta(k)));
                phi(:,k,t) = S'*phi_m(:,k,t) + Su'*phi_u;
                phi(:,k,t) = mynormc(phi(:,k,t));
            end
        end
    else
        phi = phi_m;
    end

    % p(kappat|phi,lambdap,kappa,P,alpha,beta)
    A = zeros(Nd*Nm+Nk,Nk+1);
    A(1:Nk,:) = [sqrP sqrP*kappa];
    for t = 1:T
        for k = 1:Nm
            K_phi = zeros(Nd,Nk);
            for kk = 1:Nk
                K_phi(:,kk) = subKp(:,:,kk)*phi(:,k,t);
            end
            M_phi = (lambdap(k,t)*Mp-K0p)*phi(:,k,t);
            A((Nk+1+(k-1)*Nd):(Nk+k*Nd),:) = [K_phi M_phi]*sqrt(beta(k));
        end
        Rtilt = qr(A,0);
        Rtilt11 = triu(Rtilt(1:Nk,1:Nk)); Rtilt12 = Rtilt(1:Nk,end);
        kappat_cand = Rtilt11\(Rtilt12+randn(Nk,1));

        % acceptance - rejection
        dKp = bsxfun(@times,subKp,reshape(kappat_cand,[1 1 Nk]));  % deviation of stiffness
        Kp_cand = K0p+sum(dKp,3); % normalized stiffness
        D_cand = 0; 
        for k = 1:Nm
            K_M_cand = Kp_cand-lambdap(k,t)*Mp;
            [~, S0, ~] = svd(K_M_cand);
            sigma_min = S0(end, end);
            f_lambda = - 0.5*alpha(k)*sigma_min^2 ;
            D_cand = D_cand+f_lambda;
        end

        if min(0,D_cand-D(t)) > log(rand)
            accept = accept + 1;
            kappat(:,t) = kappat_cand;
            D(t) = D_cand;
            dKmp = bsxfun(@times,subKmp,reshape(kappat_cand,[1 1 Nk]));  % deviation of stiffness
            Kmpt(:,:,t) = Km0p+sum(dKmp,3); % stiffness
            if ~isempty(Su)
                dKup = bsxfun(@times,subKup,reshape(kappat_cand,[1 1 Nk]));  % deviation of stiffness
                Kupt(:,:,t) = Ku0p+sum(dKup,3); % stiffness
            end
            Kpt(:,:,t) = Kp_cand;
        end
    end


    % p(kappa,P|kappat)
    B11 = p0+T;
    B21 = p0*kappa0+sum(kappat,2);
    B22 = p0*(kappa0*kappa0')+kappat*kappat'+iE0;
    L = chol([B11 B21';B21 B22],'lower');
    L11 = L(1,1); L21 = L(2:end,1);
    L22 = L(2:end,2:end); iL22 = inv(L22);
    e = e0 + T;
    P = wishrnd(iL22'*iL22,e,iL22);
    R22 = chol(P);sqrP = R22';
    kappa = (L21+R22\(randn(Nk,1)))/L11;

    % p(alpha|lambdap,kappat) & p(beta|lambdap,phi,kappat)
    u = Theta.prior.alpha.u0+T/2;
    r = Theta.prior.beta.r0+T*(Nd-1)/2;

    h = zeros(Nm,1); s = h;
    for k = 1:Nm
        for t = 1:T
            K_M = Kpt(:,:,t)-lambdap(k,t)*Mp;
            S0 = svd(K_M);
            sigma_min0 = S0(end);
            phi_K_M = K_M*phi(:,k,t);
            squ_phi_K_M = phi_K_M'*phi_K_M;

            if iter <= 100    % annealing
                if sigma_min0 < 1e-2; sigma_min0 = 1e-2;end
                if squ_phi_K_M < 1e-2; squ_phi_K_M = 1e-2; end
            elseif iter <= 300
                if sigma_min0 < 1e-3; sigma_min0 = 1e-3;end
                if squ_phi_K_M < 1e-3; squ_phi_K_M = 1e-3; end
            elseif iter <= 500
                if sigma_min0 < 1e-4; sigma_min0 = 1e-4;end
                if squ_phi_K_M < 1e-4; squ_phi_K_M = 1e-4; end
            elseif iter <= 600
                if sigma_min0 < 1e-5; sigma_min0 = 1e-5;end
                if squ_phi_K_M < 1e-5; squ_phi_K_M = 1e-5; end % the specific annealing scheme can be adjusted.
            end
            h(k) = h(k)+.5*sigma_min0^2;
            s(k) = s(k)+.5*squ_phi_K_M;
        end
    end
    if iter > 0
        alpha = gamrnd(u,1./h);
        beta = gamrnd(r,1./s);
    end


    % log likelihood
    kt_k = bsxfun(@minus,kappat,kappa);
    pdf_alpha = log(gampdf(alpha,u,1./h));
    pdf_beta = log(gampdf(beta,r,1./s));
    loglik = (e0-Nk+T)*sum(log(diag(R22)))-trace(iE0*P)/2+...
        Nk*log(p0)/2-(kappa-kappa0)'*P*(kappa-kappa0)*p0/2-...
        trace(kt_k'*P*kt_k)/2+sum(pdf_alpha)+sum(pdf_beta)+sum(D);

    %     Theta.posterior.phi(:,:,iter) = phi;
    if ~mod(iter,50); fprintf('   iter = %d \n',iter); end
    iter = iter + 1;
    Theta.loglik(:,iter) = loglik;
    Theta.posterior.kappa(:,iter) = kappa;
    Theta.posterior.kappat(:,:,iter) = kappat;
    Theta.posterior.P(:,:,iter) = P;
    Theta.posterior.alpha(:,iter) = alpha;
    Theta.posterior.beta(:,iter) = beta;
end

%% Parallel Multi-Chain
%%
kl = 3400;ber = size(subK,3);
number = 4;
loglikelihood = zeros(1,kl+1,number);kappat2 = zeros(ber,T*kl,4);
kap2 = zeros(ber,kl,4);al2 = zeros(Nm,kl,4);bl2 = zeros(Nm,kl,4);pp2 = zeros(ber,ber*kl,4);

accept = zeros(1,4);
parfor i = 1:4
    kappa1 = kappa;
    xx = randn(ber,1);
    for j = 1:ber
        kappa1(j) = kappa1(j)*(1+0.05*xx(j));
    end
    kappat1 = kappat;D1 = D;Kpt1 = Kpt;llh = zeros(1,kl+1);alpha1 = alpha;
    beta1 = beta;Kmpt1 = Kmpt; kap = zeros(ber,kl+1);kapt = zeros(ber,T,kl+1);
    pp = zeros(ber,ber,kl+1); al = zeros(Nm,kl+1);bl = zeros(Nm,kl+1);sqrP1 = sqrP;
    iter = 1;accept2(i) = 0;en_avg = 0;alpha_en = 0;Mup1 = [];Kupt1 = [];subKup1 = [];
    if ~isempty(Su)
        Mup1 = Mup;Kupt1 = Kupt;subKup1 = subKup;
    end

    if flag
        en_avg = 0.2; alpha_en = 0.3;gamma0 = 2.38/sqrt(Nk);
    else
        gamma0 = 1;
    end

    while iter <= kl
        % p(phi_u|phi_m,lambdap,kappat,beta)
        if ~isempty(Su)
            phi = zeros(Nd,Nm,T);         % samples
            for k = 1:Nm
                for t = 1:T
                    Mm_Km = (lambdap(k,t)*Mmp-Kmpt1(:,:,t))*phi_m(:,k,t);
                    Ku_Mu = Kupt1(:,:,t)-lambdap(k,t)*Mup1;
                    R = qr([Ku_Mu Mm_Km],0);
                    R11 = triu(R(1:Nu,1:Nu)); R12 = R(1:Nu,end);
                    phi_u = R11\(R12+randn(Nu,1)/sqrt(beta1(k)));
                    phi(:,k,t) = S'*phi_m(:,k,t) + Su'*phi_u;
                    phi(:,k,t) = mynormc(phi(:,k,t));

                end
            end
        else 
            phi = phi_m;
        end

        A = zeros(Nd*Nm+Nk,Nk+1);
        A(1:Nk,:) = [sqrP1 sqrP1*kappa1];
        %     kappat(:,iter+1) = zeros(Nk,T);
        for t = 1:T
            for k = 1:Nm
                K_phi = zeros(Nd,Nk);
                for kk = 1:Nk
                    K_phi(:,kk) = subKp(:,:,kk)*phi(:,k,t);
                end
                M_phi = (lambdap(k,t)*Mp-K0p)*phi(:,k,t);
                A((Nk+1+(k-1)*Nd):(Nk+k*Nd),:) = [K_phi M_phi]*sqrt(beta1(k));
            end
            Rtilt = qr(A,0);
            Rtilt11 = triu(Rtilt(1:Nk,1:Nk)); Rtilt12 = Rtilt(1:Nk,end);
            kappat_cand = Rtilt11\(Rtilt12+gamma0*randn(Nk,1));

            % acceptance - rejection
            dKp = bsxfun(@times,subKp,reshape(kappat_cand,[1 1 Nk]));  % deviation of stiffness
            Kp_cand = K0p+sum(dKp,3); % normalized stiffness
            D_cand = 0;
            nfora = 1.1; % numeric protection for alpha
            alpha1 = min(nfora*alpha,alpha1); % numerical protection strategy can be adjusted according to specific cases

            for k = 1:Nm
                K_M_cand = Kp_cand-lambdap(k,t)*Mp;
                [~, S0, ~] = svd(K_M_cand);
                sigma_min = S0(end, end);
                f_lambda =  - 0.5*(alpha1(k))*sigma_min^2;
                D_cand = D_cand+f_lambda;
            end
            if min(0,D_cand-D1(t)) > log(rand)
                accept(i) = accept(i) + 1;
                kappat1(:,t) = kappat_cand;
                D1(t) = D_cand;
                dKmp = bsxfun(@times,subKmp,reshape(kappat_cand,[1 1 Nk]));  % deviation of stiffness
                Kmpt1(:,:,t) = Km0p+sum(dKmp,3); % stiffness
                if ~isempty(Su)
                    dKup = bsxfun(@times,subKup1,reshape(kappat_cand,[1 1 Nk]));  % deviation of stiffness
                    Kupt1(:,:,t) = Ku0p+sum(dKup,3); % stiffness
                end
                Kpt1(:,:,t) = Kp_cand;
            end
        end

        if flag
            if mod(iter,100) == 0&&iter<1000
                en = (accept(i)-accept2(i))/100/T;
                en_avg = (1 - alpha_en) * en_avg + alpha_en * en;
                gamma0 = gamma0 * exp(0.5 * (en_avg - 0.25));
                gamma0 = min(max(gamma0, 1e-3), 10);
                accept2(i) = accept(i);
            end
        end

        % p(kappa,P|kappat)
        B11 = p0+T;
        B21 = p0*kappa0+sum(kappat1,2);
        B22 = p0*(kappa0*kappa0')+kappat1*kappat1'+iE0;
        L = chol([B11 B21';B21 B22],'lower');
        L11 = L(1,1); L21 = L(2:end,1);
        L22 = L(2:end,2:end); iL22 = inv(L22);
        e = e0 + T;
        P1 = wishrnd(iL22'*iL22,e,iL22);
        R22 = chol(P1);sqrP1 = R22';
        kappa1 = (L21+R22\(randn(Nk,1)))/L11;

        % p(alpha|lambdap,kappat) & p(beta|lambdap,phi,kappat)
        u = Theta.prior.alpha.u0+T/2;
        r = Theta.prior.beta.r0+T*(Nd-1)/2;
        h = zeros(Nm,1); s = h;
        for k = 1:Nm
            for t = 1:T
                K_M = Kpt1(:,:,t)-lambdap(k,t)*Mp;
                [~, S0, ~] = svd(K_M);
                sigma_min0 = S0(end, end);
                phi_K_M = K_M*phi(:,k,t);
                squ_phi_K_M = phi_K_M'*phi_K_M;
                h(k) = h(k)+.5*(sigma_min0^2);
                s(k) = s(k)+.5*squ_phi_K_M;
            end
        end
        alpha1 = gamrnd(u,1./h);
        beta1 = gamrnd(r,1./s);

        % log likelihood
        kt_k = bsxfun(@minus,kappat1,kappa1);
        pdf_alpha = log(gampdf(alpha1,u,1./h));
        pdf_beta = log(gampdf(beta1,r,1./s));
        loglik1 = (e0-Nk+T)*sum(log(diag(R22)))-trace(iE0*P1)/2+...
            Nk*log(p0)/2-(kappa1-kappa0)'*P1*(kappa1-kappa0)*p0/2-...
            trace(kt_k'*P1*kt_k)/2+sum(pdf_alpha)+sum(pdf_beta)+sum(D1);

        iter = iter + 1;
        llh(:,iter) = loglik1;
        kap(:,iter) = kappa1;
        kapt(:,:,iter) = kappat1;
        pp(:,:,iter) = P1;
        al(:,iter) = alpha1;
        bl(:,iter) = beta1;
        if ~mod(iter,50); fprintf('   iter = %d \n',iter); end
    end

    loglikelihood(:,:,i) = llh;kap = kap(:,2:end);al = al(:,2:end);bl = bl(:,2:end);
    kapt = kapt(:,:,2:end); kapt1 = zeros(ber,T*kl);pp = pp(:,:,2:end);pp1 = zeros(ber,ber*kl);
    for j = 1:kl
        kapt1(:,T*(j-1)+1:T*j) = kapt(:,:,j);pp1(:,ber*(j-1)+1:ber*j) = pp(:,:,j);
    end
    kappat2(:,:,i) = kapt1;kap2(:,:,i) = kap;pp2(:,:,i) = pp1;al2(:,:,i) = al;bl2(:,:,i) = bl;
end
for i = 1:number
    loglikelihood(:,1,i) = loglik;
end
a3 = Theta.loglik(1,2:600);
for i = 1:4
    a4(i,:) = [a3 loglikelihood(:,:,i)];
end
end