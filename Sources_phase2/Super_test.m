%%% Set check to true for validation %%%
check=false;

% Number of simulations
%%%
% TODO : assess the impact fot this parameter on the results
% of the reconstruction
Nens = 50;

% Wind parameter
GW = 1;

% TODO : assess the impact of this parameter on the results
% of the reconstruction
percentInfo = 0.95;

mean_time = 0;
mean_error = 0;

for i=1:100
    % Generate the simulations
    F = Model(GW,Nens);

    % Ensemble mean
    muF = mean(F,2);
    % Compute the anomaly matrix
    Z   = F - repmat(muF,1,Nens);

    %%%%%%%  Compute the SVD of A    %%%%%%%
    if (check)
        tic;
        [U,S,~] = svd(Z,0);
        D = diag(S);
        if (D(1)==0)
            disp('Alert: the matrix is null')
            return
        end
        converged=1;
        while (D(converged)/D(1)>1-percentInfo) 
            converged=converged+1;
        end
        converged=converged-1;
        U = U(:,1:converged);
        time = toc;
        fprintf('dimension of the subspace: %d\n temps: %7.3f seconds',converged,time);
    else
        addpath('Fortran') 
        p=1;
        m=10;
        maxit=300;
        eps=1.e-8;
        [sup_m,sup_n] = size(F);
        tic;
        [U,d] = fortran_subspace_iter_sv(Z(43:46:sup_m,:),m,p,percentInfo,eps,maxit);
  % [U,d] = fortran_subspace_iter_sv(Z,m,p,percentInfo,eps,maxit);
        [m,n] = size(U);
        U = [U;zeros(sup_m-m,n)];
        time=toc;
        converged=size(d,1);
        condition=1-d(end)/d(1);
        % fprintf(['%d sing. values found in %7.3f seconds ; they ' ...
        %           'provide %3.2f%% variability. \n'], converged, time,condition);
 
    end


    %%%%%%%       Reconstruction        %%%%%%%
    [X, ns, nt] = Model(GW,1);
    X0 = X(1:ns,:);
    %%%%
    Zp = zeros(size(X)); 
    %%%%
    Z0=X0-muF(1:ns);
    alpha=(U(1:ns,:)'*U(1:ns,:))\(U(1:ns,:)'*Z0);
    Zp=U*alpha;
    
    %%%% Compute the error %%%%
    Xp = Zp + muF;
    error=norm(Xp-X)/norm(X);
    fprintf('error = %f\n',error);

    mean_time = mean_time + time;
    mean_error = mean_error + error;
    
end

mean_time = mean_time/100
mean_error = mean_error/100

fprintf('mean_error = %f\n', mean_error);
fprintf('mean_time = %7.3f',mean_time);