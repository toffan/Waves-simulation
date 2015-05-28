function [] = Reco_param(n_test, Nens_i, percentInfo_i, Nens, percentInfo)

  local_results = zeros(n_test, 2);

  % Wind parameter
  GW = 1;

  for i=1:n_test

    % Generate the simulations
    F = Model(GW,Nens);

    % Ensemble mean
    muF = mean(F,2);
    % Compute the anomaly matrix
    Z   = F - repmat(muF,1,Nens);

    %%%%%%%  Compute the SVD of A    %%%%%%%
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

    local_results(i,1) = time;
    local_local_results1 = zeros(n_test, 1);

    for j=1:n_test

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

      %%% Save results
      local_local_results(j,:) = [d_time d_error]

    end

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
    converged=size(d,1);
    time=toc;
    condition=1-d(end)/d(1);

    local_results(i,1) = time;
    local_local_results2 = zeros(n_test, 1);

    for j=1:n_test

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

      %%% Save results
      local_local_results(j,:) = [d_time d_error]

    end

    local_local_results = local_local_results2 - local_local_results1;

    local_results(i,2) = mean(local_local_results);

  end

  line = mean(local_results);

  %%% Save infos %%%
  try
    load('resultats.mat');
  catch fnf
    % 25 Nens x 25 percentInfo
    resultats = zeros(25, 25, 2);
  end

  resultats(Nens_i,percentInfo_i,:) = line;

  save('resultats.mat', 'resultats');
end
