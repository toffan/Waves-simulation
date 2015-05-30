function [] = Reco_param2(n_test, Nens_i, percentInfo_i, Nens, percentInfo)

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
    U1 = U(:,1:converged);
    time1 = toc;

    %%% Fortran algorithm

    addpath('Fortran')
    p=1;
    m=10;
    maxit=300;
    eps=1.e-8;
    [sup_m,sup_n] = size(F);
    tic;
    [U2,d] = fortran_subspace_iter_sv(Z(43:46:sup_m,:),m,p,percentInfo,eps,maxit);
    [m,n] = size(U2);
    super_U = zeros(sup_m-42, n);
    for i=1:m
      if (mod(i,3) < 1)
        super_U(3*(i-1) + 1,:) = U2(i,:);
      end
    end
    U2 = [zeros(42,n) ; super_U];
    time2=toc;

    local_results(i,1) = (time2 - time1)/time1;
    local_local_results = zeros(n_test, 1);

    for j=1:n_test

      %%%%%%%       Reconstruction        %%%%%%%
      [X, ns, nt] = Model(GW,1);
      X0 = X(1:ns,:);
      %%%%
      %Zp = zeros(size(X));
      %%%%
      Z0=X0-muF(1:ns);

      %%% svd
      alpha=(U1(1:ns,:)'*U1(1:ns,:))\(U1(1:ns,:)'*Z0);
      Zp=U1*alpha;
      %%%% Compute the error %%%%
      Xp = Zp + muF;
      error1=norm(Xp-X)/norm(X);

      %%% fortran
      alpha=(U2(1:ns,:)'*U2(1:ns,:))\(U2(1:ns,:)'*Z0);
      Zp=U2*alpha;
      %%%% Compute the error %%%%
      Xp = Zp + muF;
      error2=norm(Xp-X)/norm(X);

      %%% Save results
      local_local_results(j,1) = (error2 - error1)/error1;

    end

    local_results(i,2) = mean(local_local_results);

  end

  line = mean(local_results);

  %%% Save infos %%%
  try
    load('resultats2.mat');
  catch fnf
    resultats2 = zeros(13, 13, 2);
  end

  resultats2(Nens_i,percentInfo_i,:) = line;

  save('resultats2.mat', 'resultats2');
end
