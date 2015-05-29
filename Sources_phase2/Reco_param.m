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
    U1 = U(:,1:converged);
    time1 = toc;

    %%% Fortran algorithm

    addpath('Fortran')
    p=1;
    m=10;
    maxit=300;
    eps=1.e-8;
    tic;
    [U2,d] = fortran_subspace_iter_sv(Z,m,p,percentInfo,eps,maxit);
    %converged2=size(d,1);
    time2=toc;
    %condition=1-d(end)/d(1);

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
    load('resultats.mat');
  catch fnf
    resultats = zeros(13, 13, 2);
  end

  resultats(Nens_i,percentInfo_i,:) = line;

  save('resultats.mat', 'resultats');
end
