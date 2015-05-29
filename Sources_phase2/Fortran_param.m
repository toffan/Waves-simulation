function [] = Fortran_param(n_test, p_i, m_i, p, m)

  fprintf('tout debut');

  local_results = zeros(n_test, 2);

  % Wind parameter
  GW = 1;
  Nens = 50;
  percentInfo = 0.95;

  for i=1:n_test

    fprintf('debut');

    % Generate the simulations
    F = Model(GW,Nens);

    fprintf('génération effectuée');

    % Ensemble mean
    muF = mean(F,2);
    % Compute the anomaly matrix
    Z   = F - repmat(muF,1,Nens);

    %%%%%%%  Compute the SVD of A    %%%%%%%
    addpath('Fortran')
    maxit=300;
    eps=1.e-8;
    tic;
    [U2,d] = fortran_subspace_iter_sv(Z,m,p,percentInfo,eps,maxit);
    %converged2=size(d,1);
    time=toc;
    %condition=1-d(end)/d(1);

    local_results(i,1) = time;
    local_local_results = zeros(n_test, 1);

    for j=1:n_test

      %%%%%%%       Reconstruction        %%%%%%%
      [X, ns, nt] = Model(GW,1);
      X0 = X(1:ns,:);
      %%%%
      %Zp = zeros(size(X));
      %%%%
      Z0=X0-muF(1:ns);

      %%% fortran
      alpha=(U2(1:ns,:)'*U2(1:ns,:))\(U2(1:ns,:)'*Z0);
      Zp=U2*alpha;

      %%%% Compute the error %%%%
      Xp = Zp + muF;
      error=norm(Xp-X)/norm(X);

      %%% Save results
      local_local_results(j,1) = error;

    end

    fprintf('fin');
    local_results(i,2) = mean(local_local_results);

  end

  line = mean(local_results);

  %%% Save infos %%%
  try
    load('resultats_f.mat');
  catch fnf
    resultats_f = zeros(9, 9, 2);
  end

  resultats_f(p_i,m_i,:) = line;

  save('resultats_f.mat', 'resultats_f');
end
