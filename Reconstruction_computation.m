function [rel_error,dim,Xp] = Reconstruction_computation(Nens, percentInfo, use_svd, X, ns, GW)

  % Generate the simulations
  F = Model(GW, Nens);

  % Ensemble mean
  muF = mean(F, 2);
  % Compute the anomaly matrix
  Z = F - repmat(muF, 1, Nens);

  % Computation of the dominant left eigenvectors of Z accordingly to the
  % targeted quality of the subspace approximation
  if (use_svd)
    % Computation of the singular value decomposition
    [U,S,~] = svd(Z, 0);

    % Computation of the column vector of the main diagonal elements of A
    D = diag(S);

    % Verification if D (Z) is not null
    if (D(1) == 0)
      disp('Alert: the matrix is null')
      return
    end

    % Compute the number of dominant eigenvalues
    k=2;
    while (D(k)/D(1) > 1-percentInfo)
      k = k+1;
    end
    k = k-1;

    % Only keep dominant left eigenvectors
    U = U(:,1:k);
  else
    [U,k] = Power_iteration(Z, percentInfo);
  end
  dim = k;

  % reconstruction of X with X0 using the steepest descent method to solve:
  % U0t*U0*α = U0t*Z0 (normal equation)
  X0 = X(1:ns,:);        % initial vector
  Z0 = X0 - muF(1:ns,:); % initial anomaly vector
  epsilon = 1e-14;       % precision for α computation
  U0 = U(1:ns,:);
  U0t = U0.';
  Usdp = U0t*U0;
  Z0p = U0t*Z0;
  alpha = zeros(k,1); % initialisation
  r = Z0p - Usdp*alpha;
  p = 0;
  while (norm(r) / norm(Z0p) > epsilon)
    lambda = (r.'*r) / (r.'*Usdp*r);
    alpha = alpha + lambda*r;
    r = Z0p - Usdp*alpha;
    p = p+1;
  end
  Zp = U*alpha;
  %%%%
  
  %%%% Compute the error %%%%
  Xp = Zp + muF;
  rel_error=norm(Xp-X)/norm(X);

end
