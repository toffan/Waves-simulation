function [U,k] = Power_iteration(Z, percentInfo)

  % Stopping criterions
  maxIter=100;
  epsilon = 1e-3;

  % Initialization
  k = 0;

  p = 1;
  [m,n] = size(Z);
  V = zeros(m,n);
  for i=1:n
    V(i,i) = 1;
  end
  niter = 0;
  PercentReached = 0;
  normeA = norm(Z.'*Z);
  Gamma = zeros(n);
  VapTrouvees = zeros(1,n);
  while k == 0 || (sqrt(VapTrouvees(1,k) / VapTrouvees(1,1)) > 1 - percentInfo && niter < maxIter)
    for i=1:p
      V = (Z.')*V;
      V = Z*V;
    end
    for i=1:n
      for j=1:(i-1)
        V(:,i) = V(:,i) - (dot(V(:,i),V(:,j))/norm(V(:,j)))*V(:,j);
      end
      V(:,i) = (1/norm(V(:,i)))*V(:,i);
    end
    H = (Z.')*V;
    H = (H.')*H;
    [X,Gamma] = eig(H);
    [tri,indices] = sort(diag(Gamma),'descend');
    Gamma = diag(tri);
    X = X(:,indices);
    V = V*X;
    for i=k+1:n
      temp = (Z.')*V(:,i);
      if (norm(Z*temp - Gamma(i,i)*V(:,i))/normeA <= epsilon)
        k = k+1;
        VapTrouvees(1,i) = Gamma(i,i);
        break;
      else
        break;
      end
    end
    niter = niter + 1;
  end
  k = k-1;
  U = V(:,1:k);

end
