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
percentInfo = 0.99;

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
  super_U = zeros(sup_m-42,n);
  for i=1:m
    if (mod(i,3) < 1)
      super_U(3*(i-1) + 1,:) = U(i,:);
    end
  end
  U = [zeros(42,n);super_U];
  time=toc;
  converged=size(d,1);
  condition=1-d(end)/d(1);
  fprintf(['%d sing. values found in %7.3f seconds ; they ' ...
           'provide %3.2f%% variability. \n'], converged, time,condition);

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

%%%% Display %%%%
global Lx Ly Nx Ny;

% draw result
    x = linspace(0,Lx,Nx);     %  Independent variable x
    y = linspace(0,Ly,Ny);     %  Independent variable y
    [Mx, My] = meshgrid(x,y);  %  2D arrays, mainly for plotting
    Mx = Mx'; My = My';        %  MatLab is strange!

      figure(2)
 for tt=1:nt
          subplot(1,2,1);
          z = X((tt-1)*ns+1:tt*ns,1);
          z = reshape(z,Nx,Ny);
          surf(Mx,My,z); shading('interp');
          axis([0,Lx,0,Ly ,5000,6000]);
	  pbaspect([3 1 3])
	  title('Solution')

          subplot(1,2,2);
          zappr = Xp((tt-1)*ns+1:tt*ns,1);
          zappr = reshape(zappr,Nx,Ny);
          surf(Mx,My,zappr); shading('interp');
	  axis([0,Lx,0,Ly ,5000,6000]);
	  pbaspect([3 1 3])
	  title('Reconstruction')
         drawnow
  end
