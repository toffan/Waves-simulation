clear all;
close all;

%%% Set check to true for validation %%%
check=true;

% Number of simulations
% TODO : assess the impact fot this parameter on the results
% of the reconstruction
Nens = 50;

% Wind parameter
GW = 1;
% Stopping criterions
maxIter=100;
% TODO : assess the impact of this parameter on the results
% of the reconstruction
percentInfo = 0.95;

% Generate the simulations
F = Model(GW,Nens);

% Ensemble mean
muF = mean(F,2);
% Compute the anomaly matrix
Z  = F - repmat(muF,1,Nens);

%%%%%%%  Compute the SVD of A   %%%%%%%
if (check)
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
    fprintf('dimension of the subspace: %d\n',converged);
else
    p = 1;
    epsilon = 1e-6;
    [m,n] = size(Z);
    V = zeros(m,n);
    for i=1:n
        V(i,i) = 1;
    end
    niter = 0;
    converged = 0;
    PercentReached = 0;
    normeA = norm(Z*Z.');
    while PercentReached < percentInfo || niter < MaxIter
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
        [X,Gamma] = sortem(H);
        V = V*X;
        for i=converged+1:n
            if (norm(Z*(Z.')*V(:,i) - Gamma(i,i)*V(:,i))/normeA <= epsilon)
                converged = converged + 1;
                Zk = zeros(m,n);
                for j=1:i
                    Zk = Zk + nthroot(Gamma(j,j),4)*V(:,i)*((V(:,i).')*(Z.'));
                end
                PercentReached = 1 - (norm(Z - Zk)/sqrt(normeA));
            else
                break;
            end
        end
        niter = niter + 1;
    end
    U = V(:,1:converged);
end

%%%%%%%       Reconstruction        %%%%%%%
[X, ns, nt] = Model(GW,1);
X0 = X(1:ns,:);
%%%%

% reconstruction of X with X0 using the steepest descent method to solve:
% U0t*U0*α = U0t*Z0 (normal equation)
Z0 = X0 - muF(1:ns,:); % initial anomaly vector
epsilon = 1e-6;        % precision for α computation
U0 = U(1:ns,:);
U0t = transpose(U0);
Usdp = U0t*U0;
Z0p = U0t*Z0;
alpha = zeros(converged,1); % initialisation
r = Z0p - Usdp*alpha;
nr = norm(r);
pl = 0;
while (nr / norm(Z0) > epsilon)
  lambda = nr^2 / (norm(U0*r)^2);
  alpha = alpha + lambda*r;
  r = Z0p - Usdp*alpha;
  nr = norm(r);
  pl = pl + 1;
end
fprintf('nombre itérations : %d\n', pl);
Zp = U*alpha;
%%%%

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
