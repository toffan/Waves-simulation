clear all;
close all;

% Use matlab svd function instead of power iteration method
use_svd = true;
% Basic try / parameter impact
basic = true;
% Number of simulations
Nens = 50;
% Wind parameter
GW = 1;
% Quality of the subspace approximation
percentInfo = 0.95;

if (basic)
  % Reconstruction
  [X, ns, nt] = Model(GW, 1);

  [rel_error,dim,Xp] = Reconstruction_computation(Nens, percentInfo, use_svd, X, ns, GW);
  
  fprintf('error = %f\n',rel_error);
  
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
end
