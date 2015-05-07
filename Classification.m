clear all;
close all;

% Use matlab svd function instead of power iteration method
use_svd = true;
% Basic try / parameter impact
basic = false;
% Number of simulations
Nens = 20;
% Quality of the subspace approximation
percentInfo = 0.95;

if (basic)
  % Reading the observations
  load('observation-G1.mat');

  [Pi,~] = Classification_computation(Nens, percentInfo, use_svd, Fobs);
  figure(1)
  bar(Pi)
else
  number = 5;
  % results.mat contain series of number tests for parameters
  Nens = [5:15:65]; Nens = Nens(:);
  percentInfo = [0.90:0.02:0.98]; percentInfo = percentInfo(:);
  % results = zeros(number, number, 3); has to be executed before
  load('results.mat');

  % Draw results
  figure(1)

  subplot(2, 2, 1);
  surf(percentInfo(1:number), Nens(1:number), results(:,:,1)); shading('interp');
  axis([percentInfo(1), percentInfo(number), Nens(1), Nens(number), min(min(results(:,:,1))), 1 + max(max(results(:,:,1)))]);
  title('GWi')

  subplot(2, 2, 2);
  surf(percentInfo(1:number), Nens(1:number), results(:,:,2)); shading('interp');
  axis([percentInfo(1), percentInfo(number), Nens(1), Nens(number), min(min(results(:,:,2))), 1 + max(max(results(:,:,2)))]);
  title('Pi')

  subplot(2, 2, 3);
  surf(percentInfo(1:number), Nens(1:number), results(:,:,3)); shading('interp');
  axis([percentInfo(1), percentInfo(number), Nens(1), Nens(number), min(min(results(:,:,3))), 1 + max(max(results(:,:,3)))]);
  title('dim')
end
