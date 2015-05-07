clear all;
close all;

% Use matlab svd function instead of power iteration method
use_svd = true;
% Basic try / parameter impact
basic = true;
% Number of simulations
Nens = 20;
% Quality of the subspace approximation
percentInfo = 0.95;

% Reading the observations
load('observation-G1.mat');

if (basic)
  [Pi,~] = Classification_computation(Nens, percentInfo, use_svd, Fobs);
  figure(1)
  bar(Pi)
else
  % Inatialization of parameters
  number = 3; % max 5
  Nens = [5:20:95];
  Nens = Nens(:);
  percentInfo = [0.90:0.02:0.99];
  percentInfo = percentInfo(:);

  % Initialization of results matrix
  results = zeros(number, number, 3);

  for Nens_i=1:number
    for percentInfo_i=1:number
      local_results = zeros(number, 3);
      for n_try=1:number
        time = cputime;
        [Pi,dim] = Classification_computation(Nens(Nens_i), percentInfo(percentInfo_i), use_svd, Fobs);
        time = cputime - time
        dim = mean(dim);
        Pi = min(Pi);
        local_results(n_try,:) = [Pi dim time];
      end
      results(Nens_i,percentInfo_i,:) = mean(local_results);
    end
  end

  % Draw results
  figure(1)

  subplot(2, 2, 1);
  surf(Nens(1:number), percentInfo(1:number), results(:,:,1)); shading('interp');
  axis([Nens(1), Nens(number), percentInfo(1), percentInfo(number), min(min(results(:,:,1))), max(max(results(:,:,1)))]);
  title('Pi')

  subplot(2, 2, 2);
  surf(Nens(1:number), percentInfo(1:number), results(:,:,2)); shading('interp');
  axis([Nens(1), Nens(number), percentInfo(1), percentInfo(number), min(min(results(:,:,2))), max(max(results(:,:,2)))]);
  title('dim')

  subplot(2, 3, 5);
  surf(Nens(1:number), percentInfo(1:number), results(:,:,3)); shading('interp');
  axis([Nens(1), Nens(number), percentInfo(1), percentInfo(number), min(min(results(:,:,3))), max(max(results(:,:,3)))]);
  title('time')

end
