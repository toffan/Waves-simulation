function [] = Classification_param(number, Nens_i, percentInfo_i, Nens, percentInfo, use_svd)
  load('observation-G1.mat');
  load('results.mat');

  local_results = zeros(number, 3);
  for n_try=1:number
    [Pi,dim] = Classification_computation(Nens, percentInfo, use_svd, Fobs))
    dim = mean(dim);
    Pi = min(Pi);
    local_results(n_try,:) = [GWi, Pi, dim];
  end
  if (local_results(:,1) == repmat(local_results(1,1), size(local_results(:,1)), 1))
    results(Nens_i,percentInfo_i,:) = mean(local_results);
  end

  save('results.mat', results);
end
