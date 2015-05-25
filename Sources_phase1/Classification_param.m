function [] = Classification_param(number, Nens_i, percentInfo_i, Nens, percentInfo, use_svd)
  load('observation-G1.mat');

  local_results = zeros(number, 3);
  for n_try=1:number
    [Pi,dim] = Classification_computation(Nens, percentInfo, use_svd, Fobs);
    dim = mean(dim);
    [Pi,GWi] = min(Pi);
    local_results(n_try,:) = [GWi, Pi, dim];
  end
  if (local_results(:,1) == repmat(local_results(1,1), size(local_results(:,1),1), 1))
    load('results.mat');
    fprintf('%d, %d\n', Nens, percentInfo);
    results
    results(Nens_i,percentInfo_i,:) = mean(local_results);
    results
    save('results.mat', 'results');
  end
end
