function [Pi,dim] = Classification_computation(Nens, percentInfo, use_svd, Fobs)

  % Initiatlization
  Pi = zeros(3, 1);
  dim = zeros(3, 1);

  for GWi = 1:3
    % Generation of the data set
    Fi = Model(GWi, Nens);

    % Computation of the mean and anomalies
    mFi = mean(Fi, 2);
    Zi = Fi-repmat(mFi, 1, Nens);

    % Computation of the dominant left eigenvectors of Zi accordingly to the
    % targeted quality of the subspace approximation
    if (use_svd)
      % Computation of the singular value decomposition
      [Ui,Si,~] = svd(Zi, 0);

      % Computation of the column vector of the main diagonal elements of A
      Di = diag(Si);

      % Verification if Di (Zi) is not null
      if (Di(1) == 0)
        disp('Alert: the matrix is null')
        return
      end

      % Compute the number of dominant eigenvalues
      k=2;
      while (Di(k)/Di(1) > 1-percentInfo)
        k = k+1;
      end
      k = k-1;
      dim(GWi) = k;

      % Only keep dominant left eigenvectors
      Ui = Ui(:,1:k);
    else
      % TODO: power iteration method
    end

    % Computation of the anomaly vector of observations
    Zobsi = Fobs - mFi;

    % Computation of Pi = ||(I - Ui.Uit).Zobsi||, distance of Zobsi to the
    % subspace associated to the dominant singular values
    % Computation of the projection of Zobsi on the subspace
    Intermediate = Ui.' * Zobsi;
    Projected = Ui * Intermediate;

    % Computation of the distance of Zobsi to the subspace
    Pi(GWi) = norm(Zobsi - Projected);
  end
end
