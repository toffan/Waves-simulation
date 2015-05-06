clear all;
close all;

%Use matlab svd function instead of power iteration method
use_svd = true;
%Number of simulations
Nens=20;
%Quality of the subspace approximation
percentInfo=0.95;

%Reading the observations
load('observation-G1.mat');

% Initiatlization
Pi=zeros(3,1);

for GWi = 1:3
    %Generation of the data set
    Fi = Model(GWi,Nens);

    %Computation of the mean and anomalies
    mFi= mean(Fi,2);
    Zi=Fi-repmat(mFi,1,Nens);

    %Computation of the dominant left eigenvectors of Zi accordingly to the
    %targeted quality of the subspace approximation
    if (use_svd)
      %Computation of the singular value decomposition
      [Ui,Si,~] = svd(Zi,0);
      %Computation of the column vector of the main diagonal elements of A
      Di = diag(Si);
      %Verification if Di (Zi) is not null
      if (Di(1) == 0)
        disp('Alert: the matrix is null')
        return
      end
      k=2;
      %Compute the number of dominant eigenvalues
      while (Di(k)/Di(1) > 1-percentInfo)
        k = k+1;
      end
      k = k-1;
      %Only keep dominant left eigenvectors
      Ui = Ui(:,1:k);
      fprintf('dimension of the subspace: %d\n',k);
    else
      % TODO: power iteration method
    end

    %Computation of the anomaly vector of observations
    Zobsi = Fobs - mFi;
    %Computation of Pi = ||(I - Ui.Uit).Zobsi||, distance of Zobsi to the
    %subspace associated to the dominant singular values
    %Computation of the projection of Zobsi on the subspace
    Uit = transpose(Ui);
    Intermediate = Uit*Zobsi;
    Projected = Ui*Intermediate;
    %Computation of the distance of Zobsi to the subspace
    Pi(GWi) = norm(Zobsi - Projected);
end

figure(1)
bar(Pi)
