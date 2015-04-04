
%Number of simulations
Nens=20;
%Quality of the subspace approximation
percentInfo=0.95;

%Reading the observations
%%%
% TODO: select the file associated to your group number
load('observation-G1.mat');

% Initiatlization
Pi=zeros(3,1);

for GWi = 1:3
    %Generation of the data set
    Fi = Model(GWi,Nens);

    %Computation of the mean and anomalies
    mFi= mean(Fi,2);
    Zi=Fi-repmat(mFi,1,Nens);

    %%%
    % TODO: compute the dominant left eigenvectors of Zi
    % accordingly to the targeted quality of the subspace approximation

    %%%
    % TODO: compute Pi
    %
end

figure(1)
bar(Pi)



