clear all;
close all;

Nens = [20:10:110 125 150 200]; Nens = Nens(:);
percentInfo = [0.87:0.01:0.99]; percentInfo = percentInfo(:);
load('resultats.mat');

figure(1)

subplot(1, 2, 1);
surf(percentInfo(1:end), Nens(1:end), resultats(:,:,1)); shading('interp');
axis([percentInfo(1), percentInfo(end), Nens(1), Nens(end), min(min(resultats(:,:,1))) - 1, 1 + max(max(resultats(:,:,1)))]);
title('time')

subplot(1, 2, 2);
surf(percentInfo(1:end), Nens(1:end), resultats(:,:,2)); shading('interp');
axis([percentInfo(1), percentInfo(end), Nens(1), Nens(end), min(min(resultats(:,:,2))) - 1, 1 + max(max(resultats(:,:,2)))]);
title('error')
