clear all;
close all;

Nens = [20:10:110 125 150 200]; Nens = Nens(:);
percentInfo = [0.87:0.01:0.99]; percentInfo = percentInfo(:);
load('resultats2.mat');

figure(1)

subplot(1, 2, 1);
surf(percentInfo(1:end), Nens(1:end), resultats2(:,:,1)); shading('interp');
axis([percentInfo(1), percentInfo(end), Nens(1), Nens(end), min(min(resultats2(:,:,1))) - 1, 1 + max(max(resultats2(:,:,1)))]);
title('time')

subplot(1, 2, 2);
surf(percentInfo(1:end), Nens(1:end), resultats2(:,:,2)); shading('interp');
axis([percentInfo(1), percentInfo(end), Nens(1), Nens(end), min(min(resultats2(:,:,2))) - 1, 1 + max(max(resultats2(:,:,2)))]);
title('error')
