clear all;
close all;

number = 5;
Nens = [5:3:77]; Nens = Nens(:);
percentInfo = [0.87:0.005:0.99]; percentInfo = percentInfo(:);
load('resultats.mat');

figure(1)

subplot(2, 1, 1);
surf(percentInfo(1:number), Nens(1:number), resultats(:,:,1)); shading('interp');
axis([percentInfo(1), percentInfo(number), Nens(1), Nens(number), min(min(resultats(:,:,1))), 1 + max(max(resultats(:,:,1)))]);
title('time')

subplot(2, 1, 2);
surf(percentInfo(1:number), Nens(1:number), results(:,:,2)); shading('interp');
axis([percentInfo(1), percentInfo(number), Nens(1), Nens(number), min(min(resultats(:,:,2))), 1 + max(max(resultats(:,:,2)))]);
title('error')
