clear all;
close all;

p = [1:1:9]; p = p(:);
m = [5:3:29]; m = m(:);
load('resultats_f.mat');

figure(1)

subplot(1, 2, 1);
surf(m(1:end), p(1:end), resultats_f(:,:,1)); shading('interp');
axis([m(1), m(end), p(1), p(end), min(min(resultats_f(:,:,1))) - 1, 1 + max(max(resultats_f(:,:,1)))]);
title('time')

subplot(1, 2, 2);
surf(m(1:end), p(1:end), resultats_f(:,:,2)); shading('interp');
axis([m(1), m(end), p(1), p(end), min(min(resultats_f(:,:,2))) - 1, 1 + max(max(resultats_f(:,:,2)))]);
title('error')
