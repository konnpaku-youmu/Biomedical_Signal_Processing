clear;
clc;

load ex1data.mat;

A = 1/(sqrt(5)) * [1, -2; -2, 1];

sigma = 0.05;
s = rand(2, 800) * 2 - 1;
[x, n] = noisy(A*s, 10);

SNR = 10*log10(frob(A*s).^2 / frob(n).^2);

%% ICA

[F, delta] = aci(x);
z = pinv(F) * x;
[sir, ~, D] = sir(z', s');

%% PCA
[P, V] = pca(x');



