clc;
clear;
close all;

load eegdata_artifacts.mat

eegdata = eegdata';
eegdata_delay = delayseq(eegdata, 1);

[A, B, r, U, V] = canoncorr(eegdata_delay, eegdata);

figure;
eegplot_simple(U');

U(:, 39:end) = 0;
eegdata_clean = U * inv(A);

% figure;
% eegplot_simple(eegdata');
figure;
eegplot_simple(eegdata_clean');


