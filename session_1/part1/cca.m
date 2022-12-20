clc;
clear;
close all;

load eegdata_artifacts.mat

%% Plot original EEG data
eegplot_simple(eegdata, fs);

%% BSS using CCA
eegdata_d = delayseq(eegdata', 1)';
[Wx, Wy, r, x, y] = canoncorr(eegdata', eegdata_d');

figure;
plot(r, 'LineWidth', 1);

n_srcs = 30;
w_inv = inv(Wx);
eeg_clear = x(:, 1:n_srcs) * w_inv(1:n_srcs, :);

figure;
eegplot_simple(eeg_clear', fs);
