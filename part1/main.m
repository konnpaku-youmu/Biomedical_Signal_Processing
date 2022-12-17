close all;
clc;
clear;

load eegdata_artifacts.mat

%% eye-blink artifact removal
% eyeblink_mask = mwf_getmask(eegdata, fs);
load eyeblink_arti.mat

% [n, d, W, SER, ARR, p] = mwf_process(eegdata, eyeblink_mask, 2);
% 
% eegplot_simple(n, fs);
% title(sprintf("SER = %1.2f", SER))
% 
% figure;
% plot(eegdata(3, :), 'DisplayName', "Raw");
% hold on;
% plot(n(3, :), 'DisplayName', "Clean");
% plot(d(3, :), 'DisplayName', "Artifact");
% legend;

%% Muscle artifact removal
% muscle_mask = mwf_getmask(eegdata, fs);
load muscle_arti.mat

% [n, d, W, SER, ARR, p] = mwf_process(eegdata, muscle_mask, 2);
% 
% eegplot_simple(n, fs);
% 
% figure;
% plot(eegdata(3, :), 'DisplayName', "Raw");
% hold on;
% plot(n(3, :), 'DisplayName', "Clean");
% plot(d(3, :), 'DisplayName', "Artifact");
% legend;

%% Muscle and eye-blink
mask_union = eyeblink_mask + muscle_mask;

[n, d, W, SER, ARR, p] = mwf_process(eegdata, mask_union, 2);

eegplot_simple(eegdata, fs);
figure;
eegplot_simple(n, fs);
title(sprintf("SER = %1.2f", SER))

figure;
plot(eegdata(3, :), 'DisplayName', "Raw");
hold on;
plot(n(3, :), 'DisplayName', "Clean");
plot(d(3, :), 'DisplayName', "Artifact");
legend;
