clc;
clear;
close all;

load eegdata_artifacts.mat

figure;
eegplot_simple(eegdata);

eegdata = eegdata';
eegdata_delay = delayseq(eegdata, 1);

[A, B, r, U, V] = canoncorr(eegdata_delay, eegdata);

U(:, 39:end) = 0;
eegdata_clean = U * inv(A);

%%
figure;

t = 11000/fs:1/fs:13000/fs;
plot(t, eegdata(11000:13000, 41), DisplayName="Raw");
hold on
plot(t, eegdata_clean(11000:13000, 41), DisplayName="Clean");
xlabel("Time(s)");
legend;

title(sprintf("Muscle artifacts removal: Channel %s", channelnames{41}));

%%
figure;
hold on;
plot(r);
plot(39, r(39), Marker="*", MarkerSize=8, LineWidth=1);
ylabel("Correlation");
xlabel("N-th source");
title("Sample canonical correlation");

