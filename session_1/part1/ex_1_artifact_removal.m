close all;
clc;
clear;

%% 1.2 CCA
load eegdata_artifacts.mat

figure;
eegplot_simple(eegdata);

eegdata = eegdata';
eegdata_delay = delayseq(eegdata, 1);

[A, B, r, U, V] = canoncorr(eegdata_delay, eegdata);

U(:, 39:end) = 0;
eegdata_clean = U * inv(A);

figure;

t = 11000/fs:1/fs:13000/fs;
plot(t, eegdata(11000:13000, 41), DisplayName="Raw");
hold on
plot(t, eegdata_clean(11000:13000, 41), DisplayName="Clean");
xlabel("Time(s)");
legend;

title(sprintf("Muscle artifacts removal: Channel %s", channelnames{41}));

figure;
hold on;
plot(r);
plot(39, r(39), Marker="*", MarkerSize=8, LineWidth=1);
ylabel("Correlation");
xlabel("N-th source");
title("Sample canonical correlation");


%% 1.3 MWF: eye-blink artifact removal
% eyeblink_mask = mwf_getmask(eegdata, fs);
% eyeblink_mask(isnan(eyeblink_mask)) = 0;

load eegdata_artifacts.mat
load blink_artifacts.mat
load muscle_arti.mat

delay = 3;
[n, d, W, SER, ARR, p] = mwf_process(eegdata, eyeblink_mask, delay);
single_ch_plot(eegdata, n, d, fs, 15, 2, delay, SER, ARR);

% Muscle artifact removal
% muscle_mask = mwf_getmask(eegdata, fs);
% muscle_mask(isnan(muscle_mask)) = 0;

% CCA
eegdata_delay = delayseq(eegdata', 1);
[A, B, r, U, V] = canoncorr(eegdata_delay, eegdata');

U1 = U;
U2 = U;
U1(:, 39:end) = 0;
U2(:, 1:38) = 0;
eeg_cca = U1 * inv(A);
eeg_cca_arti = U2 * inv(A);

% MWF
[n, d, W, SER, ARR, p] = mwf_process(eegdata, muscle_mask, 3);
[SER_cca, ARR_cca] = mwf_performance(eeg_cca', eeg_cca_arti', muscle_mask);

figure;
range = 11000:13000;
t = range / fs;
plot(t, eegdata(41, range), 'DisplayName', "Raw");
hold on;
plot(t, d(41, range), 'DisplayName', "Artifact");
plot(t, n(41, range), 'DisplayName', "Clean");
legend;
xlabel("Time (s)");
title(sprintf("MWF artifact removal: Delay = %d, SER=%1.2f, ARR=%1.2f", 3, SER, ARR));

figure;
plot(t, d(41, range), 'DisplayName', "Artifact", "Color", [0.7 0.7 0.7]);
hold on;
plot(t, eeg_cca(range, 41), 'DisplayName', "CCA");
plot(t, n(41, range), 'DisplayName', "MWF");
legend;
xlabel("Time (s)");
title(sprintf("CCA v.s. MWF: $SER_{CCA}$ = %1.2f, $ARR_{CCA}$ = %1.2f, $SER_{MWF}$ = %1.2f, $ARR_{MWF}$ = %1.2f", SER_cca, ARR_cca, SER, ARR));

% Muscle and eye-blink
mask_union = eyeblink_mask + muscle_mask;

[n, d, W, SER, ARR, p] = mwf_process(eegdata, mask_union, 2);

figure;

range = 1000:2000;
t = range / fs;
subplot(121);
plot(t, eegdata(3, range), 'DisplayName', "Raw");
hold on;
plot(t, n(3, range), 'DisplayName', "Clean");
plot(t, d(3, range), 'DisplayName', "Artifact");
legend;
xlabel("Time (s)");
title("Eye-blink segment");

range = 9000:11000;
t = range / fs;
subplot(122);
plot(t, eegdata(3, range), 'DisplayName', "Raw");
hold on;
plot(t, n(3, range), 'DisplayName', "Clean");
plot(t, d(3, range), 'DisplayName', "Artifact");
legend;
xlabel("Time (s)");
title("Muscle segment");

sgtitle(sprintf("MWF artifact removal on Channel 3: SER=%1.2f, ARR=%1.2f", SER, ARR));


function single_ch_plot(sig_org, sig_filt, arti_esti, fs, length, i_ch, d, SER, ARR)
    L = fs * length;
    figure("Name", "MWF artifact removal");
    hold on;
    plot(linspace(0, length, L), sig_org(i_ch, 1:L), 'DisplayName', 'Raw signal', 'LineWidth', 0.75);
    plot(linspace(0, length, L), arti_esti(i_ch, 1:L), 'DisplayName', 'Artifacts', 'LineWidth', 0.75);
    plot(linspace(0, length, L), sig_filt(i_ch, 1:L), 'DisplayName', 'Filtered', 'LineWidth', 0.75);
    legend;
    title(sprintf("MWF artifact removal: Delay = %d, SER=%1.2f, ARR=%1.2f", d, SER, ARR));
end
