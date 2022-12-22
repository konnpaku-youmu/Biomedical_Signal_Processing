clc;
clear;
close all;

load eegdata_artifacts.mat
load blink_artifacts.mat

% % mask_blink = mwf_getmask(eegdata, fs);
% delays = 0:8;
% 
% SERs = [];
% ARRs = [];
% for delay = delays
%     [n, d, W, SER, ARR, p] = mwf_process(eegdata, mask_blink, delay);
%     SERs = [SERs SER];
%     ARRs = [ARRs ARR];
% end
% 
% figure;
% hold on;
% plot(delays, SERs, DisplayName="SER");
% ylabel("SER (dB)");
% yyaxis right;
% plot(delays, ARRs, DisplayName="ARR");
% ylabel("ARR (dB)");
% xlabel("Delay (taps)");
% legend;

delay = 3;
[n, d, W, SER, ARR, p] = mwf_process(eegdata, eyeblink_mask, delay);
single_ch_plot(eegdata, n, d, fs, 15, 2, delay, SER, ARR);

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
