clc;
clear;
close all;

load eegdata_artifacts.mat
load blink_artifacts.mat

% mask_blink = mwf_getmask(eegdata, fs);
[n, d, W, SER, ARR, p] = mwf_process(eegdata, mask_blink, 0);
single_ch_plot(eegdata, n, d, fs, 60, 1);


function single_ch_plot(sig_org, sig_filt, arti_esti, fs, length, i_ch)
    L = fs * length;
    figure("Name", "MWF artifact removal");
    hold on;
    plot(linspace(0, length, L), sig_org(i_ch, 1:L), 'DisplayName', 'Raw signal', 'Color', '#837EE6', 'LineWidth', 0.75);
    plot(linspace(0, length, L), arti_esti(i_ch, 1:L), 'DisplayName', 'Artifacts', 'Color', '#E68981', 'LineWidth', 0.75);
    plot(linspace(0, length, L), sig_filt(i_ch, 1:L), 'DisplayName', 'Filtered', 'Color', '#7EE6C1', 'LineWidth', 0.75);
    legend;
    title("MWF artifact removal");
end
