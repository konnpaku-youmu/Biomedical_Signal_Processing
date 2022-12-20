clc;
clear;
close all;

load artifactData.mat

events_start = events(1:2:end);
events_end = events(2:2:end);

n_ch = size(data, 1);
idx_lb = 1;
idx_ub = 100;
t = idx_lb/fs:1/fs:idx_ub/fs;

%% plot the signal
hold on;
eeg_plt = plot(t, data(:, idx_lb:idx_ub)');
e_plt = plot(events(events >= idx_lb & events <= idx_ub) / fs, 0, 'LineStyle','none', 'Marker','*','LineWidth',1, 'Color', 'g', 'DisplayName', 'Events');
l_plt = plot(labels(labels >= idx_lb & labels <= idx_ub) / fs, 0, 'LineStyle','none', 'Marker','>','LineWidth',1, 'Color', 'k', 'DisplayName', 'Labels');
xlabel("t");
legend([eeg_plt(1); e_plt(1); l_plt(1)], {'EEG signal', 'Events', 'Labels'});
title(sprintf("Sample segment of EEG: $t=%1.2f - %1.2f$s", t(1), t(end)));

%% Solve for the filter
L = 25;

for ch = 1:n_ch
    % find the non-adjacent channels
    nadj_ch = setdiff([1:n_ch], adjacent_channels{ch});
    data_nadj_ch = data(nadj_ch, :);
    N_k = size(data_nadj_ch, 1);
    
    % construct X
    X_nk = zeros(length(events_end), N_k*L);
    idx_event = 1;
    for event = events_end
        x_nk = reshape(data_nadj_ch(:, event-L+1:event).', N_k*L, 1);
        X_nk(idx_event, :) = x_nk;
        idx_event = idx_event + 1;
    end
    
    x_k = data(ch, events_end)';
    
    % solve LS
    w_k = X_nk \ x_k;
    
    
    
end
