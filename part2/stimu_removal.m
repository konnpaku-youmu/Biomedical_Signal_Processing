clc;
clear;
close all;

load artifactData.mat

n_ch = size(data, 1);
idx_lb = 108900;
idx_ub = 109100;

% plot the signal
hold on;
eeg_plt = plot(data(:, idx_lb:idx_ub)');
events_plt = plot(events(events >= idx_lb & events <= idx_ub) - idx_lb, 0, 'LineStyle','none', 'Marker','*','LineWidth',2, 'Color', 'g', 'DisplayName', 'Events');
labels_plt = plot(labels(labels >= idx_lb & labels <= idx_ub) - idx_lb, 0, 'LineStyle','none', 'Marker','>','LineWidth',2, 'Color', 'r', 'DisplayName', 'Labels');
legend([eeg_plt; events_plt; labels_plt], {'EEG signal', 'Events', 'Labels'});

for ch = 1:n_ch
    % find the non-adjacent channels
    nadj_ch = setdiff([1:n_ch], adjacent_channels{ch});
    
    

end