clc;
clear;
close all;

load artifactData.mat

events_start = events(3:2:end);
events_end = events(4:2:end);

n_ch = size(data, 1);
idx_lb = 73800;
idx_ub = 74100;
range = idx_lb:idx_ub;
t = range ./ fs;

% %% plot the signal
% eeg_plt = plot(t, data(:, idx_lb:idx_ub)');
% hold on;
% e_plt = plot(events(events >= idx_lb & events <= idx_ub) / fs, 0, 'LineStyle','none', 'Marker','*','LineWidth',1, 'Color', 'g', 'DisplayName', 'Events');
% l_plt = plot(labels(labels >= idx_lb & labels <= idx_ub) / fs, 0, 'LineStyle','none', 'Marker','>','LineWidth',1, 'Color', 'k', 'DisplayName', 'Labels');
% xlabel("t");
% legend([eeg_plt(1); e_plt(1); l_plt(1)], {'EEG signal', 'Events', 'Labels'});
% title(sprintf("Sample segment of EEG: $t=%1.2f - %1.2f$s", t(1), t(end)));

%% Solve for the filter
mask = zeros(1, size(data, 2));
for event_int = [events_start;events_end]
    mask(event_int(1):event_int(end)) = 1;
end

Ls = 1:10;
SERs = zeros(10, 1);
ARRs = zeros(10, 1);

for L = Ls

    t_M = [];
    for event_int = [events_start;events_end]
        t_M = [t_M event_int(1):event_int(end)];
    end
    
    y = zeros(size(data));
    arti = zeros(size(data));
    
    for ch = 1:n_ch
        % find the non-adjacent channels
        nadj_ch = setdiff([1:n_ch], adjacent_channels{ch});
        data_nadj_ch = data(nadj_ch, :);
        N_k = size(data_nadj_ch, 1);
        
        % construct X
        X_nk = zeros(length(t_M), N_k*L);
        
        k_t = 1;
        for t_m = t_M
            x_nk = reshape(data_nadj_ch(:, t_m-L+1:t_m).', N_k*L, 1);
            X_nk(k_t, :) = x_nk;
            k_t = k_t + 1;
        end
        
        x_k = data(ch, t_M)';
        
        % solve LS
        w_k = X_nk \ x_k;
        
        % apply filter to the data
        y_k = data(ch, :);
        arti_k = zeros(size(y_k));
    
        for event_int = [events_start;events_end]
            t_M = event_int(1):event_int(end);
            for t_m = t_M
                x_nk = reshape(data_nadj_ch(:, t_m-L+1:t_m).', N_k*L, 1);
                arti_k(t_m) = w_k'*x_nk;
                y_k(t_m) = y_k(t_m) - arti_k(t_m);
            end
        end

        y(ch, :) = y_k;
        arti(ch, :) = arti_k; 
    end

    [SERs(L), ARRs(L)] = mwf_performance(data(:, events_start(1):end), arti(:, events_start(1):end), mask(range));
end

%%
labels_eeg = labels(labels >= idx_lb & labels <= idx_ub) / fs - idx_lb / fs;

figure;
eegplot_simple(data(:, range), fs);
hold on;
label1 = xline([labels_eeg - 2e-4; labels_eeg + 1e-3], LineWidth=2, Color='r', LineStyle='--');
legend([label1(1)], {'Neuron Spikes'});
figure;
label2 = xline([labels_eeg - 2e-4; labels_eeg + 1e-3], LineWidth=2, Color='r', LineStyle='--');
hold on;
eegplot_simple(y(:, range), fs);
legend([label2(1)], {'Neuron Spikes'});
