function [Rnn] = noiseCovEsti(noiseSegments, data)
    n_ch = size(data, 2);
    noise_seg_data = data(noiseSegments, :);
    lag = 25;
    Rnn=zeros(n_ch*lag);
%     Rnn=zeros(n_ch);
    
    for i = 1:10:size(noise_seg_data, 1)-lag
        data_esti = noise_seg_data(i:i+lag-1, :);
        data_esti = mat2stacked(data_esti)';
        Rnn = Rnn + (data_esti' * data_esti);
    end
    Rnn = Rnn./(size(noise_seg_data, 1)-lag);
end