function [Rnn] = noiseCovEsti(noiseSegments, data)
    Rnn=zeros(size(data,2));
    noise_seg_data = data(noiseSegments, :);
    lag = 25;
    
    for i = 1:size(noise_seg_data, 1)-lag
        data_esti = noise_seg_data(i:i+lag-1, :);
        Rnn = Rnn + data_esti' * data_esti;
    end
    Rnn = Rnn./(size(noise_seg_data, 1)-lag);

%     lag = 25;
%     cov=zeros(size(data,2));
%     for i=1:size(noiseSegments,1)-lag
%         segm = data(noiseSegments(i):noiseSegments(i+lag),:);
%         cov = cov + segm'*segm;
%     end
%     
%     Rnn = cov./(size(noiseSegments,1)-lag);
end