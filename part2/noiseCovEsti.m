function [cov] = noiseCovEsti(noiseSegment, data)
    cov = 0;
    noise_data = data(noiseSegment);
    lag = 25;
    
    for i = 1:size(noise_data, 1)-lag
        data_esti = noise_data(i:i+lag-1, :);
        cov = cov + data_esti * data_esti';
    end

end