function [] = plotMultiChannel(t, data, spacing, varargin )
%PLOTMULTICHANNEL Plot multi-channel data in a single figure with given
%spacing between individual channels.
%
% inputs:
%   data : multi-channel data (numberSamples, numberChannels)
%   spacing : a number which is multiplied with the mean of the individual
%               channel standard deviations (e.g., 2)
%   varargins : variable number of inputs that are directly passed to the
%               MATLAB plot()
%

% derive parameters from given data
nbChannels = size(data, 2);
stdData = mean(std(data));

% construct spacing vector
spacingData = 0:(nbChannels-1);
spacingData = spacingData * stdData * spacing;

% apply spacing
spacedData = bsxfun(@plus, data, spacingData);

if isempty(t)
    t = (1:size(data, 1))';
end

plot(t, spacedData, varargin{:});

if isempty(t)
    xlabel('discrete time [samples]');
else
    xlabel('Time [s]');
end

ylabel('amplitude [arb. unit]')

end

