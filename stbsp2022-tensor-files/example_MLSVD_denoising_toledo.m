%% data

t = linspace(0,50,200);
V{1} = [sin(2*pi*t/50)', sin(2*pi*t/25)'];
V{2} = V{1};
V{3} = V{1}(1:5:end,:);
Y = cpdgen(V);
Yn = noisy(Y,10);           % SNR = 10 dB

figure(1), surf3(Y)         % "noisefree data", vary i
title('true tensor')
figure(2), surf3(Yn)        % "noisy data", vary i
title('noisy tensor')

%% inspection ML sing value spectrum --> estimation ML rank

[UYn,SYn,svYn] = mlsvd(Yn);
disp('mode-1 singular values:')
disp(svYn{1}(1:10))         % 2 significant sing values

%% denoising by MLSVD truncation

% truncated MLSVD
[UYntrunc,SYntrunc] = mlsvd(Yn,[2 2 2]);
Yntrunc = lmlragen(UYntrunc,SYntrunc);

figure(3), surf3(Yntrunc)                       % "after noise reduction", vary i
title('denoised tensor')
figure(4), visualize(Y,'original',Yntrunc)      % result vs noisefree data
title('true vs denoised tensor')                % e.g. fiber (:,40,15)
figure(5), visualize(Yntrunc,'original',Yn)     % result vs noisy data
title('noisy vs denoised tensor')  

% accuracy 
disp('relative difference noisy vs denoised tensor [dB]:')
esterrYn = 20*log10(frob(Yn-Yntrunc)/frob(Yn))  % ~ -10dB, in line with noise level
disp('relative difference true vs denoised tensor [dB]:')
esterrY = 20*log10(frob(Y-Yntrunc)/frob(Y))     % << -10dB, good approximation of true tensor 

%% comparison truncation and optimal approximation

[UYnopt,SYnopt] = lmlra(Yn, [2 2 2]);           % optimal approximation
Ynopt = lmlragen(UYnopt,SYnopt);
disp('relative difference true vs optimally denoised tensor [dB]:')
esterrYopt = 20*log10(frob(Y-Ynopt)/frob(Y))    % little difference between opt and trunc


