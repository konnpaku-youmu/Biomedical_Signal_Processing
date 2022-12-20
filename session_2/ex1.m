clear;
clc;

load ex1data.mat;

%% 2.1
T1_rank = mlrankest(T1);
[U, S, sv] = mlsvd(T1);

semilogy(sv{1}, LineStyle="--", Marker="*", LineWidth=1);
hold on
semilogy(sv{2}, LineStyle="--", Marker="diamond", LineWidth=1);
semilogy(sv{3}, LineStyle="--", Marker="hexagram", LineWidth=1);

%% 2.2
A = 1/(sqrt(5)) * [1, -2; -2, 1];
sigma = 0.05;

SNRs = 0:5:50;
N = 200;

SIR_ICA = [];
SIR_PCA = [];

%%
for SNR = SNRs
    sir_ica = 0;
    sir_pca = 0;
    for k = 1:N
        s = rand(2, 800) * 2 - 1;
        [x, n] = noisy(A*s, SNR);
        
        % ICA
        [F, delta] = aci(x);
        z = pinv(F) * x;
        [si, ~, ~] = sir(z', s');
        sir_ica = sir_ica + si;
        % PCA
        [P, V] = pca(x');
        z_pca = pinv(P) * x;
        [si, ~, D] = sir(z_pca', s');
        sir_pca = sir_pca + si;
    end
    
    SIR_ICA = [SIR_ICA sir_ica/N];
    SIR_PCA = [SIR_PCA sir_pca/N];
end

figure;
hold on;
plot(SNRs, SIR_ICA, LineWidth=1, DisplayName="ICA");
plot(SNRs, SIR_PCA, LineWidth=1, DisplayName="PCA");
ylabel("SIR");
xlabel("SNR");
title("Source separation quality: ICA v.s. PCA");
legend


