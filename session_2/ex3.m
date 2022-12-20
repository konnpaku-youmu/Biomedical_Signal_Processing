clc;
clear;

load ex3data.mat

%% 2.3
% Plotting
figure,
subplot(131);
plot(A);
title('Components of A');
legend('Comp.1','Comp.2','Comp.3','Comp.4');

subplot(132);
plot(B);
title('Components of B');
legend('Comp.1','Comp.2','Comp.3','Comp.4');

subplot(133);
plot(C);
title('Components of C');
legend('Comp.1','Comp.2','Comp.3','Comp.4');

%% 2.3.1
%PCA
T_slice = squeeze(T3(1, :, :));
T_sliceA = squeeze(T3(:, :, 1));
[P, V] = pca(T_slice');
[PA, VA] = pca(T_sliceA');

figure;
proj_A = T_sliceA * VA(:, 1:4);
subplot(131);
hold on
plot(proj_A);
title('Estimation of A');
legend('Comp.1','Comp.2','Comp.3','Comp.4');

proj_B = T_slice * V(:, 1:4);
subplot(132);
plot(proj_B);
title('Estimation of B');
legend('Comp.1','Comp.2','Comp.3','Comp.4');

proj_C = (P(:, 1:4)' * T_slice)';
subplot(133);
plot(proj_C);
title('Estimation of C');
legend('Comp.1','Comp.2','Comp.3','Comp.4');

% sgtitle("Estimation of components in $T_3$ using PCA");

%% 2.3.2
T_noisy = noisy(T3, 25);
[Uml, Sml, sv] = mlsvd(T_noisy);

figure;
semilogy(sv{1}, DisplayName="$\sigma^{(1)}$", LineStyle="none", Marker="*", LineWidth=1);
hold on
semilogy(sv{2}, DisplayName="$\sigma^{(2)}$", LineStyle="none", Marker="x", LineWidth=1);
semilogy(sv{3}, DisplayName="$\sigma^{(3)}$", LineStyle="none", Marker="+", LineWidth=1);
title("Multilinear singular values of tensor $T_{noisy}$");
legend;

% rank = 4
figure;

for rank=1:4
    options = struct;
    options.Algorithm = @cpd_nls;
    options.ExploitStructure = true;
    [U_NLS,output_NLS] = cpd(T_noisy,rank,options);

    options.Algorithm = @cpd_als;
    options.ExploitStructure = true;
    [U_ALS,output_ALS] = cpd(T_noisy,rank,options);
    
    subplot(2, 2, rank);
    hold on;
    plot(output_NLS.Algorithm.fval);
    plot(output_ALS.Algorithm.fval);
    hold off;
    
    ylabel('Loss');
    xlabel('iter');
    plt_title = sprintf('rank = %d',rank);
    title(plt_title);
    legend('CPD: NLS','CPD: ALS');
    set(gca, 'yscale', 'log');

end

sgtitle("Convergence: NLS v.s. ALS");


