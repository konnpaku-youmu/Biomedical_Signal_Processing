clc;
clear;

load ex4data.mat

fs = 75;
t = linspace(0, length(x) / 75, length(x));
t_pred = linspace(0, length(x) / 75 + 1, length(x)+75);

H_x = hankelize(x, 'order', 3);

%% Harmonic retrieval: CPD
U_cpd = cpd_rnd(H_x, 6);
[Uhat,output] = cpd(H_x,U_cpd);

Utop = Uhat{1}(1:100,:);
Ubottom = Uhat{1}(2:101,:);

z_cpd = eig(Utop \ Ubottom);

disp(z_cpd);

%% Harmonic retrieval: SVD
[Utrunc,Strunc]=mlsvd(H_x,[6 6 6],0);

% pole estimates
U1toptrunc = Utrunc{1}(1:100,:); 
U1bottomtrunc = Utrunc{1}(2:101,:);

z_svd = eig(pinv(U1toptrunc)*U1bottomtrunc);   % pole estimates truncation ~ z


%% pzmap

figure;

[U, S, V] = mlsvd(H_x);

figure
semilogy(V{1}, 'LineStyle', 'none', 'Marker','x', DisplayName="$\sigma^{(1)}$", LineWidth=1);
title("Singular values of $\mathcal{H}_x$");
legend

figure
hold on;
plot(z_cpd, LineStyle="none", Marker="x", LineWidth=1);
plot(z_svd, LineStyle="none", Marker="diamond", LineWidth=1);
xline(0, LineStyle=":", Color='b');
yline(0, LineStyle=":", Color='b');
plot(cos(linspace(0, 2*pi, 1000)), sin(linspace(0, 2*pi, 1000)), LineStyle=":", Color='b');

axis equal
xlim([-1.5, 1.5]);
ylim([-1.5, 1.5]);
xlabel("$Re(z)$");
ylabel("$Im(z)$");
legend("Poles: CPD", "Poles: ESPRIT");
title("Pole-Zero map");

%% Prediction: CPD, SVD, Completion

% Prediction with poles: CPD
z = z_cpd;
pred = (z.^[0:300]).';
coeff = pred\x.';

pred = [pred; (z.^[301:375]).'];

xest = pred * coeff;

figure;

hold on;
plot(t, x, DisplayName="Given signal");
plot(t_pred, xest,'r', DisplayName="Prediction");
title("Harmonic retrieval by CPD: predicting with poles");
xlabel("t");
legend;


% Prediction with poles: SVD
z = z_svd;
pred = (z.^[0:300]).';
coeff = pred\x.';

pred = [pred; (z.^[301:375]).'];

xest = pred * coeff;

figure;
hold on;
plot(t, x, DisplayName="Given signal");
plot(t_pred, xest,'r', DisplayName="Prediction");
title("Harmonic retrieval by ESPRIT: predicting with poles");
xlabel("t");
legend;

% completion with LMLRA
H_x_incomp = hankelize([x zeros(1, fs)], 'order', 3);
[UYinc,SYinc] = lmlra(H_x_incomp,[6 6 6]); 
Hlmlra = lmlragen(UYinc,SYinc);

x_pred = dehankelize(Hlmlra, 'order', 3);

figure;
hold on;
plot(t, x, DisplayName="Given signal");
plot(t_pred, x_pred,'r', DisplayName="Prediction");
title("Predicting using tensor completion");
xlabel("t");
legend;



