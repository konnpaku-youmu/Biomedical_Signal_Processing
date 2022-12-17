% clc;
clear;

%% Harmonic retrieval

load ex4data.mat

fs = 75;

H_x = hankelize(x, 'order', 3);
[U, S, V] = mlsvd(H_x);

plot(V{1}, 'LineStyle', 'none', 'Marker','x');

[Utrunc,Strunc]=mlsvd(H_x,[7 7 7],0);

% pole estimates
U1toptrunc = Utrunc{1}(1:100,:); 
U1bottomtrunc = Utrunc{1}(2:101,:);

disp('ztrunc:'), disp(eig(pinv(U1toptrunc)*U1bottomtrunc))   % pole estimates truncation ~ z

figure;
plot(x)

%% Prediction
options.Initialization = @mlsvd;
[U_lmlra, S_lmlra] = lmlra(H_x, [7 7 7], options);
