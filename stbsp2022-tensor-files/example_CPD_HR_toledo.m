

close all; clc;

% Hx1 from example Hankelization
t = 0:0.01:4.98;
s1 = 2*exp(-t);
s2 = exp(-0.3*t).*cos(3*t);
x = s2-s1;
%xn = noisy(x,10);
z = [exp(-0.01), exp(-0.003+1i*0.03), exp(-0.003-1i*0.03)];
disp(z); % poles
figure(1); plot(t,x)

%%
Hx1 = hankelize(x,'order',3); 
figure(2); surf3(Hx1)

%%
% noisy version

xn = noisy(x,10);
Hx1n = hankelize(xn,'order',3); 
figure(1), hold, plot(t,xn,'xr')
figure(4), surf3(Hx1n)

%%
% Hx1n from example low ML rank approximation
% Hx1n = Hx1 + 0.3*randn(167,167,167); % add noise

U0{1} = randn(167,3)+ 1i*randn(167,3);
U0{2} = randn(167,3)+ 1i*randn(167,3);
U0{3} = randn(167,3)+ 1i*randn(167,3);
[Uhat,output] = cpd(Hx1n,U0);
figure(5)
semilogy(0:output.Algorithm.iterations,sqrt(2*output.Algorithm.fval));
xlabel('iteration');
ylabel('frob(T-cpdgen(U))');
% ideally, the columns of Uhat{1} are scaled Vandermonde vectors
grid on;
U1top = Uhat{1}(1:166,:); U1bottom = Uhat{1}(2:167,:);
% pole estimates zrest; solve U1top(:,r) zrest = U1bottom(:,r)
disp(pinv(U1top(:,1))*U1bottom(:,1)),
disp(pinv(U1top(:,2))*U1bottom(:,2))
disp(pinv(U1top(:,3))*U1bottom(:,3))