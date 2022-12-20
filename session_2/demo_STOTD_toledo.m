
T = 200; R= 2; I = 4;
X = sign(randn(R,T)) + i * sign(randn(R,T));
%figure(1), plot(X(1,:),'x'), axis([-2 2 -2 2])

%%
M = randn(I,R) + i * randn(I,R);
cond(M)
Y = M * X;
[UY,SY,VY] = svd(Y,'econ');
figure(2), plot(Y(1,:),'x')

%%
Me = ica_sttdc1_toledo(Y,2);
Xe = pinv(Me) * Y;
%figure(3), plot(Xe(1,:),'x')

%%

N = randn(I,T) + i * randn(I,T);
Yn = Y / norm(Y,'fro') + 0.1 * N / norm(N,'fro'); % 0.1 = 20 dB
[UYn,SYn,VYn] = svd(Yn);
disp(diag(SYn))
figure(4), plot(Yn(1,:),'x'), title('data 1')

%%
figure(5), plot(VYn(:,1),'x'), title('PC 1')

%%
Men = ica_sttdc1_toledo(Yn,2);
Xen = pinv(Men) * Yn;
figure(6), plot(Xen(1,:),'x'), title('IC 1')
