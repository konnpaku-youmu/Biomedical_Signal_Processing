%% data

t = linspace(0,50,200);
V{1} = [sin(2*pi*t/50)', sin(2*pi*t/25)'];
V{2} = V{1};
V{3} = V{1}(1:5:end,:);
Y = cpdgen(V);

Yinc = Y;
completionfrac = 0.8;           % completionfrac*100 = percentage missing
completionperc = completionfrac*100;
for i=1:200, for j=1:200, for k=1:40
          %a = rand(1); 
          if rand(1) < completionfrac, Yinc(i,j,k) = NaN; end
end, end, end

figure(1), surf3(Y)             % "full data"
title('true tensor')

figure(2), surf3(Yinc)          % "incomplete data"
title(['incomplete tensor: ', num2str(completionperc), '% missing'])


%% completion by low ML rank approximation

[UYinc,SYinc] = lmlra(Yinc,[2 2 2]); 
Ylmlra = lmlragen(UYinc,SYinc);

figure(3), surf3(Ylmlra)
title('completed tensor')
disp('relative error completed vs true tensor:')
frob(Ylmlra-Y)/frob(Y)

