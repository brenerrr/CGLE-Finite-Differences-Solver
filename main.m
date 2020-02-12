

clear all 
clc
%% Preproc 

% Inputs
N = 512; 
l = 0.05; 
R = 1;
alpha = 0; 
beta = -4; 
tN = 100;
theta = 0.5; % 1 for fully implicit and 0 for fully explicit
x0 = -100; 
xN = 100;

% Initial condition
x = linspace(x0, xN, N+2);
h = x(2) - x(1);
rng default % Reset random numbers seed
% A0 = x*0 + 0.01*rand(1,length(x)) - 0.01/2;
% A0 = x*0 + 0.01*rand(1,length(x)) - 0.01/2 + 1;
% A0 = sqrt(1 - (20*pi/100)^2) * exp(complex(0,20*pi*x/100)) + 0.01*rand(1, length(x)) - 0.01/2;
% A0 =  sech((x+10).^2) + 0.8 * sech((x-30).^2) + 0.01*rand(1,length(x)) - 0.01/2;
A0 =  sech((x+50).^2) + 0.8 * sech((x-50).^2) + 0.01*rand(1,length(x)) - 0.01/2;

% Setup
Nt = round(tN / l);
A = zeros(length(x),Nt);
A(:,1) = A0;
p = l / h^2;

%% Proc
for i = 1:Nt-1 
   b = rhs(A(2:end-1,i), alpha, theta, p, l, R); 
   E = lhs(A(2:end-1,i), alpha, theta, p, l, beta);
   
   A(2:end-1,i+1) = linsolve(E, b);
   A(1, i+1) = 1/6 * ( 4 * A(2,i+1) - A(3,i+1) - A(N,i+1) + 4*A(N+1,i+1) ); 
   A(end, i+1) = A(1, i+1);
   
   if mod(i, round(Nt/10)) == 0 
      disp([num2str(i/Nt*100), '%']) 
   end
end
disp('100%') 

%% Plot
close all 

[T, X] = meshgrid(linspace(0, tN, Nt), x); 
figure(1)
plot(x,real(A(:,end)), '-')
title('Last time step')
xlabel('x') 
ylabel('Re(A)')

figure(2)
plot(x,real(A(:,1)), '-')
title('Initial condition')
xlabel('x') 
ylabel('Re(A)')

figure(3)
h = contourf(X, T, real(A),'LineStyle','none');
colormap(jet)
colorbar();
caxis([-1 1]);
title('Re(A)')
xlabel('x') 
ylabel('t')

figure(5)
h = contourf(X(1:2:end,1:2:end), T(1:2:end,1:2:end), abs(A(1:2:end,1:2:end)),'LineStyle','none');
colormap(jet)
colorbar();
caxis([0 1]);
title('Abs(A)') 
xlabel('x') 
ylabel('t')

%% Animation
% figure(10) 
% Areal = real(A);
% Aabs = abs(A);
% for i = 1:Nt 
%     plot(x,Aabs(:,i))
%     ylim([-1, 1]) 
%     pause(0.01)
% end

