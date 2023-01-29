clear all; close all; clc;

t = 1:1:10;

tmax_acc = 2*max(t);
Err = 1e-6;
Nk = 40;

numt = numel(t);

N0 = 20;
lambda = 1/5;
lap_sol = @(s) N0./(s+lambda);

num_sol = zeros(numt,1);
for i = 1:1:numt
    tic
    num_sol(i) = invlaplace_dehoog2(lap_sol,t(i),tmax_acc,Err,Nk);
    toc
end

exact_sol = N0*exp(-lambda*t);

plot(t,exact_sol,'r')
hold on
plot(t,num_sol,'ob')
hold off

%%

%Test on an example solution from Bodvarsson - temperature in rock at a
%point over time

t = 0:0.01:1;

tmax_acc = 2*max(t);
Err = 1e-6;
Nk = 40;

numt = numel(t);

th = 100;
xi = 0.1;
eta = 0.8;

clear lap_rock;
lap_rock_bodvar = @(p) 1/p*exp(-((th*p + 2*p^.5*tanh(p^.5))*xi/(2+th)))...
    *(cosh(p^.5*eta) - sinh(p^.5*eta)*tanh(p^.5));

num_sol_rock = zeros(numt,1);
for i = 1:1:numt
    tic
    num_sol_rock(i) = invlaplace_dehoog2(lap_rock_bodvar,t(i),tmax_acc,Err,Nk);
    toc
end

plot(t,num_sol_rock)

%%


%Test on an example solution from Gringarten - temperature in rock at a
%point over time

t = logspace(-1,3,30);

tmax_acc = 2*max(t);
Err = 1e-6;
Nk = 40;

numt = numel(t);

z_D = 1;
alpha = 2;
figure(1); clf

for x_ED = [0.5 1 2 4 8 1000]

    clear lap_rock_grin;
    lap_rock_grin = @(s) 1/s*exp(-z_D*s^.5*tanh((x_ED)/alpha*s^.5));
    
    num_sol_rock = zeros(numt,1);
    for i = 1:1:numt
        tic
        tmax_acc = 2*t(i);
        num_sol_rock(i) = invlaplace_dehoog2(lap_rock_grin,t(i),tmax_acc,Err,Nk);
        toc
    end

    hold on
    pl1 = plot(t,num_sol_rock);
    hold off
    set(gca,'XScale','log')
    pause
end
