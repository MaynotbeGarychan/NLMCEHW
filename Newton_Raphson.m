%%
% This is the MATLAB Code for Newton Raphson algorithm by CHEN Jiawei
%%
clear
clc
close all
format longEng
%% plot the analytical solution
k_bar = 1.0;
x=0:-0.01:-2.5;
y=x + 1.5*x.^2 + 0.5*x.^3 + k_bar * x;
figure(1)
set(gca,'Fontsize',16)
fig1 = plot(x,y,'b'); hold on
axis([-2.5 0 -2.5 0.5]);
xlabel('Displacement');
ylabel('Force');
title('Newton Raphson method with k_b_a_r = 1.0');
%% Increment
% set param for increment
F_ext=-2;
F_ini=0;
fig1 = plot(x,F_ext*ones(1,length(x)),'b--'); hold on
%
delta_f = -0.05;
num_step = F_ext/-0.02;
tol=1e-6;
max_iter = 100;
% init increment
for step = 1:1:(num_step+1)
    if step==1
        u(step)=0;
        F_ini(step)=0;
        F_ext(step)=F_ini(step)+delta_f;
    else
        u(step)=u(step-1)+delta_u;
        F_ini(step)=F_ini_j;
        F_ext(step)=F_ini(step)+delta_f;
    end
    %init iter
    iter=1;
    delta_u = 0;
    Residual = delta_f;
    k = k_stif(u(step),k_bar);
    %k = 2;
    while(iter<max_iter && abs(Residual)>=tol)
        ddelta_u=Residual/k;
        delta_u = delta_u + ddelta_u;
        F_ini_j=F_inter(u(step)+delta_u,k_bar);
        Residual = F_ext(step) - F_ini_j;
        iter = iter +1;
    end
    fig3 = plot(u(step),F_ini(step),'go','Markersize',5,'Markerfacecolor','g'); hold on
    u(step+1) = u(step)+delta_u;
    
    
end


legend([fig1 fig3],'Analytical solution','Newton Raphson alogorithm','location','best');

%%
%function definition
function y=F_inter(x,k_bar)
y=x + 1.5*x.^2 + 0.5*x.^3 + k_bar * x;
end
function y=k_stif(x,k_bar)
y = (1 + 3 * x + 1.5 * x.^2) + k_bar;
end
    