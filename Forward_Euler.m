%%
% This is the MATLAB Code for forward euler algorithm by CHEN Jiawei
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
fig1 = plot(x,y,'b'); hold on
axis([-2.5 0 -2.5 0.5]);
xlabel('Displacement');
ylabel('Force');
title('Forward Euler method with k_b_a_r = 1.0');
%% set the para for iter
%
F_ext=-2;
F_ini=0;
u = 0;
fig1 = plot(x,F_ext*ones(1,length(x)),'b--'); hold on
%
delta_f = -0.05;
num_step=F_ext./delta_f;

%
for step = 1:1:(num_step+1)
    if step==1
        u(step)=0;
        F_ini(step)=0;
        F_ext(step)=F_ini(step)+delta_f;
        k = k_stif(u(step),k_bar);
        delta_u = delta_f./k;
        u(step+1)=u(step)+delta_u;
        F_ini(step+1)=F_ini(step)+delta_f;
    else
        F_ext(step)=F_ini(step)+delta_f;
        k = k_stif(u(step),k_bar);
        delta_u = delta_f./k;
        u(step+1)=u(step)+delta_u;
        F_ini(step+1)=F_ini(step)+delta_f;
        plot(u(step+1),F_ini(step+1),'ro','Markersize',5,'Markerfacecolor','r'); hold on
    end
end

fig2 = plot(u,F_ini,'ro','Markersize',5,'Markerfacecolor','r'); hold on
legend([fig1 fig2],'Analytical solution','Forward euler alogorithm','location','best');


%%
%function definition
function y=F_inter(x,k_bar)
y=x + 1.5*x.^2 + 0.5*x.^3 + k_bar * x;
end
function y=k_stif(x,k_bar)
y = (1 + 3 * x + 1.5 * x.^2) + k_bar;
end