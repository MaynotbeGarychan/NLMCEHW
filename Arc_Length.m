%%
% This is the MATLAB Code for Arc-length algorithm by CHEN Jiawei
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
fig1=plot(x,y,'b'); hold on
axis([-2.5 0 -2.5 0.5]);
xlabel('Displacement');
ylabel('Force');
title(['Arc length method with k_b_a_r =',num2str(k_bar)]);

%% arc length method
% set param for increment
F_ext=2;
F_ini=0;
fig1 = plot(x,-F_ext*ones(1,length(x)),'b--'); hold on
%
tol=1e-3;
max_iter = 100;
step=1;
max_step = 100;
Residual_inc = 1;
%
deltalL = 0.05;
beta = 1.0;

while(step<=max_step)
    if step==1
        u(step)=0;
        lambda(step)=0;
    else
        u(step) = u(step-1) + delta_u(iter);
        lambda(step) = lambda(step-1) + delta_lambda(iter);
        if F_inter(u(step),k_bar)<-2
            break;
        end
    end
    %init iter
    iter=1;
    delta_u=0;
    delta_lambda=0;
    k = k_stif(u(step),k_bar);
    Residual=1;
    
    while(iter<=max_iter && abs(Residual)>=tol)
        if iter==1
            ddelta_u2=0;
            ddelta_u1=F_ext/k;
            
            a1=ddelta_u1*ddelta_u1+beta*F_ext*F_ext;
            a2=2*ddelta_u1*(delta_u(iter)+ddelta_u2)+2*delta_lambda(iter)*(beta^2)*F_ext*F_ext;
            a3=(delta_u(iter)+ddelta_u2)*(delta_u(iter)+ddelta_u2)-deltalL^2+(delta_lambda(iter)^2)*(beta^2)*F_ext*F_ext;
            
            if step==1
                fac=a2*a2-4*a1*a3;
                ddelta_lambda1=(-a2+sqrt(fac))/(2*a1);
                ddelta_lambda2=(-a2+sqrt(fac))/(2*a1);
                ddelta_lambda = max(ddelta_lambda1,ddelta_lambda2);
            else
                ddelta_lambda=lambda_solve(a1,a2,a3,delta_u(iter),ddelta_u,delta_lambda(iter),F_ext);
            end
            
            ddelta_u = ddelta_lambda*ddelta_u1 + ddelta_u2;
            delta_lambda(iter+1)=delta_lambda(iter)+ddelta_lambda;
            delta_u(iter+1)=delta_u(iter)+ddelta_u;
            
            %fig2 = plot(u(step)+delta_u(iter+1),(lambda(step)+delta_lambda(iter+1))*F_ext,'go','Markersize',5,'Markerfacecolor','g'); hold on
            %fig4 = plot(u(step)+delta_u(iter+1),F_inter(u(step)+delta_u(iter+1),k_bar),'ro','Markersize',5,'Markerfacecolor','r'); hold on
            Residual = -F_inter(u(step)+delta_u(iter+1),k_bar)+(lambda(step)+delta_lambda(iter+1))*F_ext;
            
        else
            ddelta_u2=Residual./k;%bar
            ddelta_u1=F_ext/k;%t
            
            a1=ddelta_u1*ddelta_u1+beta*F_ext*F_ext;
            a2=2*ddelta_u1*(delta_u(iter)+ddelta_u2)+2*delta_lambda(iter)*(beta^2)*F_ext*F_ext;
            a3=(delta_u(iter)+ddelta_u2)*(delta_u(iter)+ddelta_u2)-deltalL^2+(delta_lambda(iter)^2)*(beta^2)*F_ext*F_ext;
            
            ddelta_lambda=lambda_solve(a1,a2,a3,delta_u(iter),ddelta_u,delta_lambda(iter),F_ext);
            ddelta_u = ddelta_lambda*ddelta_u1 + ddelta_u2;
            delta_lambda(iter+1)=delta_lambda(iter)+ddelta_lambda;
            delta_u(iter+1)=delta_u(iter)+ddelta_u;
            
            %fig2 = plot(u(step)+delta_u(iter+1),(lambda(step)+delta_lambda(iter+1))*F_ext,'go','Markersize',5,'Markerfacecolor','g'); hold on
            %fig4 = plot(u(step)+delta_u(iter+1),F_inter(u(step)+delta_u(iter+1),k_bar),'ro','Markersize',5,'Markerfacecolor','r'); hold on
            Residual = -F_inter(u(step)+delta_u(iter+1),k_bar)+(lambda(step)+delta_lambda(iter+1))*F_ext;
        end
        %fig2 = plot(u(step)+delta_u(iter+1),(lambda(step)+delta_lambda(iter+1))*F_ext,'go','Markersize',5,'Markerfacecolor','g'); hold on
        fig4 = plot(u(step)+delta_u(iter+1),F_inter(u(step)+delta_u(iter+1),k_bar),'bo','Markersize',5,'Markerfacecolor','b'); hold on
        iter = iter+1;
    end
    step = step+1;
        
    %legend([fig1 fig2 fig4],'Analytical solution','Arc length point','Arc length prediction','location','best');
    legend([fig1 fig4],'Analytical solution','Arc length prediction','location','best');
    
end


%% function definition
function y=F_inter(x,k_bar)
y=x + 1.5*x.^2 + 0.5*x.^3 + k_bar * x;
end
function y=k_stif(x,k_bar)
y = (1 + 3 * x + 1.5 * x.^2) + k_bar;
end
function delta_lambda=lambda_solve(a1,a2,a3,delta_u,ddelta_u,delta_lambda,F_ext)
if (a1==0 && a2==0)
    disp('No roots found')
elseif (a1==0 && a2~=0)
    delta_lambda=-a3./a2;
else
    fac=a2*a2-4*a1*a3;
    if fac<0
       disp('No root found, return')
    end
    lambda1=(-a2+sqrt(fac))/(2*a1);
    dot1 = ddelta_u*delta_u+lambda1*delta_lambda*F_ext*F_ext;
    lambda2=(-a2-sqrt(fac))/(2*a1);
    dot2 = ddelta_u*delta_u+lambda2*delta_lambda*F_ext*F_ext;
    if dot2>=dot1
        delta_lambda=lambda2;
    else
        delta_lambda=lambda1;
    end
end
end
