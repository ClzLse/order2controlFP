function f_time = time_semi_implicit_CC_1d_o2(f0,B,contr,D,dD,x,dt,Nt)
% This function approximates the solution of the FP equation with the 
% structure-preserving scheme. The input is the initial condition (f0), the 
% drift function (B), the approximation of the control (contr), the diffution
% function (D), the derivative of the diffusion function (dD), the space
% grid (x), the time step (dt) and the number of time steps (Nt).

    N = length(x);
    f_time = zeros(N,Nt+1);
    f_time(:,1) = f0;
    for n = 1:Nt
        [~,~,operatore_Q1,operatore_Q2] = function_F_1d_o2(f_time(:,n),B,contr(:,n),D,dD,x,dt);
        F_new = f_time(:,n) + dt/2*(operatore_Q1 + operatore_Q2);
        f_time(:,n+1) = F_new;
    end
