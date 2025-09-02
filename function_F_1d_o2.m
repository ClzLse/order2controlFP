function [F,g,operatore_Q1,operatore_Q2] = function_F_1d_o2(f,B,contr,D,dD,x,dt)
% Auxiliary function for the structure-preserving scheme. The input is the 
% function at the previous time step (f), the drift (B) the control
% (contr), the diffusion function (D), the derivative of the diffusion
% function (dD), the space grid (x) and the time step (dt).

    [g,operatore_Q1] = function_G_1d_o2(@(xx)B(xx,f),contr,D,dD,x,f,dt);
    % if ~isempty(find(g<0))
    %     keyboard
    % end
    N = length(x);
    dx = x(2)-x(1);
    
    % semi-implicit scheme
    [deltap,deltam,Cp,Cm,Dp,Dm] = ...
        CC_scheme_1d_o2(@(xx)B(xx,g),contr,D,dD,x);
        
    c = 1 + dt/(2*dx)*((Dp+Dm)/dx + Cm.*(1-deltam) - Cp.*deltap);
    b = -dt/(2*dx)*(Cp.*(1-deltap) + Dp/dx);
    % bnew = b(1:end-1);
    d = dt/(2*dx)*(deltam.*Cm - Dm/dx);
    % dnew = d(2:end);
    
    % M = diag(c) + diag(bnew,+1) + diag(dnew,-1);
    % F = M\reshape(f+dt*operatore_Q1/2,N*nc,1);% Tridiag(reshape(a,1,N*nc),reshape(b,1,N*nc),reshape(c,1,N*nc),reshape(f,N*nc,1));
    F = Tridiag(d,c,b,f+dt*operatore_Q1/2);%F = M\(f+dt*operatore_Q1/2);
    F = F(:);
    % F = reshape(F,N,nc);
    
    % if ~isempty(find(F<0))
    %     keyboard
    % end

    FvP = zeros(size(F));
    FvM = FvP;
    

    % deltap = reshape(deltap,N,nc);
    % deltam = reshape(deltam,N,nc);
    % Cp = reshape(Cp,N,nc);
    % Cm = reshape(Cm,N,nc);
    % Dp = reshape(Dp,N,nc);
    % Dm = reshape(Dm,N,nc);
    
    FvP(1:N-1) = Cp(1:N-1).*((1-deltap(1:N-1)).*F(2:N) + ...
        deltap(1:N-1).*F(1:N-1)) + Dp(1:N-1).*(F(2:N)-F(1:N-1))./dx;
    FvM(2:N) = Cm(2:N,:).*((1-deltam(2:N)).*F(2:N) + ...
        deltam(2:N).*F(1:N-1)) + Dm(2:N).*(F(2:N)-F(1:N-1))./dx;

    operatore_Q2 = (FvP - FvM)./dx;

 