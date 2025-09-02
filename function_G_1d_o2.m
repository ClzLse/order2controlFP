function [g,operatore_Q] = function_G_1d_o2(B,contr,D,dD,x,f,dt)

    [deltap,deltam,Cp,Cm,Dp,Dm] = CC_scheme_1d_o2(B,contr,D,dD,x);
    N = length(x); 
    dx = x(2) - x(1);
    
    FvP = zeros(size(f));
    FvM = FvP;

    FvP(1:N-1) = Cp(1:N-1).*((1-deltap(1:N-1)).*f(2:N) + ...
        deltap(1:N-1).*f(1:N-1)) + Dp(1:N-1).*(f(2:N)-f(1:N-1))./dx;
    FvM(2:N) = Cm(2:N).*((1-deltam(2:N)).*f(2:N) + ...
        deltam(2:N).*f(1:N-1)) + Dm(2:N).*(f(2:N)-f(1:N-1))./dx;

    operatore_Q = (FvP - FvM)./dx;
    g = f + dt*operatore_Q;
    g(g<0) = 0;
