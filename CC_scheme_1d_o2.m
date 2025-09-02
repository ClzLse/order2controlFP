function [deltap,deltam,Cp,Cm,Dp,Dm] = CC_scheme_1d_o2(B,contr,D,dD,x)
% Auxiliary function for the structure-preserving scheme. Using the notations 
% in the paper, the outputs are delta_{i+1/2} (deltap), delta_{i-1/2} (deltam),
% A_{i+1/2} (Cp), A_{i-1/2} (Cm), D(v_{i+1/2}) (Dp) and D(v_{i-1/2}) (Dm).
        
    N = length(x); 
    
    dx = x(2)-x(1);
    xp = x(:)+dx/2;
    xm = x(:)-dx/2;
    lambdap = zeros(N,1);
    lambdam = zeros(N,1);
    Cp = zeros(N,1);
    Cm = zeros(N,1);
    Intp = zeros(N,1);
    Intm = zeros(N,1);
    deltap = zeros(N,1);
    deltam = zeros(N,1);
     
    Bp = B(xp')-interp1(x,contr,xp,'cubic');
    Bm = B(xm')-interp1(x,contr,xm,'cubic');

    Intp(1:N-1) = dx*(Bp(1:N-1) + dD(xp(1:N-1)))./D(xp(1:N-1));
    Intm(2:N) = dx*(Bm(2:N) + dD(xm(2:N)))./D(xm(2:N));

    % Bp = B(xp,c);
    % Bm = B(xm,c);
    Cp(1:N-1) = Bp(1:N-1) + dD(xp(1:N-1));%D(xp(1:N-1)).*Intp(1:N-1)/dx;
    Cm(2:N) = Bm(2:N) + dD(xm(2:N));%D(xm(2:N)).*Intm(2:N)/dx;
    
    lambdap(1:N-1) = Intp(1:N-1);%dx.*Cp(1:N-1)./D(xp(1:N-1));
    lambdam(2:N) =  Intm(2:N);%dx.*Cm(2:N)./D(xm(2:N));
    
    deltap(1:N-1) = 1./(lambdap(1:N-1)) + 1./(1-exp(lambdap(1:N-1)));
    deltam(2:N) = 1./(lambdam(2:N)) + 1./(1-exp(lambdam(2:N)));
    toll = 1e-10;
    I = find(abs(lambdap)<toll);
    deltap(I) = Bp(I)<0;%(-sign(Bp(I))+1)/2;
    I = find(deltap==Inf);
    deltap(I) = Bp(I)<0;%(-sign(Bp(I))+1)/2;
    I = find(abs(lambdam)<toll);
    deltam(I) = Bm(I)<0;%(-sign(Bm(I))+1)/2;
    I = find(deltam==Inf);
    deltam(I) = Bm(I)<0;% (-sign(Bm(I))+1)/2;

    Dp = D(xp);
    Dm = D(xm);
    Dp = Dp(:);
    Dm = Dm(:);

   