clear
close all

set(0,'DefaultAxesTickLabelInterpreter','latex')
set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultAxesFontSize',28)
set(0,'DefaultLineLineWidth',1);

a = -1;
b = 1;
T = 4;

NN = 3;
dv = 2/2^NN;   % space step
dt = dv/10;    % time step
Nt = round(T/dt);   % number of time steps
T = Nt*dt;

v = a+dv/2:dv:b-dv/2;

nv = length(v);

delta = .1;
alpha = 1;
eps = 5; % diffusion coefficientb(sigma^2/2) 
sigma = sqrt(2*eps);
gamma = .001;   % penalization of the control in the cost

f0 = .5*ones(length(v),1);
f0 = f0/(dv*sum(f0));   

[Vt,Tt] = ndgrid(v,0:dt:T);

D = @(x) sigma^2*(1-x.^2)/2; % diffusion function
dD = @(x) -sigma^2*x;   % derivative of the diffusion function
P = @(x,y) alpha*(abs(x'-y)<=delta); % interaction function


Efp_arg = @(x,ff) -(P(x,v).*(-(x'-v)))*ff;
Efp = @(x,ff) (Efp_arg(x,ff))*dv;   % FP drift
Ehjb = @(x,ff) -Efp(x,ff);      % HJB drift
Q_arg = @(x,ff,dpsi) (P(v,x).*(x'-v))*(ff.*dpsi);
Q = @(x,ff,dpsi) (Q_arg(x,ff,dpsi))*dv;     % integral term appearing in the HJB
Sfun = @(x) 1;
% initialization
f_time = repmat(f0,1,Nt+1);
Du =  zeros(size(f_time));
contr = zeros(size(f_time));
psi_time = zeros(size(f_time));
psi_time_old = psi_time;

x_d = 0.3;      % goal opinion in the cost
Maxiter = 500;
cost_vector = zeros(1,Maxiter);

k = 1;
mk = 1;
lambda = 0.1;

DU_Lold = gamma.*contr + Du;
contr_n = contr;

err_cost = 10;

%first iteration
[psi_time,Du] = HJB_ordine2_1d(dt,dv,Nt,v,f_time,Ehjb,D,Q,gamma,x_d,contr,Du);   
f_time_old = f_time;
lambda_app = lambda;
contr = contr - lambda_app*(gamma.*contr + Du);
DU_L = gamma.*contr + Du;
f_time = time_semi_implicit_CC_1d_o2(f0,Efp,contr,D,dD,v,dt,Nt);%FP_ordine1_2d(dt,Nt,V,C,Du,E,f0,D,gamma,v,c,nv,nc,contr,tipo_iter);
cost_vector(k) = cost_1d(v(:),x_d,contr,gamma,f_time,Sfun,dv,dt);
contr_nm1 = contr_n;
contr_n = contr;
old_cost = cost_vector(k);
k = k+1;

% subsequent iterations
while (err_cost > 1e-5 && k<=Maxiter)
    psi_time_old = psi_time;
    controld = contr;
    [psi_time,Du] = HJB_ordine2_1d(dt,dv,Nt,v,f_time,Ehjb,D,Q,gamma,x_d,contr,Du);   
    f_time_old = f_time;
    lambda_app = abs(sum(reshape(contr_n - contr_nm1,1,nv*(Nt+1)).*reshape(DU_L - DU_Lold,1,nv*(Nt+1)))./(sum(reshape(DU_L - DU_Lold,1,nv*(Nt+1)).^2)));
    contr = contr - lambda_app*(gamma.*contr + Du);
    contr_nm1 = contr_n;
    contr_n = contr;
    f_time = time_semi_implicit_CC_1d_o2(f0,Efp,contr,D,dD,v,dt,Nt);%FP_ordine1_2d(dt,Nt,V,C,Du,E,f0,D,gamma,v,c,nv,nc,contr,tipo_iter);
    cost_vector(k) = cost_1d(v(:),x_d,contr,gamma,f_time,Sfun,dv,dt);
    iter_min = 0;

    % loop in case the Barzilai-Borwein fails
    while cost_vector(k) > cost_vector(k-1) && iter_min < 100
        lambda_app = lambda_app/2;
        contr = controld;
        contr = contr - lambda_app*(gamma.*contr + Du);
        contr_nm1 = contr_n;
        contr_n = contr;
        f_time = time_semi_implicit_CC_1d_o2(f0,Efp,contr,D,dD,v,dt,Nt);%FP_ordine1_2d(dt,Nt,V,C,Du,E,f0,D,gamma,v,c,nv,nc,contr,tipo_iter);
        cost_vector(k) = cost_1d(v(:),x_d,contr,gamma,f_time,Sfun,dv,dt);
        iter_min = iter_min + 1;

    end

    err_cost = abs(old_cost - cost_vector(k));
    DU_Lold = DU_L;
    DU_L = gamma.*contr + Du;
    old_cost = cost_vector(k);
    k=k+1;
end

cost_vector = cost_vector(1:k-1);
set(0,'DefaultLineLineWidth',4);
semilogy(1:k-1,cost_vector,'o')
xlabel('$k$'), ylabel('$\mathcal{J}(\cdot;f_0)$')
