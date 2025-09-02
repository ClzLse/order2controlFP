function [solHJ,gradu] = HJB_ordine2_1d(dt,dv,Nt,v,f_time,E,D,Q,gamma,x_d,u,gradold)
% This function computes the solution to the HJB equation with the SL scheme,
% given the time step (dt), the space step (dx), the number of time steps (Nt),
% the grip points (v), the solution to the FP equation in the spacextime 
% grid (f_time), the drift (E), the diffusion (D), the nonlocal integral 
% term in the HJB equation (Q), the penalization of the control (gamma), 
% the objective opinion in the cost (x_d), the approximation of the control 
% (u) and the previous approximation of the gradient of the solution to the HJB
% (gradold).

    solHJ = zeros(size(f_time));
    gradu = zeros(size(f_time));
    
    DHJ = grad_1d(solHJ(:,end),dv);
    gradu(:,end) = DHJ; 
    for k = Nt:-1:1 
        y0 = HJ(v(:),solHJ(:,k+1),dt,D,E(v,f_time(:,k)),E(v,f_time(:,k+1)),u(:,k),u(:,k+1),@(xx)Q(xx,f_time(:,k+1),gradold(:,k+1)),gamma,x_d) + ...
            dt*0.5*(0.5*(v(:)-x_d).^2 + Q(v,f_time(:,k),DHJ) + .5*gamma*u(:,k).^2);
        solHJ(:,k) = y0;
        DHJ = grad_1d(y0,dv);
        gradu(:,k) = DHJ;
    end
end

% function that computes the characteristics
function y=HJ(V,fun,dt,D,H,Hp1,contr,contrp1,Q,gamma,x_d)
    xp1 = controlloNeumann_1d(V + (H + contr)*dt + sqrt(dt*6*D(V)),V,dt);
    xp1 = controlloNeumann_1d(V + dt*(H + contr + interp1(V,Hp1+contrp1,xp1,'cubic'))/2 + sqrt(dt*6*D(V)),V,dt);

    xp2 = controlloNeumann_1d(V + (H + contr)*dt - sqrt(dt*6*D(V)),V,dt);
    xp2 = controlloNeumann_1d(V + dt*(H + contr + interp1(V,Hp1+contrp1,xp2,'cubic'))/2 - sqrt(dt*6*D(V)),V,dt);

    xp3 = controlloNeumann_1d(V + (H + contr)*dt,V,dt);
    xp3 = controlloNeumann_1d(V + dt*(H + contr + interp1(V,Hp1+contrp1,xp3,'cubic'))/2,V,dt);

    integ_term = @(xx) 0.5*dt*(0.5*(xx-x_d).^2 + Q(xx) + .5*gamma*interp1(V,contr.^2,xx));
    
    y = (interp1(V,fun,xp1,'cubic') + integ_term(xp1))/6 + (interp1(V,fun,xp2,'cubic') + integ_term(xp2))/6 + 2*(interp1(V,fun,xp3,'cubic') + integ_term(xp3))/3;

end

