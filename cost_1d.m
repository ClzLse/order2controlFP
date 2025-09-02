function J = cost_1d(V,V_obj,u,gamma,f,dv,dt)
% Value of the cost dunctional (J). The inputs are the grid (V), the objective
% value (V_obj), the control (u), the penalization of the cost (gamma), the
% values of the function f (f), the space step (dv) and the time step (dt).  
    J = 0.5*sum(sum((repmat(abs(V-V_obj).^2,1,size(f,2)) + gamma.*(u).^2).*f))*dv*dt;