function Dv1 = grad_1d(v,dx)
% First order derivatives with accuracy 2

    Dv1 = v;
    Dv1(1,:) = (v(2,:) - v(1,:))/(dx);
    Dv1(2:end-1,:) = (v(3:end,:) - v(1:end-2,:))/(2*dx);
    Dv1(end,:) = (v(end,:) - v(end-1,:))/dx;
