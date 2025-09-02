function xp = controlloNeumann_1d(xp,x,dt)
% Function thta implements Neumann zero boundary conditiond. xp is the foot
% of the characteristic, x is the space grid, dt is the space step. The
% output is given by the reflected characteristic.
    xp(xp<x(1)) = x(1)+.1*sqrt(dt);
    xp(xp>x(end)) = x(end)-.1*sqrt(dt);
    %size(find(xp<x(1) | xp>x(end)))
