function [ T1_real ] = T1_calc(r1, t, dr, dt, A, p, T, T_interp)

if r1/dr - round(r1)/dr ~= 0
T1 = T_interp(1,1);   % T(r1)
T2 = T_interp(2,1);   % T(r1 + dr)
end

if round(r1) - r1 == 0
T1 = T(r1/dr,t/dt);
T2 = T((r1 + dr)/dr,t/dt);
end

dT = T2 - T1;

%B = - (A*T1)/((dt*dT/dr) - A);

%T1_real = B - (A*p*dr^2*(T1 - T0))/(dt*(T1 - T0) - A);
T1_real = ((-A/dt)*(T1 + p*dr^2*dT/dr))/(dT/dr - A/dt);
%T1_real = (B - (A*p*dr*dT)/(dt*dT/dr - A));
end

