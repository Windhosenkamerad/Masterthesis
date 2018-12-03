function [ p ] = P_calc2(r1,t, T, T_interp, dr, dt, A, alpha, cp, rho)

G = 6.67*10^-11; 
g = 4/3*pi*G*rho*r1;

if r1/dr - round(r1)/dr ~= 0
%T0 = T_interp(1,1);   % T((r1-dr)/dr,t/dt);     %für Tspline
T1 = T_interp(1,1);   % T(r1)
T2 = T_interp(2,1);   % T(r1 + dr)
T3 = T_interp(3,1);   % T(r1 + 2dr)
end

if round(r1) - r1 == 0
T1 = T(r1/dr,t/dt);
T2 = T((r1 + dr)/dr,t/dt);
T3 = T((r1 + 2*dr)/dr,t/dt);
%T0 = T((r1-dr)/dr,t/dt);
end

% C = (A*T1 - T0)/(dt*(T1 - T0) - A);
% B = -(A*T1 - T0)/(dt*(T1 - T0) - A);
%C = (A*T2 - T1)/(dt*(T2 - T1) - A);
%B = -(A*T2 - T1)/(dt*(T2 - T1) - A);

M = dr*(T2 - T3/8 - T1*9/8);
N = T3/4 - T2 + 5*T1/4 - rho*T1*dr^2*alpha*4*pi*G/3*cp;
R = (1/dr)*(2*T2 + T3/2 - 3*T1/2) + alpha*g*T1/cp;
  
 p = - N/(2*M) - sqrt(N.^2 / (4*M.^2) - R/M);

end