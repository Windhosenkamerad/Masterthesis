function [ p ] = P_calc(r1,t, T, T_interp, dr, dt, A, alpha, cp, rho)

G = 6.67*10^-11; 
g = 4/3*pi*G*rho*r1;

% if r1(t/dt)/dr - round(r1(t/dt))/dr ~= 0
% T0 = T_interp(1,1);   % T((r1-dr)/dr,t/dt);     %für Tspline
% T1 = T_interp(2,1);   % T(r1)
% T2 = T_interp(3,1);   % T(r1 + dr)
% T3 = T_interp(4,1);   % T(r1 + 2dr)
% end
% 
% if round(r1(t/dt)) - r1(t/dt) == 0
% T1 = T(r1(t/dt)/dr,t/dt);
% T2 = T((r1(t/dt) + dr)/dr,t/dt);
% T3 = T((r1(t/dt) + 2*dr)/dr,t/dt);
% T0 = T((r1(t/dt)-dr)/dr,t/dt);
% end

if r1/dr - round(r1)/dr ~= 0
T0 = T_interp(1,1);   % T((r1-dr)/dr,t/dt);     %für Tspline
T1 = T_interp(2,1);   % T(r1)
T2 = T_interp(3,1);   % T(r1 + dr)
T3 = T_interp(4,1);   % T(r1 + 2dr)
end

if round(r1) - r1 == 0
T1 = T(r1/dr,t/dt);
T2 = T((r1 + dr)/dr,t/dt);
T3 = T((r1 + 2*dr)/dr,t/dt);
T0 = T((r1-dr)/dr,t/dt);
end

C = (A*T1 - T0)/(dt*(T1 - T0) - A);
B = -(A*T1 - T0)/(dt*(T1 - T0) - A);

M = T2*dr - T3*dr/8 - 9*B*dr/8;
N = T3/4 - T2 + 5*B/4 + 3*C*dr/2 - C*alpha*g*T1*dr^2/cp - (B*alpha/cp)*(T1*4*pi*rho/3 + g*(T2-T1)/dr);%(B*dr^2*alpha*4*pi*rho*G*T1)/(3*cp);
R = 1/dr*(2*T2 + T3/2 - 3*B/2) + B*alpha*g*T1/cp;
  
% p = - N/(2*M) - sqrt(N^2 / (4*M^2) - R/M);
  p = - N/(2*M) + sqrt(N.^2 / (4*M.^2) - R/M);

end

