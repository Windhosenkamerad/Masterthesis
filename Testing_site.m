close all;
clear all;
clc;

% q = grad(T) * k;

%% Initiale Parameter
cp = 800;                %in J/K

%R = 2440000;            %Planetenradius in m
rho_o = 7000;            %Dichte äußerer Kern in kg/m3
rho_i = 8000;            %Dichte des inneren Kerns in kg/m3

%r_i = R - 5100000;      %Radius innerer Kern
r_cmb = 2*10^6;          %Radius CMB

T_cmb = 1800;            % in K
alpha = 9.0*10^-5;       % Therm. Diffusivität in 1/K
G = 6.67*10^-11;         % Gravitationskonstante
R = 6371000;             % Erdradius in m
k = 50;                  % Wärmeleitfähigkeit in W/m K

%Auflösung und Dimensionen
yr = 60*60*24*365.5;     %Jahr in Sekunden
t_end = 90*yr;
dt = 10*yr;                  
t_begin = 10*yr;
dr = 1000;
q_cmb = 0.0000000014;     % W/m^2

%Gravitationsbeschleunigung berechnen
for r = 1000 : dr : r_cmb
    rho(r/dr,1) = rho_o;
    g(r/dr,1) = 4*pi*G*rho(r/dr,1)*r/3;
end

N = 2000;

% preallocate the T-arrays
T = zeros(N,t_end/dt);
T_nonlayer = zeros(N,t_end/dt);
T_layer = zeros(N,t_end/dt);

%initiales adiabatisches T-Profil
grid_initial = linspace(0,r_cmb,N)';
T(:,1) = T_cmb*exp((-alpha*G*2*pi*rho_i*(grid_initial.^2 - r_cmb*r_cmb))/(3*cp));    
         

        
 %% Zeitschleife
for t = t_begin : dt : t_end
        %Suche nach erstem Auftreten von r1 an CMB + inital adiabatic temperature profile
        if t > t_begin && exist('r1(t/dt)','var') == 0
            
        T(:,t/dt) = T_cmb*exp((-alpha*G*2*pi*rho_i*(grid_initial.^2 - r_cmb*r_cmb))/(3*cp));
        r1(t/dt) = r1_detect_cmb(q_cmb, r_cmb, t, dr, g, T, cp, alpha);
        
        end
        
       %Wenn r1 existiert und einen Wert ungleich NaN hat
        if exist('r1(t)','var')== 1 && isnan(r1(t/dt)) ~= 1
 
        % creating two equidistant grids for regions on top and below of r1   
         grid_nonlayer = linspace(0,round(r1(t/dt)), r1(t/dt)/dr)';
         grid_layer = linspace(round(r1(t/dt)),r_cmb, N - r1(t/dt)/dr)';
         dr_nonlayer = grid_nonlayer(3) - grid_nonlayer(2);
         dr_layer = grid_layer(3) - grid_layer(2);
         
         if t/dt >2 && exist('T(r1/dr,t/dt)', 'var') == 0
         new_grid = linspace(r1(t/dt), r_cmb, 3)'; 
         % interpolating T on "uneven" grid to get an "even" grid
         T_spline = spline(grid_layer,T_layer(r1(t/dt)-p*dr^2:end,t/dt),new_grid)';
         end
         
        %Integral rho*T*dV von 0 bis r_cmb und Teilen durch 4Pi*r1^2*k --> A für T1 & p
         A = (4/3)*pi*cp*rho_o*T_cmb*r_cmb^3*(1 + (2*r_cmb*2*pi*alpha*G*rho_o)/(3*cp))/(4*pi*r1(t/dt)^2*k);
         
        %% T1 und P berechnen für WLG
          U = -t/dt*20000;
         %p = P_calc(r1(t/dt),t, r_cmb, T, dr, dt, A, alpha, cp, rho_o);
         T1 = T1_calc(r1(t/dt),t, r_cmb, dr, dt, A, p, T);
%          if t/dt >= 2
%          T1 = T(r1/dr,t/dt) + 30*t/dt;
%          end
        
         
         T_nonlayer(1:size(grid_nonlayer,1),t/dt) = T1*exp((-alpha*G*2*pi*rho_i*(grid_nonlayer.^2 - r1(t/dt).^2))/(3*cp));  %adiabatic temperature profile on lower grid if r1 and T1 are existing  
         
         
        %% Wärmeleitungsgleichung (für konduktive Schicht wenn r1 existiert) und r1(t+dt)  

         %T_layer(:,(t+dt)/dt) = WLG(T(r1/dr:r_cmb/dr,t/dt),r_cmb, dt, q_cmb, T1 ,r1, rho_o, cp, k);
         T_layer(r1(t/dt)/dr:end,t/dt) = WLG(T(r1(t/dt)/dr:r_cmb/dr,t/dt),r_cmb, dt, q_cmb, T1 ,r1(t/dt), rho_o, cp, k);
        
         % merge the two T vectors
         T_nonlayer(r1(t/dt)/dr,t/dt) = 0;
         T(:,t/dt) = T_nonlayer(:,t/dt) + T_layer(:,t/dt);
        
        r1(t/dt+1) = r1(t/dt) + p*dr^2; %r1 für nächsten Schritt berechnen
        %r1(t/dt +1) = r1(t/dt) - 2543.1;
        end
       disp('timestep:')
       disp(t/dt)
       
       if exist('r1(t/dt)','var') ~= 0
       disp(r1(t/dt)/dr)
       disp(p)
       end
end

disp('Calculations completed without error.')
 contourf(T)
 xlabel('timestep');
 ylabel('depth');
 colorbar