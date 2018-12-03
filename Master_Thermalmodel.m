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

dr = 1000;
dV = 4*pi*dr^3/3;        % Volumen einer Schale mit Mächtigkeit dr

%Zeiteinstellungen
yr = 60*60*24*365.5;     %Jahr in Sekunden
t_end = 60*yr;
dt = 10*yr;                  
t_begin = 10*yr;

q_cmb = 0.00000014;     % W/m^2

%Dichteverteilung innerer und äußerer Kern
% for r = 1000:dr:r_cmb 
% %     if r > 0 && r <= r_i
% %         rho(r/dr,1) = rho_i;
% %     end
%     
% %     if r > r_i 
% %         rho(r/dr,1) = rho_o;
% %     end
%         rho(r/dr,1) = rho_o;
% 
% end

%Gravitationsbeschleunigung berechnen
for r = 1000 : dr : r_cmb
    rho(r/dr,1) = rho_o;
    g(r/dr,1) = 4*pi*G*rho(r/dr,1)*r/3;
%     if r <= r_i
%         rhog = rho_i;
%        g(r/dr,1) = 4*pi*G*rhog*r/3; 
%     end
%     
%     if r > r_i
%         rhog = rho_o;
%         g(r/dr,1) =  g(r_i/dr,1) + 4*pi*G*rhog*r/3;
%         %g(r/dr,1) = SUMME(rho,r_cmb/dr)*dV;
%     end
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
        if t > t_begin && exist('r1','var') == 0
            
        T(:,t/dt) = T_cmb*exp((-alpha*G*2*pi*rho_i*(grid_initial.^2 - r_cmb*r_cmb))/(3*cp));
            
        r1 = r1_detect_cmb(q_cmb, r_cmb, t, dr, g, T, cp, alpha);
        
        end
       
        if exist('r1','var')== 1 && isnan(r1) ~= 1  %wenn r1 existiert und einen Wert ungleich NaN hat

 
            
        %Integral rho*T*dV von 0 bis r_cmb und Teilen durch 4Pi*r1^2*k --> A für T1 & p
         A = (4/3)*pi*cp*rho_o*T_cmb*r_cmb^3*(1 + (2*r_cmb*2*pi*alpha*G*rho_o)/(3*cp))/(4*pi*r1^2*k);
        
        %% T1 und P berechnen für WLG
         p = P_calc(r1,t, r_cmb, T, dr, dt, A, alpha, cp, rho_o);
         T1 = T1_calc(r1,t, r_cmb, dr, dt, A, p, T);
        
        
        % creating two equidistant grids for regions on top and below of r1   
         %grid_low = linspace(0,r1, round(N - (r_cmb - r1/1000)));
         %grid_up = linspace(r1,r_cmb,round(r_cmb - r1/1000));
         grid_nonlayer = linspace(0,r1, N - (r_cmb - r1)/1000)';
         grid_layer = linspace(r1,r_cmb,(r_cmb - r1)/1000)';
         %grid_whole = [grid_nonlayer(1:end -1 ,1)' , grid_layer(1:end ,1)' ]';
         
         T_nonlayer(1:size(grid_nonlayer,1),t/dt) = T1*exp((-alpha*G*2*pi*rho_i*(grid_nonlayer.^2 - r1*r1))/(3*cp));  %adiabatic temperature profile on lower grid if r1 and T1 are existing  
         
         
        %% Wärmeleitungsgleichung (für konduktive Schicht wenn r1 existiert) und r1(t+dt)  

         %T_layer(:,(t+dt)/dt) = WLG(T(r1/dr:r_cmb/dr,t/dt),r_cmb, dt, q_cmb, T1 ,r1, rho_o, cp, k);
         T_layer(r1/dr:end,t/dt) = WLG(T(r1/dr:r_cmb/dr,t/dt),r_cmb, dt, q_cmb, T1 ,r1, rho_o, cp, k);
        
         T(:,t/dt) = [T_nonlayer(1:r1/dr-1,t/dt)' , T_layer(r1/dr:end,t/dt)']';  % merge the two T vectors
        
        %T(r1/dr:r_cmb/dr,t/dt) = T_layer(1:end, t/dt);
        r1 = r1 + p*dr^2; %r1 für nächsten Schritt berechnen
        
        end
       disp('timestep:')
       disp(t/dt)
end

disp('Calculations completed without error.')
 contourf(T)
 colorbar
%plot(1:1:size(g,1),g)