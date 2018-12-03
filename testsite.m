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
%R = 6371000;             % Erdradius in m
k = 50;                  % Wärmeleitfähigkeit in W/m K

%Auflösung und Dimensionen
yr = 60*60*24*365.5;     %Jahr in Sekunden
t_end = 200000*yr;
dt = 1000*yr;                  
t_begin = 1000*yr;
dr = 1000;
q_cmb = 0.0005;     % W/m^2

%Gravitationsbeschleunigung berechnen
for r = 1000 : dr : r_cmb
    rho(r/dr,1) = rho_o;
    g(r/dr,1) = 4*pi*G*rho(r/dr,1)*r/3;
end

N = 2000;

% preallocate the arrays
%T = NaN(N + 200 ,t_end/dt);
%T_nonlayer = zeros(N,t_end/dt);
%T_layer = zeros(N,t_end/dt);
%r1 = zeros(t_end/dt);

%meshgrid für Zeit --> damit man später contourplot machen kann
for i = 1:1:t_end/dt-2
    for ii = 1:1:1997
    t_plot(ii,i) = i;
    end
end

%initiales adiabatisches T-Profil
grid_initial = linspace(0,r_cmb,N)';
T(:,1) = T_cmb*exp((-alpha*G*2*pi*rho_i*(grid_initial.^2 - r_cmb*r_cmb))/(3*cp));    
         

        
 %% Zeitschleife
for t = t_begin : dt : t_end
        %Suche nach erstem Auftreten von r1 an CMB + inital adiabatic temperature profile
        if t > t_begin && exist('r1') == 0
            
        T(:,t/dt) = T_cmb*exp((-alpha*G*2*pi*rho_i*(grid_initial.^2 - r_cmb*r_cmb))/(3*cp));
        r1(t/dt + 1) = r1_detect_cmb(q_cmb, r_cmb, t, dr, g, T, cp, alpha);
        
        end
        
       %Wenn r1 existiert und einen Wert ungleich NaN hat
        if exist('r1')== 1 && isnan(r1(t/dt -1 )) ~= 1 && r1(t/dt) ~= 0
 
        %% creating two equidistant grids for regions on top and below of r1
         grid_nonlayer = linspace(0,r1(t/dt), r1(t/dt)/dr)';
         grid_layer = linspace(r1(end),r_cmb, N - r1(t/dt)/dr)';
         
         %dr_nonlayer = grid_nonlayer(3) - grid_nonlayer(2);
         if size(grid_layer,1) > 1
         dr_layer = grid_layer(2) - grid_layer(1);
         end         
         
        %% creating new grids if T must be interpolated + interpolation
         if t/dt >2 && exist('T(r1/dr,t/dt)', 'var') == 0 && r1(t/dt)/dr - round(r1(t/dt))/dr ~= 0
             
         N_grid_interp = floor((r_cmb - r1(end))/1000); %Anzahl Gridpunkte für Interpolierten T-Vektor
                                                        % proportional zu Abstand CMB - r1
        
         T_layer_clean = snip(T_layer(:,t/dt-1) , '0'); %Vektor(t/dt -1) von Nullen säubern
         %T_layer_clean(1,t/dt) = T1(end);    %update T1 at new r1      
         
         new_grid = linspace(r1(end), r_cmb, N_grid_interp)'; 
         
         old_grid = linspace(r1(end - 1), r_cmb, size(T_layer_clean,1))';
         
         T_interp = spline(old_grid, T_layer_clean, new_grid);  %interpolation
         end
         
         % wenn keine Interpolation nötig weil r1 integer, setze T_interp
         if exist('T(r1/dr,t/dt)', 'var') == 0 && r1(t/dt)/dr - round(r1(t/dt))/dr == 0
           T_interp = T;  
         end
         
        %% Integral rho*T*dV von 0 bis r_cmb und Teilen durch 4Pi*r1^2*k --> A für T1 & p
         A = (4/3)*pi*cp*rho_o*T_cmb*r_cmb^3*(1 + (2*r_cmb*2*pi*alpha*G*rho_o)/(3*cp))/(4*pi*r1(end)^2*k);
         
        %% T1 und P berechnen für WLG
         if (r1(t/dt)/dr - round(r1(t/dt))/dr) == 0 %wenn r1 integer
         p = P_calc2(r1(end),t, T, T_interp, dr, dt, A, alpha, cp, rho_o);
         T1(t/dt + 1) = T1_calc(r1(t/dt),t, dr, dt, A, p, T);        
         end
         
         if r1(t/dt)/dr - round(r1(t/dt))/dr ~= 0  %wenn r1 kein integer
         p = P_calc2(r1(end), t, T, T_interp, dr_layer, dt, A, alpha, cp, rho_o);
         T1(t/dt + 1) = T1_calc(r1(t/dt),t, dr_layer, dt, A, p, T, T_interp);        
         end
         
         %% adiabatic temperature profile on lower grid if r1 and T1 are existing
         T_nonlayer(1:size(grid_nonlayer,1),t/dt) = T1(end)*exp((-alpha*G*2*pi*rho_i*(grid_nonlayer.^2 - r1(end).^2))/(3*cp));    
         
         %% Wärmeleitungsgleichung (für konduktive Schicht wenn r1 existiert) und r1(t+dt)  

         %T_layer(r1(t/dt)/dr:end,t/dt) = WLG(T(r1/dr:r_cmb/dr,t/dt),r_cmb, dt, q_cmb, T1 ,r1, rho_o, cp, k);
         if r1(t/dt)/dr - round(r1(t/dt))/dr == 0
         T_layer(:,t/dt) = WLG(T(r1(t/dt)/dr:r_cmb/dr,t/dt),r_cmb, dt,dr, q_cmb, T1 ,r1(t/dt), rho_o, cp, k);
         end
         
         if r1(t/dt)/dr - round(r1(t/dt))/dr ~= 0
         T_layer(1:1:N_grid_interp,t/dt) = WLG(T_interp,r_cmb, dt, dr_layer, q_cmb, T1 ,r1(t/dt), rho_o, cp, k);
         end      
         
        %% merge spatial grids
        if exist('new_grid') == 0
        grid_sum(1:size(vertcat(grid_nonlayer, grid_layer),1),t/dt) = vertcat(grid_nonlayer, grid_layer);
        end
        if exist('new_grid') == 1 && size(new_grid,1) > 2
        grid_sum(1:size(vertcat(grid_nonlayer, new_grid),1),t/dt) = vertcat(grid_nonlayer, new_grid);
        end
         
         %% merge the two T vectors
        T_nonlayer_clean = snip(T_nonlayer(:,t/dt) , '0');  % alle Nullen aus Vektor löschen 
        T(1:size(vertcat(T_nonlayer_clean , T_layer(:,t/dt))),t/dt) = vertcat(T_nonlayer_clean , T_layer(:,t/dt));
        
        %% Update new r1
         if r1(t/dt)/dr - round(r1(t/dt))/dr == 0
          r1(t/dt+1) = r1_update(r1(end), p, dr);
         end
         
         if r1(t/dt)/dr - round(r1(t/dt))/dr ~= 0
          r1(t/dt+1) = r1_update(r1(end), p, dr_layer);         
         end
        
        end
       disp('timestep:')
       disp(t/dt)
       if exist('p','var') ==1
       disp(p)
       end
end

 disp('Calculations completed without error.')
 
 %% Plotting stuff
 %[X , Y] = meshgrid(1:1:18,grid_sum(1:1997,3:end));
 %contourf(t_plot, grid_sum(1:1997,3:end), T(1:1997,3:end))
 xlabel('timestep');
 ylabel('depth');
 colorbar
 
 figure(2)
 plot(1:1:t_end/dt-1, T1(2:end))
 xlabel('timestep in [kYr]');
 ylabel('T1 in [K]');
 figure(3)
 plot(1:1:t_end/dt, r1(2:end))
 xlabel('timestep in [kYr]');
 ylabel('r1 in [m]');
 
 