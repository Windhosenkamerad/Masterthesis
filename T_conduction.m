function [ T ] = T_conduction( r_1, r_cmb, dr, dt, t, rhog, cp, lambda, T)

%   Calculation of the 1D - Heat transport equation from r_1 to r_cmb in
%   steps of dr

for r=r_1:dr:r_cmb
    
%Wenn Wärmequelle an unterster Schicht existiert    
        %  if r == r_start
        %            T(r,t+1) = dt*lambda*((T(r+1,t)- 2*T(r,t) + T(r-1,t))/ rho(r,t)*cp(r,t)*2*dr) + (dt / (rho(r,t)*cp(r,t))*(1+ T(r,t)));
        %         end
            
        if r > r_1 && r < r_cmb
          T(r/dr,t+1) = dt*lambda*((T(r+1/dr,t) - 2*T(r/dr,t) + T(r-1/dr,t))/ rhog*cp*2*dr) + (dt / (rhog*cp)*(T(r/dr,t)));
        end
        
        if r == r_cmb
            T(r/dr,t+1) = dt*lambda*((T(r+1/dr,t) - 2*T(r/dr,t) + T(r-1/dr,t))/ rhog*cp*2*dr) + (dt / (rhog*cp)*(-q_out + T(r/dr,t)));

        end
end
end

