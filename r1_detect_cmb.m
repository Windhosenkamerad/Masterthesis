function [r1] = r1_detect_cmb(q_cmb, r_cmb, t, dr, g, T, cp, alpha)

%   Check if the adiabatic conditions are given or not, then set r1

if exist('r1','var') == 0  

    %dT_norm = (T(r_cmb/dr,t) - T(r_cmb/dr - 1,t))/dr;
    dT_adia = g(r_cmb/dr,1)*alpha*T(r_cmb/dr,1)/cp;
    
    if  q_cmb < dT_adia
        
       r1 = r_cmb - 3*dr;
       
    end
    
    if q_cmb >= dT_adia
        
       r1 = NaN;
       
    end
    
    
end
    
end

