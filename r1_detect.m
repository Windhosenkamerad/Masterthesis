function [r1] = r1_detect(r_i ,r_cmb, t, dr, dt, g, T, cp, alpha)

%   Check if the adiabatic conditions are given or not, then set r1

for r=r_i:dr:r_cmb
    
    if (T(r-r_i+2,(t-dt)/dt) - T(r-r_i+1,(t-dt)/dt))/dr < g(r-r_i+1,(t-dt)/dt)*alpha*T(r-r_i+1,(t-dt)/dt)/cp
        
       r1 = r;
       
       %Wenn Bedingung das erste mal erfüllt ist, aufhören weiterzusuchen
         if r_cmb - r1 > 0
        
             break;
       
         end
    
       
    end
    
    %Garantieren einer Schichtmächtigkeit von mindestens 3dr wenn Bedingung
    %erfüllt ist
    if r_cmb - r1 < 3*dr 
                                
        r1 = r_cmb - 3*dr;
        
    end
        
end

