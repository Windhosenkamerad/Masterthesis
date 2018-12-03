function [ S ] = SUMME( x, r )
%UNTITLED2 Summary of this function goes here
%   Sums up the elements of vector x from element 1 to r
    for i = 1:1:r
        
        if i == 1
            S = 0;
        end
        
        S = S + x(i);
        
    end

end


%         syms k
%         f = k;
%         V = subs(f,k,1:r);
%         S_sum = sum(V);