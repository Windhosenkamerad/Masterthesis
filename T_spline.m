function [ T_spline ] = T_spline( r1, r_cmb, dr, T )
%   interpolate the vector T to get the entry at a non-integer position of
%   T


r1_floor = floor(r1/dr); %find nearest integer rounded down at entry T(r1/dr) if r1/dr is non-integer

dr_new = (r_cmb - r1)/100;

old_grid = r1_floor*dr:dr:r_cmb;
TT = T(r1_floor:r_cmb/dr,1);

new_grid = r1 - dr_new:dr_new:r_cmb;

T_spline = spline(old_grid,TT,new_grid)';

end

