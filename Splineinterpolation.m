
r1_floor = floor(r1/dr);
dr_new = (r_cmb - r1)/r1 * dr;

dr_start = r1/dr - r1_floor;
%old_grid = r1_floor*dr:dr:size(T(r1_floor:r_cmb/dr,1),1)*dr;
old_grid = r1_floor*dr:dr:r_cmb;
TT = T(r1_floor:r_cmb/dr,1);
%new_grid = dr_start:dr_new:size(T(r1_floor:r_cmb/dr,1),1)*dr;
new_grid = r1:dr_new:r_cmb;
T_spline = spline(old_grid,TT,new_grid);
plot(old_grid,TT,'red',new_grid,T_spline,'o')

size(T_spline)