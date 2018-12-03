function [ r1_new ] = r1_update(r1, p, dr )

r1_new = r1 + p*dr^2;

end

