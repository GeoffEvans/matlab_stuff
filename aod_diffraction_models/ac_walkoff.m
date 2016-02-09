function [ walkoff ] = ac_walkoff( rotation )

rotation = rotation / 180 * pi;
dTh = 1e-9;
V = teo2.find_v_ac_min( pi/2+rotation+[0,dTh], pi/4+[0,0]);
dV = V(2) - V(1);

walkoff = atan(dV / V(2) / dTh) * 180 / pi;

end

