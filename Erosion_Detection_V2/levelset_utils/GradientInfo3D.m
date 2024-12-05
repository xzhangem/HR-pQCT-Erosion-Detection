function [Gx_c, Gy_c, Gz_c, Gx_l, Gy_l, Gz_l, Gx_r, Gy_r, Gz_r] = GradientInfo3D(phi, indx_x, indx_y, indx_z)
[Xv, Yv, Zv] = size(phi);
[num, ~] = size(indx_x);

Gx_c = zeros(num, 1); Gy_c = zeros(num, 1); Gz_c = zeros(num, 1);
Gx_l = zeros(num, 1); Gy_l = zeros(num, 1); Gz_l = zeros(num, 1);
Gx_r = zeros(num, 1); Gy_r = zeros(num, 1); Gz_r = zeros(num, 1);

for i = 1 : num

gx_c = (indx_x(i) == 1) .* (phi(indx_x(i)+1, indx_y(i), indx_z(i)) - phi(indx_x(i), indx_y(i), indx_z(i))) ...
    + (indx_x(i) == Xv) .* (phi(indx_x(i), indx_y(i), indx_z(i)) - phi(indx_x(i)-1, indx_y(i), indx_z(i))) ...
    + ((indx_x(i) < Xv) & (indx_x(i) > 1)) .* (phi(indx_x(i)+1, indx_y(i), indx_z(i)) - phi(indx_x(i)-1, indx_y(i), indx_z(i))) / 2;
Gx_c(i) = gx_c;

gy_c = (indx_y(i) == 1) .* (phi(indx_x(i), indx_y(i)+1, indx_z(i)) - phi(indx_x(i), indx_y(i), indx_z(i))) ...
    + (indx_y(i) == Yv) .* (phi(indx_x(i), indx_y(i), indx_z(i)) - phi(indx_x(i), indx_y(i)-1, indx_z(i))) ...
    + ((indx_y(i) < Yv) & (indx_y(i) > 1)) .* (phi(indx_x(i), indx_y(i)+1, indx_z(i)) - phi(indx_x(i), indx_y(i)-1, indx_z(i))) / 2;
Gy_c(i) = gy_c;

gz_c = (indx_z(i) == 1) .* (phi(indx_x(i), indx_y(i), indx_z(i)+1) - phi(indx_x(i), indx_y(i), indx_z(i))) ...
    + (indx_z(i) == Zv) .* (phi(indx_x(i), indx_y(i), indx_z(i)) - phi(indx_x(i), indx_y(i), indx_z(i)-1)) ...
    + ((indx_z(i) < Zv) & (indx_z(i) > 1)) .* (phi(indx_x(i), indx_y(i), indx_z(i)+1) - phi(indx_x(i), indx_y(i), indx_z(i)-1)) / 2;
Gz_c(i) = gz_c;


% Boundary condition for its left and right side is added in this way: when
% we use its left and right side to calculate central difference in general
% cases (right - left) / 2., then the boudary condition for this case is the
% same as the standard forward one besides a factor 2, which only change
% the scale but not affect the direction of the gradient. 


gx_l = (indx_x(i) == 1) .* (phi(indx_x(i)+1, indx_y(i), indx_z(i)) - phi(indx_x(i), indx_y(i), indx_z(i))) ...
    + (indx_x(i) == 2) .* (phi(indx_x(i), indx_y(i), indx_z(i)) - phi(indx_x(i)-1, indx_y(i), indx_z(i))) ...
    + ((indx_x(i) > 2) & (indx_x(i) <= Xv)) .* (phi(indx_x(i), indx_y(i), indx_z(i)) - phi(indx_x(i)-2, indx_y(i), indx_z(i))) / 2;
Gx_l(i) = gx_l;

gy_l = (indx_y(i) == 1) .* (phi(indx_x(i), indx_y(i)+1, indx_z(i)) - phi(indx_x(i), indx_y(i), indx_z(i))) ... 
    + (indx_y(i) == 2) .* (phi(indx_x(i), indx_y(i), indx_z(i)) - phi(indx_x(i), indx_y(i)-1, indx_z(i))) ...
    + ((indx_y(i) > 2) & (indx_y(i) <= Yv)) .* (phi(indx_x(i), indx_y(i), indx_z(i)) - phi(indx_x(i), indx_y(i)-2, indx_z(i))) / 2;
Gy_l(i) = gy_l;

gz_l = (indx_z(i) == 1) .* (phi(indx_x(i), indx_y(i), indx_z(i)+1) - phi(indx_x(i), indx_y(i), indx_z(i))) ...
    + (indx_z(i) == 2) .* (phi(indx_x(i), indx_y(i), indx_z(i)) - phi(indx_x(i), indx_y(i), indx_z(i)-1)) ...
    + ((indx_z(i) > 2) & (indx_z(i) <= Zv)) .* (phi(indx_x(i), indx_y(i), indx_z(i)) - phi(indx_x(i), indx_y(i), indx_z(i)-2)) / 2;
Gz_l(i) = gz_l;


gx_r = (indx_x(i) == Xv) .* (phi(indx_x(i), indx_y(i), indx_z(i)) - phi(indx_x(i)-1, indx_y(i), indx_z(i))) ...
    + (indx_x(i) == Xv-1) .* (phi(indx_x(i)+1, indx_y(i), indx_z(i)) - phi(indx_x(i), indx_y(i), indx_z(i))) ...
    + ((indx_x(i) >= 1) & (indx_x(i) < Xv-1)) .* (phi(indx_x(i)+2, indx_y(i), indx_z(i)) - phi(indx_x(i), indx_y(i), indx_z(i))) / 2;
Gx_r(i) = gx_r;

gy_r = (indx_y(i) == Yv) .* (phi(indx_x(i), indx_y(i), indx_z(i)) - phi(indx_x(i), indx_y(i)-1, indx_z(i))) ...
    + (indx_y(i) == Yv-1) .* (phi(indx_x(i), indx_y(i)+1, indx_z(i)) - phi(indx_x(i), indx_y(i), indx_z(i))) ...
    + ((indx_y(i) >= 1) & (indx_y(i) < Yv-1)) .* (phi(indx_x(i), indx_y(i)+2, indx_z(i)) - phi(indx_x(i), indx_y(i), indx_z(i))) / 2;
Gy_r(i) = gy_r;

gz_r = (indx_z(i) == Zv) .* (phi(indx_x(i), indx_y(i), indx_z(i)) - phi(indx_x(i), indx_y(i), indx_z(i)-1)) ...
    + (indx_z(i) == Zv-1) .* (phi(indx_x(i), indx_y(i), indx_z(i)+1) - phi(indx_x(i), indx_y(i), indx_z(i))) ...
    + ((indx_y(i) >= 1) & (indx_y(i) < Zv-1)) .* (phi(indx_x(i), indx_y(i), indx_z(i)+2) - phi(indx_x(i), indx_y(i), indx_z(i))) / 2;
Gz_r(i) = gz_r;

end

end
