function phi = LSEdgeSmooth(phi_0, vol, mu, lambda, alpha, epsilon, timestep, iter)
% Input: 
%   phi_0: Intial levelset function got from the distance function of the
%       binary contour;
%   vol: Original binary volume got from last step;
%   mu: Weight of the distance regularization term;
%   lambda: Weight of the length term;
%   alpha: Weight of the data fitting term;
%   epsilon: Width of the Dirac Delta function and Heaviside function;
%   timestep: time step;
%   iter: number of iterations;
% Output:
%   phi: The updated leveset function according to Energy functional 
%       E = \int alpha * (vol - H(phi))^2 + lambda * |\downtriangle H(phi)|
%           + mu * R(phi)

small_number = 1e-9;
phi = phi_0;
vol = NeumannBoundCond(vol);
for k = 1 : iter
    phi = NeumannBoundCond(phi);
    [phi_x, phi_y, phi_z] = gradient(phi);
    s = sqrt(phi_x.^2 + phi_y.^2 + phi_z.^2);
    Nx = phi_x ./ (s + small_number);
    Ny = phi_y ./ (s + small_number);
    Nz = phi_z ./ (s + small_number);
    m_curvature = Div(Nx, Ny, Nz);
    dist_regu = DistReg(phi);
    dirac_phi = SmoothDirac(phi, epsilon);
    heavi_phi = SmoothHeavi(phi, epsilon);
    edge_term = dirac_phi .* m_curvature;
    data_term = dirac_phi .* (heavi_phi - vol);
    phi = phi + timestep * (mu * dist_regu + lambda * edge_term - alpha * data_term);
end
end

function f = Div(x, y, z)
[x_x, ~, ~] = gradient(x);
[~, y_y, ~] = gradient(y);
[~, ~, z_z] = gradient(z);
f = x_x + y_y + z_z;
end

function f = DistReg(phi)
[phi_x, phi_y, phi_z] = gradient(phi);
s = sqrt(phi_x.^2 + phi_y.^2 + phi_z.^2);
a = (s>=0) & (s<=1);
b = (s>1);
ps = a .* sin(2 * pi * s) / (2 * pi) + b .* (s-1);
dps=((ps~=0) .* ps + (ps==0)) ./ ((s~=0) .* s + (s==0));
f = Div(dps .* phi_x - phi_x, dps .* phi_y - phi_y, dps .* phi_z - phi_z) + 4 * del2(phi);
end

function f = SmoothDirac(phi, epsilon)
f = (1 / (2 * epsilon)) * (1 + cos(pi * phi / epsilon));
b = (phi <= epsilon) & (phi >= -epsilon);
f = f .* b;
end

function f = SmoothHeavi(phi, epsilon)
f = (1 / 2) * (1 + phi / epsilon + (1 / pi) * sin(pi * phi / epsilon));
b = (phi <= epsilon) & (phi >= -epsilon);
a = (phi > epsilon);
c = (phi < -epsilon);
f = f .* b + 1. * a + 0. * c;
end

function f = int3D(x)
f = sum(sum(sum(x)));
end

function g = NeumannBoundCond(f)
[nx, ny, nz] = size(f);
g = f;
g([1 nx], [1 ny], [1 nz]) = g([3 nx-2], [3, ny-2], [3, nz-2]);
g([1 nx], 2:end-1, 2:end-1) = g([3 nx-2], 2:end-1, 2:end-1);
g(2:end-1, [1 ny], 2:end-1) = g(2:end-1, [3 ny-2], 2:end-1);
g(2:end-1, 2:end-1, [1 nz]) = g(2:end-1, 2:end-1, [3 nz-2]);
end

