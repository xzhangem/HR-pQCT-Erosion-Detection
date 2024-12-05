function phi = GammaLevelSet(phi_0, vol, g, mu, lambda, alpha, beta, gamma_v, epsilon, timestep, iter)
% Input: 
%   phi_0: inital level set function in this stage
%   vol: double volume for updating mean for external energy
%   template: template used in the following step to make it loyal to the
%       boundary
%   g: 3D surface indicator function
%   mu: weight of the distance regularization term
%   lambda: weight of the weighted length term 
%   alpha: weight of the weighted area term
%   beta: weight of the weighted internal Gamma
%   gamma: weight of the weighted external Gaussian
%   epsilon: width of Dirac Delta function and Heaviside function
%   timestep: time step
%   iter: number of iterations

%Output:
%   phi: the updated level set function

vol = double(vol);

% Ensure that the intensity for Inverse Gaussian is above 0 (>= 1). 
%min_val = min(min(min(vol)));

%vol = vol - min(min_val, 0) + 1;
%log_vol = log(vol);

phi = double(phi_0);
[vx, vy, vz] = gradient(g);
[lx, ly, lz] = size(vol);

for i = 1 : iter
    phi = NeumannBoundCond(phi);
    [phi_x, phi_y, phi_z] = gradient(phi);
    s = sqrt(phi_x.^2 + phi_y.^2 + phi_z.^2);
    small_number = 1e-7;
    Nx = phi_x ./ (s + small_number);
    Ny = phi_y ./ (s + small_number);
    Nz = phi_z ./ (s + small_number);
    curvature = Div(Nx, Ny, Nz);
    dist_regu = DistReg(phi);
    dirac_phi = SmoothDirac(phi, epsilon);
    heavi_phi = SmoothHeavi(phi, epsilon);
    %area_term = dirac_phi .* g;
    %edge_term = dirac_phi .* (vx .* Nx + vy .* Ny + vz .* Nz) + dirac_phi .* g .* curvature;
    edge_term = dirac_phi .* curvature;
    area_term = dirac_phi;   
    c_ex = (int3D(vol .* (1 - heavi_phi))) / (int3D(1 - heavi_phi));
    sigma_ex = (int3D((vol - c_ex).^2 .* (1 - heavi_phi))) / (int3D(1 - heavi_phi));
    ext_term = ((vol - c_ex).^2 / (2 * (sigma_ex + small_number)) + 0.5 * log(sigma_ex));
    
    
    % Calculate k and theta in shape-scale parameteration of Gamma
    % distribution. Use the modified closed MLE form.
    
    %in_mu = (int3D(vol .* heavi_phi)) / (int3D(heavi_phi));
    %in_lambda = 1 / (int3D(inverse_vol .* heavi_phi) / int3D(heavi_phi) - int3D(heavi_phi) / int3D(vol .* heavi_phi));
    

    %inter_term = dirac_phi .* (-log(in_lambda + small_number) + 3 * log(vol + small_number) + (in_lambda * (vol - in_mu).^2) ...
    %    ./ (in_mu^2 * (vol + small_number)));
    
    %indicator_phi = single(phi > 0);
    
    inner_shift = min(min(min(vol .* heavi_phi))) - 0.5;
    
    %fprintf('%f \n', inner_shift);
    
    inner_log_vol = LogVol(vol - inner_shift);
    
    %heavi_indicate = single(heavi_phi > 0);
    
    k = int3D((vol - inner_shift) .* heavi_phi) / (int3D((vol - inner_shift) .* inner_log_vol .* heavi_phi) - int3D(inner_log_vol .* heavi_phi) ...
        * int3D((vol - inner_shift) .* heavi_phi) / int3D(heavi_phi));
    
    theta = int3D((vol - inner_shift) .* inner_log_vol .* heavi_phi) / int3D(heavi_phi) - (int3D((vol - inner_shift) .* heavi_phi) * ...
        int3D(inner_log_vol .* heavi_phi) / (int3D(heavi_phi))^2);
    
    modify_k = k - (3 * k - 2 * k / (3 + 3*k) - 4 * k / (5 + 10 * k + 5 * k^2)) / int3D(heavi_phi);
    
    modify_theta = (int3D(heavi_phi) / (int3D(heavi_phi) - 1)) * theta;
    
    inter_term = (-(modify_k - 1) * inner_log_vol + (vol - inner_shift) / modify_theta + modify_k * ...
        log(modify_theta) + log(gamma(modify_k)));
    
    phi = phi + timestep * (mu * dist_regu + lambda * edge_term - alpha * area_term - beta * inter_term + ...
        gamma_v * ext_term);
    
    %{
    binary = single(phi > 0);
    phi_contour = single(bwperim(phi));
    unsign_phi = bwdist(phi_contour);
    phi = unsign_phi .* (2 * binary - 1);
    %}
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

function g = LogVol(x)
[vx, vy, vz] = size(x);
g = zeros(vx, vy, vz);

for i = 1 : vx
    for j = 1 : vy
        for k = 1 : vz
            if x(i, j, k) <= 0 
                g(i, j, k) = 0;
            else
                g(i, j, k) = log(x(i, j, k));
            end
        end
    end
end

end