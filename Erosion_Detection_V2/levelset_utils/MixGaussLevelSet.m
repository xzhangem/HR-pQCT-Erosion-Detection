function phi = MixGaussLevelSet(phi_0, vol, thres_seg, g, mu, lambda, alpha, beta, gamma, epsilon, timestep, iter, gauss_iter)
% Input: 
%   phi_0: inital level set function in this stage
%   vol: double volume for updating mean for external energy
%   thres_seg: segmentation from threshold levelset for initialization  
%   g: 3D surface indicator function
%   mu: weight of the distance regularization term
%   lambda: weight of the weighted length term 
%   alpha: weight of the weighted area term
%   beta: weight of the weighted internal Gaussian mixture
%   gamma: weight of the weighted external Gaussian
%   epsilon: width of Dirac Delta function and Heaviside function
%   timestep: time step
%   iter: number of iterations
%   gauss_iter: number of iterations for Gaussian mixture model

% Output:
%   phi: the updated level set function

vol = double(vol);
phi = double(phi_0);
bone_mean = int3D(vol .* thres_seg) / int3D(thres_seg);
bone_var = int3D((vol - bone_mean).^2 .* thres_seg) / int3D(thres_seg);
[vx, vy, vz] = gradient(g);
[lx, ly, lz] = size(vol);
w_fore = 0.5;
w_bck = 0.5;

for k = 1 : iter
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
    area_term = dirac_phi .* g;
    edge_term = dirac_phi .* (vx .* Nx + vy .* Ny + vz .* Nz) + dirac_phi .* g .* curvature;
    
    c_ex = (int3D(vol .* (1 - heavi_phi))) / (int3D(1 - heavi_phi));
    sigma_ex = (int3D((vol - c_ex).^2 .* (1 - heavi_phi))) / (int3D(1 - heavi_phi));
    
    if k == 1
        mean_fore = bone_mean; mean_bck = c_ex;
        var_fore = bone_var; var_bck = sigma_ex;
    end
    
    % Gaussian mixture model update in a loop
    for p = 1 : gauss_iter
        w_gauss_fore = w_fore * exp(-(vol - mean_fore).^2 / (2 * var_fore + small_number)) / (sqrt(var_fore) + small_number);
        w_gauss_bck = w_bck * exp(-(vol - mean_bck).^2 / (2 * var_bck + small_number)) / (sqrt(var_bck) + small_number);
    
        gamma_fore = w_gauss_fore ./ (w_gauss_fore + w_gauss_bck);
        gamma_bck = w_gauss_bck ./ (w_gauss_fore + w_gauss_bck);
    
        mean_fore = int3D(gamma_fore .* vol .* heavi_phi) / int3D(gamma_fore .* heavi_phi);
        mean_bck = int3D(gamma_bck .* vol .* heavi_phi) / int3D(gamma_bck .* heavi_phi);
    
        var_fore = int3D(gamma_fore .* (vol - mean_fore).^2 .* heavi_phi) / int3D(gamma_fore .* heavi_phi);
        var_bck = int3D(gamma_bck .* (vol - mean_bck).^2 .* heavi_phi) / int3D(gamma_bck .* heavi_phi);
    
        w_fore = int3D(gamma_fore .* heavi_phi) / int3D(heavi_phi);
        w_bck = 1 - w_fore;        
    end
    
    
    ext_term = (vol - c_ex).^2 / (2 * (sigma_ex + small_number)) + 0.5 * log(sigma_ex);
    gauss_mixture = w_fore * exp(-(vol - mean_fore).^2 / (2 * var_fore + small_number)) / (sqrt(var_fore) + small_number) + ...
        w_bck * exp(-(vol - mean_bck).^2 / (2 * var_bck + small_number)) / (sqrt(var_bck) + small_number);
    inter_term = (log(gauss_mixture));
    
    
    
    phi = phi + timestep * (mu * dist_regu + lambda * edge_term - alpha * area_term - beta * inter_term + ...
        gamma * ext_term);
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


