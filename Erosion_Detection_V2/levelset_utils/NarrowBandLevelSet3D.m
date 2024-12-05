function phi = NarrowBandLevelSet3D(phi_0, g, mu, lambda, alpha, epsilon, timestep, iter, r)
% This is the implementation of narrow-band "Distance Regularized Level Set Evolution and Its
% Application to Image Segmentation" 

% Input: 
%   phi_0: inital level set function in this stage
%   g: 3D surface indicator function
%   mu: weight of the distance regularization term
%   lambda: weight of the weighted length term 
%   alpha: weight of the weighted area term
%   epsilon: width of Dirac Delta function and Heaviside function
%   timestep: time step
%   iter: number of iterations
%   r: radium of the neghborhood
%Output:
%   phi: the updated narrow band level set function

small_number = 1e-9;
%[X, Y, Z] = size(g);
[vx, vy, vz] = gradient(g);
phi = phi_0;
[indx_x, indx_y, indx_z] = NarrowBand3D(phi, r);
it = 1;
while (it <= iter)
    phi_old = phi;
    indx_old_x = indx_x; indx_old_y = indx_y; indx_old_z = indx_z;
    [num, ~] = size(indx_x);
    % Find the corresponding narrow-banded part in g. 
    % Find the corresponding narrow-banded part in gradient matrix of g.  
    % Find the left and right neighbor of the corresponding narrow-banded
    % part in gradient matrix of g for computing distance regularitzation.
    g_vector  =zeros(num, 1);
    gx_vector = zeros(num, 1);
    gy_vector = zeros(num, 1);
    gz_vector = zeros(num, 1);
    %gx_left = zeros(num, 1);
    %gy_left = zeros(num, 1);
    %gz_left = zeros(num, 1);
    %gx_right = zeros(num, 1);
    %gy_right = zeros(num, 1);
    %gz_right = zeros(num, 1);
    for i = 1 : num
        g_vector(i) = g(indx_x(i), indx_y(i), indx_z(i));
        gx_vector(i) = vx(indx_x(i), indx_y(i), indx_z(i));
        gy_vector(i) = vy(indx_x(i), indx_y(i), indx_z(i));
        gz_vector(i) = vz(indx_x(i), indx_y(i), indx_z(i));
        %{
        % Intialize the left vector of the gradient matrix of g 
        if indx_x(i) == 1
            gx_left(i) = vx(indx_x(i), indx_y(i), indx_z(i));
        else
            gx_left(i) = vx(indx_x(i)-1, indx_y(i), indx_z(i));
        end
        
        if indx_y(i) == 1
            gy_left(i) = vy(indx_x(i), indx_y(i), indx_z(i));
        else
            gy_left(i) = vy(indx_x(i), indx_y(i)-1, indx_z(i));
        end
        
        if indx_z(i) == 1
            gz_left(i) = vz(indx_x(i), indx_y(i), indx_z(i));
        else
            gz_left(i) = vz(indx_x(i), indx_y(i), indx_z(i)-1);
        end
        % Initialize the right vector of the gradient matrix of g
        if indx_x(i) == X
            gx_right(i) = vx(indx_x(i), indx_y(i), indx_z(i));
        else
            gx_right(i) = vx(indx_x(i)+1, indx_y(i), indx_z(i));
        end
        
        if indx_y(i) == Y
            gy_right(i) = vy(indx_x(i), indx_y(i), indx_z(i));
        else
            gy_right(i) = vy(indx_x(i), indx_y(i)+1, indx_z(i));
        end
        
        if indx_z(i) == Z
            gz_right(i) = vz(indx_x(i), indx_y(i), indx_z(i));
        else
            gz_right(i) = vz(indx_x(i), indx_y(i), indx_z(i)+1);
        end
        %}
    end
    
    [gx_c, gy_c, gz_c, gx_l, gy_l, gz_l, gx_r, gy_r, gz_r] = GradientInfo3D(phi, indx_x, indx_y, indx_z);
    
    gc_s = sqrt(gx_c.^2 + gy_c.^2 + gz_c.^2);
    Ngx_c = gx_c ./ (gc_s + small_number);
    Ngy_c = gy_c ./ (gc_s + small_number);
    Ngz_c = gz_c ./ (gc_s + small_number);
    
    gl_s = sqrt(gx_l.^2 + gy_l.^2 + gz_l.^2);
    Ngx_l = gx_l ./ (gl_s + small_number);
    Ngy_l = gy_l ./ (gl_s + small_number);
    Ngz_l = gz_l ./ (gl_s + small_number);
    
    gr_s = sqrt(gx_r.^2 + gy_r.^2 + gz_r.^2);
    Ngx_r = gx_r ./ (gr_s + small_number);
    Ngy_r = gy_r ./ (gr_s + small_number);
    Ngz_r = gz_r ./ (gr_s + small_number);
    
    % Compute Distance Regularization given the gradient information here
    %a_RD = (gc_s >= 0) & (gc_s <= 1);
    %b_RD = (gc_s > 1);
    %ps_RD = a_RD .* sin(2 * pi * gc_s) / (2 * pi) + b_RD .* (gc_s - 1);
    %dps_RD = ((ps_RD ~= 0) .* ps_RD + (ps_RD == 0)) ./ ((gc_s ~= 0) .* gc_s + (gc_s == 0));
    
    
    a_RD_l = (gl_s >= 0) & (gl_s <= 1);
    b_RD_l = (gl_s > 1);
    ps_RD_l = a_RD_l .* sin(2 * pi * gl_s) / (2 * pi) + b_RD_l .* (gl_s - 1);
    dps_RD_l = ((ps_RD_l ~= 0) .* ps_RD_l + (ps_RD_l == 0)) ./ ((gl_s ~= 0) .* gl_s + (gl_s == 0));

    
    a_RD_r = (gr_s >= 0) & (gr_s <= 1);
    b_RD_r = (gr_s > 1);
    ps_RD_r = a_RD_r .* sin(2 * pi * gr_s) / (2 * pi) + b_RD_r .* (gr_s - 1);
    dps_RD_r = ((ps_RD_r ~= 0) .* ps_RD_r + (ps_RD_r == 0)) ./ ((gr_s ~= 0) .* gr_s + (gr_s == 0));

    dis_reg = (dps_RD_r .* gx_r - dps_RD_l .* gx_l) / 2 + (dps_RD_r .* gy_r - dps_RD_l .* gy_l) / 2 ...
        + (dps_RD_r .* gz_r - dps_RD_l .* gz_l) / 2;
    
    curvature = (Ngx_r - Ngx_l) / 2 + (Ngy_r - Ngy_l) / 2 + (Ngz_r - Ngz_l) / 2;
    dirac_phi = NBSmoothDirac(phi, indx_x, indx_y, indx_z, epsilon);
    area_term = dirac_phi .* g_vector;
    edge_term = dirac_phi .* (gx_vector .* Ngx_c + gy_vector .* Ngy_c + gz_vector .* Ngz_c) ...
        + dirac_phi .* g_vector .* curvature;
    % Compute Loss and Update phi
    %Loss = 0.0;
    for j = 1 : num
        var = timestep * (mu * dis_reg(j) + lambda * edge_term(j) - alpha * area_term(j));
        %Loss = abs(Loss) + var;
        phi(indx_x(j), indx_y(j), indx_z(j)) = phi(indx_x(j), indx_y(j), indx_z(j)) + var;
    end
    
    %if (Loss < 0.5)
    %    fprintf('Loss is %s small enough to terminate!\n', Loss);
    %    break;
    %end
    
    % Update the narrow band index given based on the updated phi
    [indx_x, indx_y, indx_z] = NarrowBand3D(phi, r);
    % Assign values to new pixels on the narrow band
    same_indx = ismember(indx_x, indx_old_x) .* ismember(indx_y, indx_old_y) .* ismember(indx_z, indx_old_z);
    [num_new, ~] = size(indx_x);
    num = max(num, num_new);
    for k = 1 : num
        if (same_indx(k) == 1)
            if (phi_old(indx_old_x(k), indx_old_y(k), indx_old_z(k)) < 0)
                phi(indx_x(k), indx_y(k), indx_z(k)) = -r - 1;
            elseif (phi_old(indx_old_x(k), indx_old_y(k), indx_old_z(k)) > 0)
                phi(indx_x(k), indx_y(k), indx_z(k)) = r + 1;
            end
        end
    end
    it = it + 1;
end
end


function f = NBSmoothDirac(phi, indx_x, indx_y, indx_z, epsilon)
[num, ~] = size(indx_x);
f = zeros(num, 1);
for i = 1 : num
    f(i) = (1/2/epsilon) * (1 + cos(pi * phi(indx_x(i), indx_y(i), indx_z(i)) / epsilon));
end
b = (f <= epsilon) & (f >= epsilon);
f = f .* b;
end

function f = NBSmoothHeavi(phi, indx_x, indx_y, indx_z, epsilon)
[num, ~] = size(indx_x);
f = zeros(num, 1);
for i = 1 : num
    f(i) = (1 / 2) * (1 + phi(indx_x(i), indx_y(i), indx_z(i)) / epsilon + (1 / pi) * sin(pi * phi(indx_x(i), indx_y(i), indx_z(i)) ...
        / epsilon));
end
b = (f <= epsilon) & (f >= -epsilon);
a = (f > epsilon);
f = f .* b + 1 .* a;
end

