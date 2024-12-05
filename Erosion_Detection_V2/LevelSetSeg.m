function seg_result = LevelSetSeg(vol, mode, iter_num, input_PDE_step)

addpath('levelset_utils');
addpath('curvature_utils');

%% Parameter Setting
% Segmentation Parameter Part
timestep = 1; % 1
mu = 0.2 / timestep; % weight of the distance regularization term
iter_dtr = iter_num; % iteration of the global distribution levelset
lambda = 6; % weight of the weighted length term
alpha = 1.5; % weight of the weighted area term
beta = 1; % weight of the foreground distribution
gamma = 1; % weight of the background distribution
epsilon = 1; % width of Dirac Delta function and Heaviside function
sigma = 0.2; % variance of Gaussian filter for extracting edge
PDE_step = input_PDE_step; % timestep for PDE-based morphology

pdf_case = mode; % specified foreground distribution 

vol = AboveZero(vol);
vol = double(vol);

vol_smooth = imgaussfilt3(vol, sigma, 'FilterSize', 11);
[Vx, Vy, Vz] = gradient(vol_smooth);
f = Vx.^2 + Vy.^2 + Vz.^2;
g = 1 ./ (1 + f);

% Fixed Initialization for Global Distribution Level Set Segmentation
c0 = 2;
[vx, vy, vz] = size(vol);
initialLSF = - c0 * ones(vx, vy, vz);
coe = 0.15;
initialLSF(floor(vx * coe) :  vx - floor(vx * coe), floor(vy * coe) : ...
    vy - floor(vy * coe), floor(vz * coe) : vz - floor(vz * coe)) = c0;
phi = initialLSF;
    
    
switch(pdf_case)
    case 'Gamma'
        phi = GammaLevelSet(phi, vol, g, mu, lambda, alpha, beta, gamma, epsilon, timestep, iter_dtr);
    case 'InverseGauss'
        phi = InverseGaussLevelSet(phi, vol, g, mu, lambda, alpha, beta, gamma, epsilon, timestep, iter_dtr);
    case 'Lognormal'
        phi = LogNormalLevelSet(phi, vol, g, mu, lambda, alpha, beta, gamma, epsilon, timestep, iter_dtr);
    case 'Gumbel'
        phi = GumbelLevelSet(phi, vol, g, mu, lambda, alpha, beta, gamma, epsilon, timestep, iter_dtr);
    case 'Weibull'
        phi = WeibullLevelSet(phi, vol, g, mu, lambda, alpha, beta, gamma, epsilon, timestep, iter_dtr);
    case 'GEV'
        phi = GEVLevelSet(phi, vol, g, mu, lambda, alpha, beta, gamma, epsilon, timestep, iter_dtr);
    case 'Thresh'
        phi = phi_0;
    case 'Gauss'
        phi = LevelSet_3D(phi, vol, g, mu, lambda, alpha, beta, gamma, epsilon, timestep, iter_dtr);
end

contour = single(phi > 0);
%seg = IsoDiffusionFill(contour, PDE_step, 1);
%refine_seg = MorphFillSmooth(seg, 1, 3, 'cuboid');
%{
refine_seg = contour;
refine_seg = padarray(refine_seg, [200,200,200], 'symmetric', 'both');
[vx_p, vy_p, vz_p] = size(refine_seg);
refine_seg = imfill(refine_seg);
seg_result = refine_seg(1+200:vx_p-200, 1+200:vy_p-200, 1+200:vz_p-200);
%}

seg_result = contour;
end
