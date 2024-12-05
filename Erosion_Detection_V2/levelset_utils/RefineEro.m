function erosion_region = RefineEro(erosion_ellp, joint)

timestep = 5; % 1
mu = 0.2 / timestep; % weight of the distance regularization term
iter_dtr = 20; % iteration of the global distribution levelset
lambda = 6; % weight of the weighted length term
alpha = 0.3; % weight of the weighted area term
epsilon = 2; % width of Dirac Delta function and Heaviside function



ellp = single(erosion_ellp > 0);
srfc = joint;
ori_srfc = joint;
se_0 = strel('sphere', 3);
ellp = imdilate(ellp, se_0);

srfc = single((ellp + srfc) > 0);

srfc = imfill(srfc);
phi_0 = SignDistance(srfc);
phi = LSEdgeSmooth(phi_0, srfc, mu, lambda, alpha, epsilon, timestep, iter_dtr);

vol = single(phi > 0);
res = single((vol - ori_srfc) ~= 0);
se_1 = strel('sphere', 2);
res = imopen(res, se_1);
se_2 = strel('sphere', 3);
res = imclose(res, se_2);
[con_lab, ~] = bwlabeln(res, 26);
erosion_region = zeros(size(con_lab));
mask = single(ellp > 0);
main_vol = mask .* con_lab;
unique_vec = unique(main_vol);

for u = 1 : length(unique_vec)
    if unique_vec(u) ~= 0
        erosion_region = unique_vec(u) * single(con_lab == unique_vec(u)) + erosion_region;
    end
end


end