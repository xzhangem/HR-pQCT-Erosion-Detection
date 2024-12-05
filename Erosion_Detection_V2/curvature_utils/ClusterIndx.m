function cluster_indx = ClusterIndx(cluster_vol, K_num)
%% Parameter Setting
% Input:
%   cluster_vol: the numbered erosion volume
%   K_num: the number of erosion in the volume
% Output:
%   cluster_indx: a cell of size (K_num, 1), K-th element of the cell is
%       the coeffecient set of the K-th erosion with format (point_num, 3)

cluster_num = zeros(K_num, 1);
[vx, vy, vz] = size(cluster_vol);

for i = 1 : vx
    for j = 1 : vy
        for k = 1 : vz
            if cluster_vol(i,j,k) ~= 0
                indx = cluster_vol(i,j,k);
                cluster_num(indx, 1) = cluster_num(indx, 1) + 1;
            end
        end
    end
end

cluster_indx = cell(K_num, 1);
K_count = ones(K_num, 1);

for i = 1 : K_num
   cluster_indx{i} = zeros(cluster_num(i,1), 3);
end

for i = 1 : vx
    for j = 1 : vy
        for k = 1 : vz
            if cluster_vol(i,j,k) ~= 0
                K = cluster_vol(i,j,k);
                cluster_indx{K}(K_count(K,1), :) = [i, j, k];
                K_count(K, 1) = K_count(K, 1) + 1;
            end
        end
    end
end


end