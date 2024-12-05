function [cluster_vol, K_num] = SpatialCluster(erosion, r)
%% Parameter Setting
% Input:
%   erosion: erosion detection 
%   r: radium for neighborhood searching
% Output:
%   cluster_vol: indexed erosion after spatial clustering
%   K_num: erosion numbers of under this predescribed r

[vx, vy, vz] = size(erosion);


cluster_mat = zeros(vx, vy, vz, 2);
cluster_mat(:, :, :, 1) = erosion; % (intensity_value, cluster)

K_num = 0;

for i = 1 + r : vx - r
    for j = 1 + r : vy - r
        for k = 1 + r : vz - r
            if cluster_mat(i,j,k,1) == 1
                if cluster_mat(i,j,k,2) == 0
                    s = sum(sum(sum(cluster_mat(i-r:i+r, j-r:j+r, k-r:k+r, 2))));
                    if s == 0
                        cluster_mat(i-r:i+r,j-r:j+r,k-r:k+r,2) = (K_num + 1) * ... 
                            cluster_mat(i-r:i+r,j-r:j+r,k-r:k+r,1);
                        K_num = K_num + 1;
                    else
                        sub_indx_vol = cluster_mat(i-r:i+r,j-r:j+r,k-r:k+r,2);
                        sub_indx_vol(find(sub_indx_vol == 0)) = [];
                        sub_mean = mean(sub_indx_vol);
                        K_value = mode(sub_indx_vol);
                        if K_value ~= sub_mean
                            fprintf('[%d, %d, %d]---mean: %f; mode: %f.\n', i,j,k,sub_mean,K_value);
                        end
                        cluster_mat(i-r:i+r,j-r:j+r,k-r:k+r,2) = K_value * ... 
                            cluster_mat(i-r:i+r,j-r:j+r,k-r:k+r,1);                        
                    end
                else
                    cluster_mat(i-r:i+r,j-r:j+r,k-r:k+r,2) = cluster_mat(i,j,k,2) * ...
                        cluster_mat(i-r:i+r,j-r:j+r,k-r:k+r,1);
                end
            end
        end
    end
end

cluster_vol = cluster_mat(:,:,:,2);


end