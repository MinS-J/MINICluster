function [K, clusterIndices] = Kmeans(merged_data, dim)
fh4 = @(X,K)(kmeans(X,K));
eva = evalclusters(merged_data,fh4,"CalinskiHarabasz","KList",1:100);
clear fh4
K = eva.OptimalK;
clusterIndices = eva.OptimalY;

figure;
if dim == 2
    scatter(merged_data(:,1), merged_data(:,2),'k.'); 
    legend('data')
    xlabel('x1'); 
    ylabel('x2');
    hold on;
else
    scatter3(merged_data(:,1), merged_data(:,2), merged_data(:,3),'k.'); 
    legend('data')
    xlabel('x1'); 
    ylabel('x2');
    zlabel('x3');
    hold on;
end

% Labeling the clusters
if dim == 2
    for i = 1:K
        scatter(merged_data(clusterIndices ==i,1), merged_data(clusterIndices ==i,2), 'filled'); 
    end
else
    for i = 1:K
        scatter3(merged_data(clusterIndices ==i,1), merged_data(clusterIndices ==i,2), merged_data(clusterIndices ==i,3), 'filled'); 
    end
end
end
