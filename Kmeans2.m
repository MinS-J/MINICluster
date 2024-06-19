function [K2, clusterIndices] = Kmeans2(transformedData, iteration)
fh4 = @(X,K2)(kmeans(X,K2));
eva = evalclusters(transformedData,fh4,"CalinskiHarabasz","KList",1:iteration);
clear fh4
K2 = eva.OptimalK;
clusterIndices = eva.OptimalY;

%%

figure;

scatter(transformedData(:,1), transformedData(:,2),'k.'); 

legend('data')

xlabel('x1'); 

ylabel('x2');


hold on;

% Labeling the clusters

for i = 1:K2

scatter(transformedData(clusterIndices ==i,1), transformedData(clusterIndices ==i,2), 'filled'); 

end
end