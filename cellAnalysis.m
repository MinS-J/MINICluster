function [frequency, MEAN, SEM] = cellAnalysis(conditionidententity, cluster, clusternum)
for ii = 1:length(conditionidententity)
    Ttemp = conditionidententity{ii};
    idx = find(isnan(Ttemp.('decay time')));
    Ttemp(idx,:) = [];
    conditionidententity{ii} = Ttemp;
end

cellidentity = [];
for ii = 1:length(conditionidententity)
    cellidentity = [cellidentity, ii.*ones(1,height(conditionidententity{ii}))];
end

cell_ = cell(length(conditionidententity),1);
for ii = 1:length(conditionidententity)
    cell_{ii} = cluster(cellidentity == ii);
end

clusterdist = cell(1,length(conditionidententity));
for ii = 1:length(conditionidententity)
    cell_cluster = [];
    for jj = 1:clusternum
        cell_cluster(jj,1) = nnz(cell_{ii}==jj);
    end
    clusterdist{ii} = cell_cluster;
end

timewindow = 150;
clusterdist_mat = cell2mat(clusterdist);
frequency = clusterdist_mat./timewindow;
MEAN = mean(frequency,2);
STD = std(frequency,0,2);
SEM = STD/sqrt(size(frequency,2));
end