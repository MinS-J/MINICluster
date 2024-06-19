clear all
clc
close all
load('Serotonin_neuron.mat')
%% Condition 1

condition1 = fields(finalStruct.Serotonin_neuron_condition_1);
Condition1_DR_Excitatory = cell(0);
Condition1_MR_Excitatory = cell(0);
Condition1_DR_Inhibitory = cell(0);
Condition1_MR_Inhibitory = cell(0);

for ii = 1:length(fields(finalStruct.Serotonin_neuron_condition_1))
    Excitatory = fields(finalStruct.Serotonin_neuron_condition_1.(condition1{ii}).Excitatory);
    ExcitatoryDR = Excitatory(contains(Excitatory,'DR'));
    ExcitatoryMR = Excitatory(contains(Excitatory,'MR'));

    Inhibitory = fields(finalStruct.Serotonin_neuron_condition_1.(condition1{ii}).Inhibitory);
    InhibitoryDR = Inhibitory(contains(Inhibitory,'DR'));
    InhibitoryMR = Inhibitory(contains(Inhibitory,'MR'));

    for jj = 1:length(ExcitatoryDR)
        Condition1_DR_Excitatory{end+1,1} = finalStruct.Serotonin_neuron_condition_1.(condition1{ii}).Excitatory.(ExcitatoryDR{jj});
    end
    for jj = 1:length(ExcitatoryMR)
        Condition1_MR_Excitatory{end+1,1} = finalStruct.Serotonin_neuron_condition_1.(condition1{ii}).Excitatory.(ExcitatoryMR{jj});
    end
    for jj = 1:length(InhibitoryDR)
        Condition1_DR_Inhibitory{end+1,1} = finalStruct.Serotonin_neuron_condition_1.(condition1{ii}).Inhibitory.(InhibitoryDR{jj});
    end
    for jj = 1:length(InhibitoryMR)
        Condition1_MR_Inhibitory{end+1,1} = finalStruct.Serotonin_neuron_condition_1.(condition1{ii}).Inhibitory.(InhibitoryMR{jj});
    end
end
%% Condition1 - DR

[Condition1_DR_Excitatory_Amp, Condition1_DR_Excitatory_RT, Condition1_DR_Excitatory_DT] = getparameters(Condition1_DR_Excitatory);

[Condition1_DR_Inhibitory_Amp, Condition1_DR_Inhibitory_RT, Condition1_DR_Inhibitory_DT] = getparameters(Condition1_DR_Inhibitory);

% %% 2D plot - Condition1 - DR - Amp, RT
% fig1 = figure;
% plot(Condition1_DR_Excitatory_Amp,Condition1_DR_Excitatory_RT,'LineStyle','none','Marker','+','MarkerSize',1);
% 
% fig2 = figure;
% plot(Condition1_DR_Inhibitory_Amp,Condition1_DR_Inhibitory_RT,'LineStyle','none','Marker','+','MarkerSize',1);

%% Condition1 - MR

[Condition1_MR_Excitatory_Amp, Condition1_MR_Excitatory_RT, Condition1_MR_Excitatory_DT] = getparameters(Condition1_MR_Excitatory);

[Condition1_MR_Inhibitory_Amp, Condition1_MR_Inhibitory_RT, Condition1_MR_Inhibitory_DT] = getparameters(Condition1_MR_Inhibitory);

% %% 2D plot - Amp, RT
% fig3 = figure;
% plot(Condition1_MR_Excitatory_Amp,Condition1_MR_Excitatory_RT,'LineStyle','none','Marker','+','MarkerSize',1);
% 
% fig4 = figure;
% plot(Condition1_MR_Inhibitory_Amp,Condition1_MR_Inhibitory_RT,'LineStyle','none','Marker','+','MarkerSize',1);

%% Condition 2

condition2 = fields(finalStruct.Serotonin_neuron_condition_2);
Condition2_DR_Excitatory = cell(0);
Condition2_MR_Excitatory = cell(0);
Condition2_DR_Inhibitory = cell(0);
Condition2_MR_Inhibitory = cell(0);

for ii = 1:length(fields(finalStruct.Serotonin_neuron_condition_2))
    Excitatory = fields(finalStruct.Serotonin_neuron_condition_2.(condition2{ii}).Excitatory);
    ExcitatoryDR = Excitatory(contains(Excitatory,'DR'));
    ExcitatoryMR = Excitatory(contains(Excitatory,'MR'));

    Inhibitory = fields(finalStruct.Serotonin_neuron_condition_2.(condition2{ii}).Inhibitory);
    InhibitoryDR = Inhibitory(contains(Inhibitory,'DR'));
    InhibitoryMR = Inhibitory(contains(Inhibitory,'MR'));

    for jj = 1:length(ExcitatoryDR)
        Condition2_DR_Excitatory{end+1,1} = finalStruct.Serotonin_neuron_condition_2.(condition2{ii}).Excitatory.(ExcitatoryDR{jj});
    end
    for jj = 1:length(ExcitatoryMR)
        Condition2_MR_Excitatory{end+1,1} = finalStruct.Serotonin_neuron_condition_2.(condition2{ii}).Excitatory.(ExcitatoryMR{jj});
    end
    for jj = 1:length(InhibitoryDR)
        Condition2_DR_Inhibitory{end+1,1} = finalStruct.Serotonin_neuron_condition_2.(condition2{ii}).Inhibitory.(InhibitoryDR{jj});
    end
    for jj = 1:length(InhibitoryMR)
        Condition2_MR_Inhibitory{end+1,1} = finalStruct.Serotonin_neuron_condition_2.(condition2{ii}).Inhibitory.(InhibitoryMR{jj});
    end
end

%% Condition2 - DR

[Condition2_DR_Excitatory_Amp, Condition2_DR_Excitatory_RT, Condition2_DR_Excitatory_DT] = getparameters(Condition2_DR_Excitatory);

[Condition2_DR_Inhibitory_Amp, Condition2_DR_Inhibitory_RT, Condition2_DR_Inhibitory_DT] = getparameters(Condition2_DR_Inhibitory);

% %% 2D plot - Condition2 - DR - Amp, RT
% fig5 = figure;
% plot(Condition2_DR_Excitatory_Amp,Condition2_DR_Excitatory_RT,'LineStyle','none','Marker','+','MarkerSize',1);
% 
% fig6 = figure;
% plot(Condition2_DR_Inhibitory_Amp,Condition2_DR_Inhibitory_RT,'LineStyle','none','Marker','+','MarkerSize',1);

%% Condition2 - MR

[Condition2_MR_Excitatory_Amp, Condition2_MR_Excitatory_RT, Condition2_MR_Excitatory_DT] = getparameters(Condition2_MR_Excitatory);

[Condition2_MR_Inhibitory_Amp, Condition2_MR_Inhibitory_RT, Condition2_MR_Inhibitory_DT] = getparameters(Condition2_MR_Inhibitory);

% %% 2D plot - Condition2 - MR - Amp, RT
% fig7 = figure;
% plot(Condition2_MR_Excitatory_Amp,Condition2_MR_Excitatory_RT,'LineStyle','none','Marker','+','MarkerSize',1);
% 
% fig8 = figure;
% plot(Condition2_MR_Inhibitory_Amp,Condition2_MR_Inhibitory_RT,'LineStyle','none','Marker','+','MarkerSize',1);

%% Remove Nan

[Condition1_DR_Excitatory_Amp_Final, Condition1_DR_Excitatory_RT_Final, Condition1_DR_Excitatory_DT_Final] = removeNan( ...
    Condition1_DR_Excitatory_Amp, Condition1_DR_Excitatory_RT, Condition1_DR_Excitatory_DT);

[Condition1_MR_Excitatory_Amp_Final, Condition1_MR_Excitatory_RT_Final, Condition1_MR_Excitatory_DT_Final] = removeNan( ...
    Condition1_MR_Excitatory_Amp, Condition1_MR_Excitatory_RT, Condition1_MR_Excitatory_DT);

[Condition1_DR_Inhibitory_Amp_Final, Condition1_DR_Inhibitory_RT_Final, Condition1_DR_Inhibitory_DT_Final] = removeNan( ...
    Condition1_DR_Inhibitory_Amp, Condition1_DR_Inhibitory_RT, Condition1_DR_Inhibitory_DT);

[Condition1_MR_Inhibitory_Amp_Final, Condition1_MR_Inhibitory_RT_Final, Condition1_MR_Inhibitory_DT_Final] = removeNan( ...
    Condition1_MR_Inhibitory_Amp, Condition1_MR_Inhibitory_RT, Condition1_MR_Inhibitory_DT);


[Condition2_DR_Excitatory_Amp_Final, Condition2_DR_Excitatory_RT_Final, Condition2_DR_Excitatory_DT_Final] = removeNan( ...
    Condition2_DR_Excitatory_Amp, Condition2_DR_Excitatory_RT, Condition2_DR_Excitatory_DT);

[Condition2_MR_Excitatory_Amp_Final, Condition2_MR_Excitatory_RT_Final, Condition2_MR_Excitatory_DT_Final] = removeNan( ...
    Condition2_MR_Excitatory_Amp, Condition2_MR_Excitatory_RT, Condition2_MR_Excitatory_DT);

[Condition2_DR_Inhibitory_Amp_Final, Condition2_DR_Inhibitory_RT_Final, Condition2_DR_Inhibitory_DT_Final] = removeNan( ...
    Condition2_DR_Inhibitory_Amp, Condition2_DR_Inhibitory_RT, Condition2_DR_Inhibitory_DT);

[Condition2_MR_Inhibitory_Amp_Final, Condition2_MR_Inhibitory_RT_Final, Condition2_MR_Inhibitory_DT_Final] = removeNan( ...
    Condition2_MR_Inhibitory_Amp, Condition2_MR_Inhibitory_RT, Condition2_MR_Inhibitory_DT);


%% Merge data

Excitatory_Amp = [Condition1_DR_Excitatory_Amp_Final, Condition2_DR_Excitatory_Amp_Final, Condition1_MR_Excitatory_Amp_Final, Condition2_MR_Excitatory_Amp_Final];
Excitatory_RT = [Condition1_DR_Excitatory_RT_Final, Condition2_DR_Excitatory_RT_Final, Condition1_MR_Excitatory_RT_Final, Condition2_MR_Excitatory_RT_Final];
Excitatory_DT = [Condition1_DR_Excitatory_DT_Final, Condition2_DR_Excitatory_DT_Final, Condition1_MR_Excitatory_DT_Final, Condition2_MR_Excitatory_DT_Final];
Excitatory_merge = [Excitatory_Amp' Excitatory_RT' Excitatory_DT']; % Amplitude | Risetime | Decaytime


Inhibitory_Amp = [Condition1_DR_Inhibitory_Amp_Final, Condition2_DR_Inhibitory_Amp_Final, Condition1_MR_Inhibitory_Amp_Final, Condition2_MR_Inhibitory_Amp_Final];
Inhibitory_RT = [Condition1_DR_Inhibitory_RT_Final, Condition2_DR_Inhibitory_RT_Final, Condition1_MR_Inhibitory_RT_Final, Condition2_MR_Inhibitory_RT_Final];
Inhibitory_DT = [Condition1_DR_Inhibitory_DT_Final, Condition2_DR_Inhibitory_DT_Final, Condition1_MR_Inhibitory_DT_Final, Condition2_MR_Inhibitory_DT_Final];
Inhibitory_merge = [Inhibitory_Amp' Inhibitory_RT' Inhibitory_DT']; % Amplitude | Risetime | Decaytime

%% Dimensionality reduction

% 주성분 분석 - Excitatory
transformedData_Excitatory = dimReduction(Excitatory_merge);

% 주성분 분석 - Inhibitory
transformedData_Inhibitory = dimReduction(Inhibitory_merge);


%% 2D kmeans

[K_2D_Excitatory, clusterIndices_2D_Excitatory] = Kmeans(transformedData_Excitatory, 2);

[K_2D_Inhibitory, clusterIndices_2D_Inhibitory] = Kmeans(transformedData_Inhibitory, 2);

%% Label map

Excitatory_label = [ones(1,numel(Condition1_DR_Excitatory_Amp_Final)), 2.*ones(1,numel(Condition2_DR_Excitatory_Amp_Final)), 3.*ones(1,numel(Condition1_MR_Excitatory_Amp_Final)), 4.*ones(1,numel(Condition2_MR_Excitatory_Amp_Final))];
Inhibitory_label = [ones(1,numel(Condition1_DR_Inhibitory_Amp_Final)), 2.*ones(1,numel(Condition2_DR_Inhibitory_Amp_Final)), 3.*ones(1,numel(Condition1_MR_Inhibitory_Amp_Final)), 4.*ones(1,numel(Condition2_MR_Inhibitory_Amp_Final))];

%% Comparison
[Condition1_DR_excitatory_cluster,Condition2_DR_excitatory_cluster,Condition1_MR_excitatory_cluster,Condition2_MR_excitatory_cluster] = SEPARATION(clusterIndices_2D_Excitatory,Excitatory_label,K_2D_Excitatory);

[Condition1_DR_inhibitory_cluster,Condition2_DR_inhibitory_cluster,Condition1_MR_inhibitory_cluster,Condition2_MR_inhibitory_cluster] = SEPARATION(clusterIndices_2D_Inhibitory,Inhibitory_label,K_2D_Inhibitory);

% 1: Condition1 DR, 2: Condition2 DR, 3: Condition1 MR, 4: Condition2 MR

%% Cell comparison

[Condition1_DR_frequency_e, Mean_Condition1_DR_cell_Cluster_e, sem_Condition1_DR_cell_Cluster_e] = cellAnalysis(Condition1_DR_Excitatory,Condition1_DR_excitatory_cluster,K_2D_Excitatory);
[Condition2_DR_frequency_e, Mean_Condition2_DR_cell_Cluster_e, sem_Condition2_DR_cell_Cluster_e] = cellAnalysis(Condition2_DR_Excitatory,Condition2_DR_excitatory_cluster,K_2D_Excitatory);
[Condition1_MR_frequency_e, Mean_Condition1_MR_cell_Cluster_e, sem_Condition1_MR_cell_Cluster_e] = cellAnalysis(Condition1_MR_Excitatory,Condition1_MR_excitatory_cluster,K_2D_Excitatory);
[Condition2_MR_frequency_e, Mean_Condition2_MR_cell_Cluster_e, sem_Condition2_MR_cell_Cluster_e] = cellAnalysis(Condition2_MR_Excitatory,Condition2_MR_excitatory_cluster,K_2D_Excitatory);

[Condition1_DR_frequency_i, Mean_Condition1_DR_cell_Cluster_i, sem_Condition1_DR_cell_Cluster_i] = cellAnalysis(Condition1_DR_Inhibitory,Condition1_DR_inhibitory_cluster,K_2D_Inhibitory);
[Condition2_DR_frequency_i, Mean_Condition2_DR_cell_Cluster_i, sem_Condition2_DR_cell_Cluster_i] = cellAnalysis(Condition2_DR_Inhibitory,Condition2_DR_inhibitory_cluster,K_2D_Inhibitory);
[Condition1_MR_frequency_i, Mean_Condition1_MR_cell_Cluster_i, sem_Condition1_MR_cell_Cluster_i] = cellAnalysis(Condition1_MR_Inhibitory,Condition1_MR_inhibitory_cluster,K_2D_Inhibitory);
[Condition2_MR_frequency_i, Mean_Condition2_MR_cell_Cluster_i, sem_Condition2_MR_cell_Cluster_i] = cellAnalysis(Condition2_MR_Inhibitory,Condition2_MR_inhibitory_cluster,K_2D_Inhibitory);
%% ttest - DR between group - Excitatory

[h1,p1] = ttest_cell(Condition1_DR_frequency_e, Condition2_DR_frequency_e, Mean_Condition1_DR_cell_Cluster_e, Mean_Condition2_DR_cell_Cluster_e, ...
    sem_Condition1_DR_cell_Cluster_e, sem_Condition2_DR_cell_Cluster_e, K_2D_Excitatory);

[h2,p2] = ttest_cell(Condition1_MR_frequency_e, Condition2_MR_frequency_e, Mean_Condition1_MR_cell_Cluster_e, Mean_Condition2_MR_cell_Cluster_e, ...
    sem_Condition1_MR_cell_Cluster_e, sem_Condition2_MR_cell_Cluster_e, K_2D_Excitatory);

[h3,p3] = ttest_cell(Condition1_DR_frequency_e, Condition1_MR_frequency_e, Mean_Condition1_DR_cell_Cluster_e, Mean_Condition1_MR_cell_Cluster_e, ...
    sem_Condition1_DR_cell_Cluster_e, sem_Condition1_MR_cell_Cluster_e, K_2D_Excitatory);

%% ttest - DR between group - Inhibitory

[h4,p4] = ttest_cell(Condition1_DR_frequency_i, Condition2_DR_frequency_i, Mean_Condition1_DR_cell_Cluster_i, Mean_Condition2_DR_cell_Cluster_i, ...
    sem_Condition1_DR_cell_Cluster_i, sem_Condition2_DR_cell_Cluster_i, K_2D_Inhibitory);

[h5,p5] = ttest_cell(Condition1_MR_frequency_i, Condition2_MR_frequency_i, Mean_Condition1_MR_cell_Cluster_i, Mean_Condition2_MR_cell_Cluster_i, ...
    sem_Condition1_MR_cell_Cluster_i, sem_Condition2_MR_cell_Cluster_i, K_2D_Inhibitory);

[h6,p6] = ttest_cell(Condition1_DR_frequency_i, Condition1_MR_frequency_i, Mean_Condition1_DR_cell_Cluster_i, Mean_Condition1_MR_cell_Cluster_i, ...
    sem_Condition1_DR_cell_Cluster_i, sem_Condition1_MR_cell_Cluster_i, K_2D_Inhibitory);

%% Excitatory
figure;
scatter(transformedData_Excitatory(:,1), transformedData_Excitatory(:,2),'k.');
xlabel('x1');
ylabel('x2');
hold on;

for i = find(cell2mat(h4)==1)
    scatter(transformedData_Excitatory(clusterIndices_2D_Excitatory==i,1), transformedData_Excitatory(clusterIndices_2D_Excitatory==i,2), 'filled');
end

%% Inhibitory
figure;
scatter(transformedData_Inhibitory(:,1), transformedData_Inhibitory(:,2),'k.');
xlabel('x1');
ylabel('x2');
hold on;

for i = find(cell2mat(h6)==1)
    scatter(transformedData_Inhibitory(clusterIndices_2D_Inhibitory==i,1), transformedData_Inhibitory(clusterIndices_2D_Inhibitory==i,2), 'filled');
end

%% 3D kmeans

[K_3D_Excitatory, clusterIndices_3D_Excitatory] = Kmeans(Excitatory_merge,3);

[K_3D_Inhibitory, clusterIndices_3D_Inhibitory] = Kmeans(Inhibitory_merge,3);

%% Label map

Excitatory_label = [ones(1,numel(Condition1_DR_Excitatory_Amp_Final)), 2.*ones(1,numel(Condition2_DR_Excitatory_Amp_Final)), 3.*ones(1,numel(Condition1_MR_Excitatory_Amp_Final)), 4.*ones(1,numel(Condition2_MR_Excitatory_Amp_Final))];
Inhibitory_label = [ones(1,numel(Condition1_DR_Inhibitory_Amp_Final)), 2.*ones(1,numel(Condition2_DR_Inhibitory_Amp_Final)), 3.*ones(1,numel(Condition1_MR_Inhibitory_Amp_Final)), 4.*ones(1,numel(Condition2_MR_Inhibitory_Amp_Final))];

%% Comparison
[Condition1_DR_excitatory_cluster,Condition2_DR_excitatory_cluster,Condition1_MR_excitatory_cluster,Condition2_MR_excitatory_cluster] = SEPARATION(clusterIndices_3D_Excitatory,Excitatory_label,K_3D_Excitatory);

[Condition1_DR_inhibitory_cluster,Condition2_DR_inhibitory_cluster,Condition1_MR_inhibitory_cluster,Condition2_MR_inhibitory_cluster] = SEPARATION(clusterIndices_3D_Inhibitory,Inhibitory_label,K_3D_Inhibitory);

% 1: Condition1 DR, 2: Condition2 DR, 3: Condition1 MR, 4: Condition2 MR
%% Cell comparison

[Condition1_DR_frequency_e, Mean_Condition1_DR_cell_Cluster_e, sem_Condition1_DR_cell_Cluster_e] = cellAnalysis(Condition1_DR_Excitatory,Condition1_DR_excitatory_cluster,K_3D_Excitatory);
[Condition2_DR_frequency_e, Mean_Condition2_DR_cell_Cluster_e, sem_Condition2_DR_cell_Cluster_e] = cellAnalysis(Condition2_DR_Excitatory,Condition2_DR_excitatory_cluster,K_3D_Excitatory);
[Condition1_MR_frequency_e, Mean_Condition1_MR_cell_Cluster_e, sem_Condition1_MR_cell_Cluster_e] = cellAnalysis(Condition1_MR_Excitatory,Condition1_MR_excitatory_cluster,K_3D_Excitatory);
[Condition2_MR_frequency_e, Mean_Condition2_MR_cell_Cluster_e, sem_Condition2_MR_cell_Cluster_e] = cellAnalysis(Condition2_MR_Excitatory,Condition2_MR_excitatory_cluster,K_3D_Excitatory);

[Condition1_DR_frequency_i, Mean_Condition1_DR_cell_Cluster_i, sem_Condition1_DR_cell_Cluster_i] = cellAnalysis(Condition1_DR_Inhibitory,Condition1_DR_inhibitory_cluster,K_3D_Inhibitory);
[Condition2_DR_frequency_i, Mean_Condition2_DR_cell_Cluster_i, sem_Condition2_DR_cell_Cluster_i] = cellAnalysis(Condition2_DR_Inhibitory,Condition2_DR_inhibitory_cluster,K_3D_Inhibitory);
[Condition1_MR_frequency_i, Mean_Condition1_MR_cell_Cluster_i, sem_Condition1_MR_cell_Cluster_i] = cellAnalysis(Condition1_MR_Inhibitory,Condition1_MR_inhibitory_cluster,K_3D_Inhibitory);
[Condition2_MR_frequency_i, Mean_Condition2_MR_cell_Cluster_i, sem_Condition2_MR_cell_Cluster_i] = cellAnalysis(Condition2_MR_Inhibitory,Condition2_MR_inhibitory_cluster,K_3D_Inhibitory);
%% ttest - DR between group - Excitatory

[h7,p7] = ttest_cell(Condition1_DR_frequency_e, Condition2_DR_frequency_e, Mean_Condition1_DR_cell_Cluster_e, Mean_Condition2_DR_cell_Cluster_e, ...
    sem_Condition1_DR_cell_Cluster_e, sem_Condition2_DR_cell_Cluster_e, K_3D_Excitatory);

[h8,p8] = ttest_cell(Condition1_MR_frequency_e, Condition2_MR_frequency_e, Mean_Condition1_MR_cell_Cluster_e, Mean_Condition2_MR_cell_Cluster_e, ...
    sem_Condition1_MR_cell_Cluster_e, sem_Condition2_MR_cell_Cluster_e, K_3D_Excitatory);

[h9,p9] = ttest_cell(Condition1_DR_frequency_e, Condition1_MR_frequency_e, Mean_Condition1_DR_cell_Cluster_e, Mean_Condition1_MR_cell_Cluster_e, ...
    sem_Condition1_DR_cell_Cluster_e, sem_Condition1_MR_cell_Cluster_e, K_3D_Excitatory);

%% ttest - DR between group - Inhibitory

[h10,p10] = ttest_cell(Condition1_DR_frequency_i, Condition2_DR_frequency_i, Mean_Condition1_DR_cell_Cluster_i, Mean_Condition2_DR_cell_Cluster_i, ...
    sem_Condition1_DR_cell_Cluster_i, sem_Condition2_DR_cell_Cluster_i, K_3D_Inhibitory);

[h11,p11] = ttest_cell(Condition1_MR_frequency_i, Condition2_MR_frequency_i, Mean_Condition1_MR_cell_Cluster_i, Mean_Condition2_MR_cell_Cluster_i, ...
    sem_Condition1_MR_cell_Cluster_i, sem_Condition2_MR_cell_Cluster_i, K_3D_Inhibitory);

[h12,p12] = ttest_cell(Condition1_DR_frequency_i, Condition1_MR_frequency_i, Mean_Condition1_DR_cell_Cluster_i, Mean_Condition1_MR_cell_Cluster_i, ...
    sem_Condition1_DR_cell_Cluster_i, sem_Condition1_MR_cell_Cluster_i, K_3D_Inhibitory);


%% Excitatory
figure;
scatter3(Excitatory_merge(:,1), Excitatory_merge(:,2), Excitatory_merge(:,3),'k.');
legend('data')
xlabel('Amplitdue (pA)');
ylabel('Rise time (ms)');
zlabel('Decay (ms)');
hold on;

for i = find(cell2mat(h9)==1)
    scatter3(Excitatory_merge(clusterIndices_3D_Excitatory==i,1), Excitatory_merge(clusterIndices_3D_Excitatory==i,2), Excitatory_merge(clusterIndices_3D_Excitatory==i,3), 'filled');
end
%% Inhibitory
figure;
scatter3(Inhibitory_merge(:,1), Inhibitory_merge(:,2), Inhibitory_merge(:,3),'k.');
legend('data')
xlabel('Amplitdue (pA)');
ylabel('Rise time (ms)');
zlabel('Decay (ms)');
hold on;

for i = find(cell2mat(h12)==1)
    scatter3(Inhibitory_merge(clusterIndices_2D_Inhibitory==i,1), Inhibitory_merge(clusterIndices_2D_Inhibitory==i,2), Inhibitory_merge(clusterIndices_2D_Inhibitory==i,3), 'filled');
end
%% Local functions

function [Amp, RT, DT] = getparameters(struct)
    
    Amp = [];
    RT = [];
    DT = [];
    for ii = 1:length(struct)
        Amp = [Amp struct{ii}.Amplitude'];
        RT = [RT struct{ii}.("10-90% RT")'];
        DT = [DT struct{ii}.("decay time")'];
    end
    
end

function [Amp, RT, DT] = removeNan(orgAmp, orgRT, orgDT)
    nanidx = isnan(orgDT);
    Amp = orgAmp(~nanidx);
    RT = orgRT(~nanidx);
    DT = orgDT(~nanidx);
end


