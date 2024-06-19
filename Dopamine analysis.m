clear all
clc
close all
load('Dopamine_neuron.mat')
%% Condition 1

condition1 = fields(finalStruct.Dopamine_neuron_condition_1);
Condition1_SNc_Excitatory = cell(0);
Condition1_VTA_Excitatory = cell(0);
Condition1_SNc_Inhibitory = cell(0);
Condition1_VTA_Inhibitory = cell(0);

for ii = 1:length(fields(finalStruct.Dopamine_neuron_condition_1))
    Excitatory = fields(finalStruct.Dopamine_neuron_condition_1.(condition1{ii}).Excitatory);
    ExcitatorySNc = Excitatory(contains(Excitatory,'SNc'));
    ExcitatoryVTA = Excitatory(contains(Excitatory,'VTA'));

    Inhibitory = fields(finalStruct.Dopamine_neuron_condition_1.(condition1{ii}).Inhibitory);
    InhibitorySNc = Inhibitory(contains(Inhibitory,'SNc'));
    InhibitoryVTA = Inhibitory(contains(Inhibitory,'VTA'));

    for jj = 1:length(ExcitatorySNc)
        Condition1_SNc_Excitatory{end+1,1} = finalStruct.Dopamine_neuron_condition_1.(condition1{ii}).Excitatory.(ExcitatorySNc{jj});
    end
    for jj = 1:length(ExcitatoryVTA)
        Condition1_VTA_Excitatory{end+1,1} = finalStruct.Dopamine_neuron_condition_1.(condition1{ii}).Excitatory.(ExcitatoryVTA{jj});
    end
    for jj = 1:length(InhibitorySNc)
        Condition1_SNc_Inhibitory{end+1,1} = finalStruct.Dopamine_neuron_condition_1.(condition1{ii}).Inhibitory.(InhibitorySNc{jj});
    end
    for jj = 1:length(InhibitoryVTA)
        Condition1_VTA_Inhibitory{end+1,1} = finalStruct.Dopamine_neuron_condition_1.(condition1{ii}).Inhibitory.(InhibitoryVTA{jj});
    end
end
%% Condition1 - SNc

[Condition1_SNc_Excitatory_Amp, Condition1_SNc_Excitatory_RT, Condition1_SNc_Excitatory_DT] = getparameters(Condition1_SNc_Excitatory);

[Condition1_SNc_Inhibitory_Amp, Condition1_SNc_Inhibitory_RT, Condition1_SNc_Inhibitory_DT] = getparameters(Condition1_SNc_Inhibitory);

% %% 2D plot - Condition1 - SNc - Amp, RT
% fig1 = figure;
% plot(Condition1_SNc_Excitatory_Amp,Condition1_SNc_Excitatory_RT,'LineStyle','none','Marker','+','MarkerSize',1);
% 
% fig2 = figure;
% plot(Condition1_SNc_Inhibitory_Amp,Condition1_SNc_Inhibitory_RT,'LineStyle','none','Marker','+','MarkerSize',1);

%% Condition1 - VTA

[Condition1_VTA_Excitatory_Amp, Condition1_VTA_Excitatory_RT, Condition1_VTA_Excitatory_DT] = getparameters(Condition1_VTA_Excitatory);

[Condition1_VTA_Inhibitory_Amp, Condition1_VTA_Inhibitory_RT, Condition1_VTA_Inhibitory_DT] = getparameters(Condition1_VTA_Inhibitory);

% %% 2D plot - Amp, RT
% fig3 = figure;
% plot(Condition1_VTA_Excitatory_Amp,Condition1_VTA_Excitatory_RT,'LineStyle','none','Marker','+','MarkerSize',1);
% 
% fig4 = figure;
% plot(Condition1_VTA_Inhibitory_Amp,Condition1_VTA_Inhibitory_RT,'LineStyle','none','Marker','+','MarkerSize',1);

%% Condition 2

condition2 = fields(finalStruct.Dopamine_neuron_condition_2);
Condition2_SNc_Excitatory = cell(0);
Condition2_VTA_Excitatory = cell(0);
Condition2_SNc_Inhibitory = cell(0);
Condition2_VTA_Inhibitory = cell(0);

for ii = 1:length(fields(finalStruct.Dopamine_neuron_condition_2))
    Excitatory = fields(finalStruct.Dopamine_neuron_condition_2.(condition2{ii}).Excitatory);
    ExcitatorySNc = Excitatory(contains(Excitatory,'SNc'));
    ExcitatoryVTA = Excitatory(contains(Excitatory,'VTA'));

    Inhibitory = fields(finalStruct.Dopamine_neuron_condition_2.(condition2{ii}).Inhibitory);
    InhibitorySNc = Inhibitory(contains(Inhibitory,'SNc'));
    InhibitoryVTA = Inhibitory(contains(Inhibitory,'VTA'));

    for jj = 1:length(ExcitatorySNc)
        Condition2_SNc_Excitatory{end+1,1} = finalStruct.Dopamine_neuron_condition_2.(condition2{ii}).Excitatory.(ExcitatorySNc{jj});
    end
    for jj = 1:length(ExcitatoryVTA)
        Condition2_VTA_Excitatory{end+1,1} = finalStruct.Dopamine_neuron_condition_2.(condition2{ii}).Excitatory.(ExcitatoryVTA{jj});
    end
    for jj = 1:length(InhibitorySNc)
        Condition2_SNc_Inhibitory{end+1,1} = finalStruct.Dopamine_neuron_condition_2.(condition2{ii}).Inhibitory.(InhibitorySNc{jj});
    end
    for jj = 1:length(InhibitoryVTA)
        Condition2_VTA_Inhibitory{end+1,1} = finalStruct.Dopamine_neuron_condition_2.(condition2{ii}).Inhibitory.(InhibitoryVTA{jj});
    end
end

%% Condition2 - SNc

[Condition2_SNc_Excitatory_Amp, Condition2_SNc_Excitatory_RT, Condition2_SNc_Excitatory_DT] = getparameters(Condition2_SNc_Excitatory);

[Condition2_SNc_Inhibitory_Amp, Condition2_SNc_Inhibitory_RT, Condition2_SNc_Inhibitory_DT] = getparameters(Condition2_SNc_Inhibitory);

% %% 2D plot - Condition2 - SNc - Amp, RT
% fig5 = figure;
% plot(Condition2_SNc_Excitatory_Amp,Condition2_SNc_Excitatory_RT,'LineStyle','none','Marker','+','MarkerSize',1);
% 
% fig6 = figure;
% plot(Condition2_SNc_Inhibitory_Amp,Condition2_SNc_Inhibitory_RT,'LineStyle','none','Marker','+','MarkerSize',1);

%% Condition2 - VTA

[Condition2_VTA_Excitatory_Amp, Condition2_VTA_Excitatory_RT, Condition2_VTA_Excitatory_DT] = getparameters(Condition2_VTA_Excitatory);

[Condition2_VTA_Inhibitory_Amp, Condition2_VTA_Inhibitory_RT, Condition2_VTA_Inhibitory_DT] = getparameters(Condition2_VTA_Inhibitory);

% %% 2D plot - Condition2 - VTA - Amp, RT
% fig7 = figure;
% plot(Condition2_VTA_Excitatory_Amp,Condition2_VTA_Excitatory_RT,'LineStyle','none','Marker','+','MarkerSize',1);
% 
% fig8 = figure;
% plot(Condition2_VTA_Inhibitory_Amp,Condition2_VTA_Inhibitory_RT,'LineStyle','none','Marker','+','MarkerSize',1);

%% Remove Nan

[Condition1_SNc_Excitatory_Amp_Final, Condition1_SNc_Excitatory_RT_Final, Condition1_SNc_Excitatory_DT_Final] = removeNan( ...
    Condition1_SNc_Excitatory_Amp, Condition1_SNc_Excitatory_RT, Condition1_SNc_Excitatory_DT);

[Condition1_VTA_Excitatory_Amp_Final, Condition1_VTA_Excitatory_RT_Final, Condition1_VTA_Excitatory_DT_Final] = removeNan( ...
    Condition1_VTA_Excitatory_Amp, Condition1_VTA_Excitatory_RT, Condition1_VTA_Excitatory_DT);

[Condition1_SNc_Inhibitory_Amp_Final, Condition1_SNc_Inhibitory_RT_Final, Condition1_SNc_Inhibitory_DT_Final] = removeNan( ...
    Condition1_SNc_Inhibitory_Amp, Condition1_SNc_Inhibitory_RT, Condition1_SNc_Inhibitory_DT);

[Condition1_VTA_Inhibitory_Amp_Final, Condition1_VTA_Inhibitory_RT_Final, Condition1_VTA_Inhibitory_DT_Final] = removeNan( ...
    Condition1_VTA_Inhibitory_Amp, Condition1_VTA_Inhibitory_RT, Condition1_VTA_Inhibitory_DT);


[Condition2_SNc_Excitatory_Amp_Final, Condition2_SNc_Excitatory_RT_Final, Condition2_SNc_Excitatory_DT_Final] = removeNan( ...
    Condition2_SNc_Excitatory_Amp, Condition2_SNc_Excitatory_RT, Condition2_SNc_Excitatory_DT);

[Condition2_VTA_Excitatory_Amp_Final, Condition2_VTA_Excitatory_RT_Final, Condition2_VTA_Excitatory_DT_Final] = removeNan( ...
    Condition2_VTA_Excitatory_Amp, Condition2_VTA_Excitatory_RT, Condition2_VTA_Excitatory_DT);

[Condition2_SNc_Inhibitory_Amp_Final, Condition2_SNc_Inhibitory_RT_Final, Condition2_SNc_Inhibitory_DT_Final] = removeNan( ...
    Condition2_SNc_Inhibitory_Amp, Condition2_SNc_Inhibitory_RT, Condition2_SNc_Inhibitory_DT);

[Condition2_VTA_Inhibitory_Amp_Final, Condition2_VTA_Inhibitory_RT_Final, Condition2_VTA_Inhibitory_DT_Final] = removeNan( ...
    Condition2_VTA_Inhibitory_Amp, Condition2_VTA_Inhibitory_RT, Condition2_VTA_Inhibitory_DT);


%% Merge data

Excitatory_Amp = [Condition1_SNc_Excitatory_Amp_Final, Condition2_SNc_Excitatory_Amp_Final, Condition1_VTA_Excitatory_Amp_Final, Condition2_VTA_Excitatory_Amp_Final];
Excitatory_RT = [Condition1_SNc_Excitatory_RT_Final, Condition2_SNc_Excitatory_RT_Final, Condition1_VTA_Excitatory_RT_Final, Condition2_VTA_Excitatory_RT_Final];
Excitatory_DT = [Condition1_SNc_Excitatory_DT_Final, Condition2_SNc_Excitatory_DT_Final, Condition1_VTA_Excitatory_DT_Final, Condition2_VTA_Excitatory_DT_Final];
Excitatory_merge = [Excitatory_Amp' Excitatory_RT' Excitatory_DT']; % Amplitude | Risetime | Decaytime


Inhibitory_Amp = [Condition1_SNc_Inhibitory_Amp_Final, Condition2_SNc_Inhibitory_Amp_Final, Condition1_VTA_Inhibitory_Amp_Final, Condition2_VTA_Inhibitory_Amp_Final];
Inhibitory_RT = [Condition1_SNc_Inhibitory_RT_Final, Condition2_SNc_Inhibitory_RT_Final, Condition1_VTA_Inhibitory_RT_Final, Condition2_VTA_Inhibitory_RT_Final];
Inhibitory_DT = [Condition1_SNc_Inhibitory_DT_Final, Condition2_SNc_Inhibitory_DT_Final, Condition1_VTA_Inhibitory_DT_Final, Condition2_VTA_Inhibitory_DT_Final];
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

Excitatory_label = [ones(1,numel(Condition1_SNc_Excitatory_Amp_Final)), 2.*ones(1,numel(Condition2_SNc_Excitatory_Amp_Final)), 3.*ones(1,numel(Condition1_VTA_Excitatory_Amp_Final)), 4.*ones(1,numel(Condition2_VTA_Excitatory_Amp_Final))];
Inhibitory_label = [ones(1,numel(Condition1_SNc_Inhibitory_Amp_Final)), 2.*ones(1,numel(Condition2_SNc_Inhibitory_Amp_Final)), 3.*ones(1,numel(Condition1_VTA_Inhibitory_Amp_Final)), 4.*ones(1,numel(Condition2_VTA_Inhibitory_Amp_Final))];

%% Comparison
[Condition1_SNc_excitatory_cluster,Condition2_SNc_excitatory_cluster,Condition1_VTA_excitatory_cluster,Condition2_VTA_excitatory_cluster] = SEPARATION(clusterIndices_2D_Excitatory,Excitatory_label,K_2D_Excitatory);

[Condition1_SNc_inhibitory_cluster,Condition2_SNc_inhibitory_cluster,Condition1_VTA_inhibitory_cluster,Condition2_VTA_inhibitory_cluster] = SEPARATION(clusterIndices_2D_Inhibitory,Inhibitory_label,K_2D_Inhibitory);

% 1: Condition1 SNc, 2: Condition2 SNc, 3: Condition1 VTA, 4: Condition2 VTA

%% Cell comparison

[Condition1_SNc_frequency_e, Mean_Condition1_SNc_cell_Cluster_e, sem_Condition1_SNc_cell_Cluster_e] = cellAnalysis(Condition1_SNc_Excitatory,Condition1_SNc_excitatory_cluster,K_2D_Excitatory);
[Condition2_SNc_frequency_e, Mean_Condition2_SNc_cell_Cluster_e, sem_Condition2_SNc_cell_Cluster_e] = cellAnalysis(Condition2_SNc_Excitatory,Condition2_SNc_excitatory_cluster,K_2D_Excitatory);
[Condition1_VTA_frequency_e, Mean_Condition1_VTA_cell_Cluster_e, sem_Condition1_VTA_cell_Cluster_e] = cellAnalysis(Condition1_VTA_Excitatory,Condition1_VTA_excitatory_cluster,K_2D_Excitatory);
[Condition2_VTA_frequency_e, Mean_Condition2_VTA_cell_Cluster_e, sem_Condition2_VTA_cell_Cluster_e] = cellAnalysis(Condition2_VTA_Excitatory,Condition2_VTA_excitatory_cluster,K_2D_Excitatory);

[Condition1_SNc_frequency_i, Mean_Condition1_SNc_cell_Cluster_i, sem_Condition1_SNc_cell_Cluster_i] = cellAnalysis(Condition1_SNc_Inhibitory,Condition1_SNc_inhibitory_cluster,K_2D_Inhibitory);
[Condition2_SNc_frequency_i, Mean_Condition2_SNc_cell_Cluster_i, sem_Condition2_SNc_cell_Cluster_i] = cellAnalysis(Condition2_SNc_Inhibitory,Condition2_SNc_inhibitory_cluster,K_2D_Inhibitory);
[Condition1_VTA_frequency_i, Mean_Condition1_VTA_cell_Cluster_i, sem_Condition1_VTA_cell_Cluster_i] = cellAnalysis(Condition1_VTA_Inhibitory,Condition1_VTA_inhibitory_cluster,K_2D_Inhibitory);
[Condition2_VTA_frequency_i, Mean_Condition2_VTA_cell_Cluster_i, sem_Condition2_VTA_cell_Cluster_i] = cellAnalysis(Condition2_VTA_Inhibitory,Condition2_VTA_inhibitory_cluster,K_2D_Inhibitory);
%% ttest - SNc between group - Excitatory

[h1,p1] = ttest_cell(Condition1_SNc_frequency_e, Condition2_SNc_frequency_e, Mean_Condition1_SNc_cell_Cluster_e, Mean_Condition2_SNc_cell_Cluster_e, ...
    sem_Condition1_SNc_cell_Cluster_e, sem_Condition2_SNc_cell_Cluster_e, K_2D_Excitatory);

[h2,p2] = ttest_cell(Condition1_VTA_frequency_e, Condition2_VTA_frequency_e, Mean_Condition1_VTA_cell_Cluster_e, Mean_Condition2_VTA_cell_Cluster_e, ...
    sem_Condition1_VTA_cell_Cluster_e, sem_Condition2_VTA_cell_Cluster_e, K_2D_Excitatory);

[h3,p3] = ttest_cell(Condition1_SNc_frequency_e, Condition1_VTA_frequency_e, Mean_Condition1_SNc_cell_Cluster_e, Mean_Condition1_VTA_cell_Cluster_e, ...
    sem_Condition1_SNc_cell_Cluster_e, sem_Condition1_VTA_cell_Cluster_e, K_2D_Excitatory);

%% ttest - SNc between group - Inhibitory

[h4,p4] = ttest_cell(Condition1_SNc_frequency_i, Condition2_SNc_frequency_i, Mean_Condition1_SNc_cell_Cluster_i, Mean_Condition2_SNc_cell_Cluster_i, ...
    sem_Condition1_SNc_cell_Cluster_i, sem_Condition2_SNc_cell_Cluster_i, K_2D_Inhibitory);

[h5,p5] = ttest_cell(Condition1_VTA_frequency_i, Condition2_VTA_frequency_i, Mean_Condition1_VTA_cell_Cluster_i, Mean_Condition2_VTA_cell_Cluster_i, ...
    sem_Condition1_VTA_cell_Cluster_i, sem_Condition2_VTA_cell_Cluster_i, K_2D_Inhibitory);

[h6,p6] = ttest_cell(Condition1_SNc_frequency_i, Condition1_VTA_frequency_i, Mean_Condition1_SNc_cell_Cluster_i, Mean_Condition1_VTA_cell_Cluster_i, ...
    sem_Condition1_SNc_cell_Cluster_i, sem_Condition1_VTA_cell_Cluster_i, K_2D_Inhibitory);

%% Excitatory
figure;
scatter(transformedData_Excitatory(:,1), transformedData_Excitatory(:,2),'k.');
xlabel('x1');
ylabel('x2');
hold on;

for i = find(cell2mat(h3)==1)
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

Excitatory_label = [ones(1,numel(Condition1_SNc_Excitatory_Amp_Final)), 2.*ones(1,numel(Condition2_SNc_Excitatory_Amp_Final)), 3.*ones(1,numel(Condition1_VTA_Excitatory_Amp_Final)), 4.*ones(1,numel(Condition2_VTA_Excitatory_Amp_Final))];
Inhibitory_label = [ones(1,numel(Condition1_SNc_Inhibitory_Amp_Final)), 2.*ones(1,numel(Condition2_SNc_Inhibitory_Amp_Final)), 3.*ones(1,numel(Condition1_VTA_Inhibitory_Amp_Final)), 4.*ones(1,numel(Condition2_VTA_Inhibitory_Amp_Final))];

%% Comparison
[Condition1_SNc_excitatory_cluster,Condition2_SNc_excitatory_cluster,Condition1_VTA_excitatory_cluster,Condition2_VTA_excitatory_cluster] = SEPARATION(clusterIndices_3D_Excitatory,Excitatory_label,K_3D_Excitatory);

[Condition1_SNc_inhibitory_cluster,Condition2_SNc_inhibitory_cluster,Condition1_VTA_inhibitory_cluster,Condition2_VTA_inhibitory_cluster] = SEPARATION(clusterIndices_3D_Inhibitory,Inhibitory_label,K_3D_Inhibitory);

% 1: Condition1 SNc, 2: Condition2 SNc, 3: Condition1 VTA, 4: Condition2 VTA
%% Cell comparison

[Condition1_SNc_frequency_e, Mean_Condition1_SNc_cell_Cluster_e, sem_Condition1_SNc_cell_Cluster_e] = cellAnalysis(Condition1_SNc_Excitatory,Condition1_SNc_excitatory_cluster,K_3D_Excitatory);
[Condition2_SNc_frequency_e, Mean_Condition2_SNc_cell_Cluster_e, sem_Condition2_SNc_cell_Cluster_e] = cellAnalysis(Condition2_SNc_Excitatory,Condition2_SNc_excitatory_cluster,K_3D_Excitatory);
[Condition1_VTA_frequency_e, Mean_Condition1_VTA_cell_Cluster_e, sem_Condition1_VTA_cell_Cluster_e] = cellAnalysis(Condition1_VTA_Excitatory,Condition1_VTA_excitatory_cluster,K_3D_Excitatory);
[Condition2_VTA_frequency_e, Mean_Condition2_VTA_cell_Cluster_e, sem_Condition2_VTA_cell_Cluster_e] = cellAnalysis(Condition2_VTA_Excitatory,Condition2_VTA_excitatory_cluster,K_3D_Excitatory);

[Condition1_SNc_frequency_i, Mean_Condition1_SNc_cell_Cluster_i, sem_Condition1_SNc_cell_Cluster_i] = cellAnalysis(Condition1_SNc_Inhibitory,Condition1_SNc_inhibitory_cluster,K_3D_Inhibitory);
[Condition2_SNc_frequency_i, Mean_Condition2_SNc_cell_Cluster_i, sem_Condition2_SNc_cell_Cluster_i] = cellAnalysis(Condition2_SNc_Inhibitory,Condition2_SNc_inhibitory_cluster,K_3D_Inhibitory);
[Condition1_VTA_frequency_i, Mean_Condition1_VTA_cell_Cluster_i, sem_Condition1_VTA_cell_Cluster_i] = cellAnalysis(Condition1_VTA_Inhibitory,Condition1_VTA_inhibitory_cluster,K_3D_Inhibitory);
[Condition2_VTA_frequency_i, Mean_Condition2_VTA_cell_Cluster_i, sem_Condition2_VTA_cell_Cluster_i] = cellAnalysis(Condition2_VTA_Inhibitory,Condition2_VTA_inhibitory_cluster,K_3D_Inhibitory);
%% ttest - SNc between group - Excitatory

[h7,p7] = ttest_cell(Condition1_SNc_frequency_e, Condition2_SNc_frequency_e, Mean_Condition1_SNc_cell_Cluster_e, Mean_Condition2_SNc_cell_Cluster_e, ...
    sem_Condition1_SNc_cell_Cluster_e, sem_Condition2_SNc_cell_Cluster_e, K_3D_Excitatory);

[h8,p8] = ttest_cell(Condition1_VTA_frequency_e, Condition2_VTA_frequency_e, Mean_Condition1_VTA_cell_Cluster_e, Mean_Condition2_VTA_cell_Cluster_e, ...
    sem_Condition1_VTA_cell_Cluster_e, sem_Condition2_VTA_cell_Cluster_e, K_3D_Excitatory);

[h9,p9] = ttest_cell(Condition1_SNc_frequency_e, Condition1_VTA_frequency_e, Mean_Condition1_SNc_cell_Cluster_e, Mean_Condition1_VTA_cell_Cluster_e, ...
    sem_Condition1_SNc_cell_Cluster_e, sem_Condition1_VTA_cell_Cluster_e, K_3D_Excitatory);

%% ttest - SNc between group - Inhibitory

[h10,p10] = ttest_cell(Condition1_SNc_frequency_i, Condition2_SNc_frequency_i, Mean_Condition1_SNc_cell_Cluster_i, Mean_Condition2_SNc_cell_Cluster_i, ...
    sem_Condition1_SNc_cell_Cluster_i, sem_Condition2_SNc_cell_Cluster_i, K_3D_Inhibitory);

[h11,p11] = ttest_cell(Condition1_VTA_frequency_i, Condition2_VTA_frequency_i, Mean_Condition1_VTA_cell_Cluster_i, Mean_Condition2_VTA_cell_Cluster_i, ...
    sem_Condition1_VTA_cell_Cluster_i, sem_Condition2_VTA_cell_Cluster_i, K_3D_Inhibitory);

[h12,p12] = ttest_cell(Condition1_SNc_frequency_i, Condition1_VTA_frequency_i, Mean_Condition1_SNc_cell_Cluster_i, Mean_Condition1_VTA_cell_Cluster_i, ...
    sem_Condition1_SNc_cell_Cluster_i, sem_Condition1_VTA_cell_Cluster_i, K_3D_Inhibitory);


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


