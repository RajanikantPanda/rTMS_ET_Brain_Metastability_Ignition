clc;
clear;
close all;
%%% Load the data and make the re-formating
% CNT = load('C:\Users\bkbme\Desktop\Sneha_program\PMT\data1\CTL_data3d_aal.mat');
% CNT=CNT.data_3d; CNT = num2cell(CNT,[2,3])';
% for i=1:size(CNT,2)
%     CNT{i} = squeeze(CNT{1,i});
% end
% %MCS = load('C:\Users\bkbme\Desktop\Sneha_program\PMT\ODN.mat');
% MCS = load('C:\Users\bkbme\Desktop\Sneha_program\PMT\data1\ODN_data3d_aal.mat');
% MCS=MCS.data_3d; MCS = num2cell(MCS,[2,3])';
% for i=1:size(MCS,2)
%     MCS{i} = squeeze(MCS{1,i});
% end
% %UWS = load('C:\Users\bkbme\Desktop\Sneha_program\PMT\ODP.mat');
% UWS = load('C:\Users\bkbme\Desktop\Sneha_program\PMT\data1\ODP_data3d_aal.mat');
% UWS=UWS.data_3d; UWS = num2cell(UWS,[2,3])';
% for i=1:size(UWS,2)
%     UWS{i} = squeeze(UWS{1,i});
% end
load('C:\Research\RAJ_RP_BlueHD\NMS\ET\fMRI\Dyn\ET_HC_ShenArt.mat')
CNT = ET_HC;
MCS = ET_Pre;
UWS = ET_Post;
%% to Compute Syncronisation and DFC matrices
[CNT_meta, CNT_sync, CNT_dFC, CNT_dFC_cos, CNT_int, CNT_seg] = dynamicMetrics_2(CNT,2.5);  
[MCS_meta, MCS_sync, MCS_dFC, MCS_dFC_cos, MCS_int, MCS_seg] = dynamicMetrics_2(MCS,2.5);
[UWS_meta, UWS_sync, UWS_dFC, UWS_dFC_cos, UWS_int, UWS_seg] = dynamicMetrics_2(UWS,2.5);
%% Ploting metastability Whole Brain
Grp_con = ones(length(CNT),1); Grp_mcs = 2*(ones(length(MCS),1)); Grp_uws = 3*(ones(length(UWS),1));
Grp = [Grp_con; Grp_mcs; Grp_uws]; %Grp = Grp';
Grp_meta = [CNT_meta'; MCS_meta'; UWS_meta'];
figure(1); notBoxPlot(Grp_meta,Grp,0.5,'patch',ones(length(Grp_meta),1));
%% Ploting wholebrain intigration Whol Brain
CNT_wb_int = mean(CNT_int,2);MCS_wb_int = mean(MCS_int,2);UWS_wb_int = mean(UWS_int,2);
Grp_integ = [CNT_wb_int; MCS_wb_int; UWS_wb_int];
figure (2); notBoxPlot(Grp_integ,Grp,0.5,'patch',ones(length(Grp_integ),1));
%% Ploting wholebrain Segrigation Whol Brain
CNT_wb_seg = mean(CNT_seg,2);MCS_wb_seg = mean(MCS_seg,2);UWS_wb_seg = mean(UWS_seg,2);
Grp_seg = [CNT_wb_seg; MCS_wb_seg; UWS_wb_seg];
figure (3); notBoxPlot(Grp_seg,Grp,0.5,'patch',ones(length(Grp_seg),1));
%% %% Plot the mean connectivity
for i=1:size(CNT_dFC,2)
%CNT_dFC_all(i,:,:) = mean(abs(CNT_dFC {1,i}),3);
CNT_dFC_all(i,:,:) = mean((CNT_dFC_cos {1,i}),3);
end
for i=1:size(MCS_dFC,2)
%MCS_dFC_all(i,:,:) = mean(abs(MCS_dFC {1,i}),3);
MCS_dFC_all(i,:,:) = mean((MCS_dFC_cos {1,i}),3);
end
for i=1:size(UWS_dFC,2)
%UWS_dFC_all(i,:,:) = mean(abs(UWS_dFC {1,i}),3);
UWS_dFC_all(i,:,:) = mean((UWS_dFC_cos {1,i}),3);
end
close all
figure; imagesc(squeeze(mean(CNT_dFC_all,1)));colorbar;caxis([-0.2 0.5])
title ('Mean DFC of Healthy Control')
xlabel('No of ROIs', 'FontName','Calibri', 'FontSize', 14);
ylabel('No of ROIs', 'FontName','Calibri', 'FontSize', 14);
set(gca,'FontSize',14, 'FontName','Calibri', 'FontWeight', 'bold')
figure; imagesc(squeeze(mean(MCS_dFC_all,1)));colorbar;caxis([-0.2 0.5])
title ('Mean DFC of PD-Normal Cognition')
xlabel('No of ROIs', 'FontName','Calibri', 'FontSize', 14);
ylabel('No of ROIs', 'FontName','Calibri', 'FontSize', 14);
set(gca,'FontSize',14, 'FontName','Calibri', 'FontWeight', 'bold')
figure; imagesc(squeeze(mean(UWS_dFC_all,1)));colorbar;caxis([-0.2 0.5])
title ('Mean DFC of PD-Severe Hyposmia')
xlabel('No of ROIs', 'FontName','Calibri', 'FontSize', 11);
ylabel('No of ROIs', 'FontName','Calibri', 'FontSize', 11);
set(gca,'FontSize',14, 'FontName','Calibri', 'FontWeight', 'bold')
CNT_dFC_all_wb=mean(mean(CNT_dFC_all,3),2);
MCS_dFC_all_wb=mean(mean(MCS_dFC_all,3),2);
UWS_dFC_all_wb=mean(mean(UWS_dFC_all,3),2);
Grp_dFC_wb = [CNT_dFC_all_wb; MCS_dFC_all_wb; UWS_dFC_all_wb];
Grp_con = ones(length(CNT_dFC_all_wb),1); Grp_mcs = 2*(ones(length(MCS_dFC_all_wb),1)); Grp_uws = 3*(ones(length(UWS_dFC_all_wb),1));
Grp = [Grp_con; Grp_mcs; Grp_uws];
figure; notBoxPlot(Grp_dFC_wb,Grp,0.5,'patch',ones(length(Grp_dFC_wb),1)); title ('Whole brain mean dFC')
xticklabels({'HC','PD NC','PD SH'}); set(gca,'FontSize',14, 'FontName','Times New Roman', 'FontWeight', 'bold')
%% select, reshape and concatinate data

all_dFC = [UWS_dFC_cos MCS_dFC_cos CNT_dFC_cos];                           % concatinate data for subject groups
%clear UWS_dFC_cos MCS_dFC_cos CNT_dFC_cos
%% Data Arrangement
n_subjects = length(all_dFC);                                              % total number of subjects
[n_regions, ~, n_images] = size(all_dFC{1});                               % dimensions of subject data
n_features = n_regions*(n_regions - 1)/2;                                  % number of pairs between ROI's

mask = tril(true(n_regions), -1);                                          % mask for selection of upper triangle (=lower triangle because of symmetry)
mask = repmat(mask, [1, 1, n_images]);                                     % extand mask for number of images

data = zeros([n_features, n_subjects*n_images]);                           % empty matrix for concatinated data
for idx = 1:n_subjects                                                     % loop over subjects
    dFC = all_dFC{idx};                                                    % select subject data
    dFC = reshape(dFC(mask), [n_features, n_images]);                      % apply mask to extract lower triangle and reshape to n_features x n_images
    data(:, (idx-1)*n_images+1:idx*n_images) = dFC;                        
end
data = data.';                                                             % observations x features

%% k-means clustering

k = 4;                                                                     % number of desired clusters/patterns

[labels, clusters] = kmeans(data, k, ...                                   % perform k-means clustering
    'Distance', 'cityblock', ...                                           % Manhattan distance
    'Replicates', 50);                                                      % number of times to repeat algorithm with random initialisation
%% Load the results
clc
clear
load('C:\Users\bkbme\Desktop\Sneha_program\PMT\PhaseBased_DFC\BrainStates_4_50_final.mat')
%% reconstruct phase-difference matrices
%
mask = tril(true(n_regions), -1);                                          % mask for selection of upper triangle

patterns = cell(1, 4);
for idx = 1:k
    pattern = zeros(n_regions);
    pattern(mask) = clusters(idx, :);
    patterns{idx} = pattern + pattern.' + eye(n_regions);                  % reconstruct symmetric matrix
end

%%% labels

labels = reshape(labels, [n_images, n_subjects]);                          % label of pattern corresponding to each image, for each subject (columns)

%% probabilities

prob = zeros(k, n_subjects);
for idx = 1:n_subjects
    prob_i = labels(:, idx) == 1:k;
    prob(:, idx) = sum(prob_i, 1) / n_images;
end
%
prob=[prob(:,1:4) prob(:,21) prob(:,6:20) prob(:,5) prob(:,22:42)]

%% plot Brain State Matrix results
% % figure
% % for idx = 1:k
% %     subplot(1, k, idx)
% %     imagesc(patterns{idx})
% % end

%% plot Brain State Matrix results
figure
subplot(1, 4, 1)
imagesc(patterns{2})
title ('Pattern A (Complex)')
xlabel('No of ROIs', 'FontName','Calibri', 'FontSize', 14);
ylabel('No of ROIs', 'FontName','Calibri', 'FontSize', 14);
set(gca,'FontSize',14, 'FontName','Calibri', 'FontWeight', 'bold')
subplot(1, 4, 2)
imagesc(patterns{4})
title ('Pattern B')
xlabel('No of ROIs', 'FontName','Calibri', 'FontSize', 14);
ylabel('No of ROIs', 'FontName','Calibri', 'FontSize', 14);
set(gca,'FontSize',14, 'FontName','Calibri', 'FontWeight', 'bold')
subplot(1, 4, 3)
imagesc(patterns{1})
title ('Pattern C')
xlabel('No of ROIs', 'FontName','Calibri', 'FontSize', 14);
ylabel('No of ROIs', 'FontName','Calibri', 'FontSize', 14);
set(gca,'FontSize',14, 'FontName','Calibri', 'FontWeight', 'bold')
subplot(1, 4, 4)
imagesc(patterns{3})
title ('Pattern D (Simple)')
xlabel('No of ROIs', 'FontName','Calibri', 'FontSize', 14);
ylabel('No of ROIs', 'FontName','Calibri', 'FontSize', 14);
set(gca,'FontSize',14, 'FontName','Calibri', 'FontWeight', 'bold')
%%
Grp_HC = ones(12,1); Grp_MCS = 2*ones(15,1); Grp_UWS = 3*ones(15,1); 
Grp = [Grp_UWS; Grp_MCS; Grp_HC]; %Grp = Gzarp';
figure; notBoxPlot(prob(2,:)',Grp,0.5,'patch',ones(length(Grp),1)); title('Pattern A (Complex)')
ylabel('Rate of Pattern Occurance (Probability)', 'FontName','Calibri', 'FontSize', 14);
xticklabels({'HC','PD-NC','PD-SH'});
set(gca,'FontSize',14, 'FontName','Calibri', 'FontWeight', 'bold')
figure; notBoxPlot(prob(4,:)',Grp,0.5,'patch',ones(length(Grp),1)); title('Pattern B')
ylabel('Rate of Pattern Occurance (Probability)', 'FontName','Calibri', 'FontSize', 14);
xticklabels({'HC','PD-NC','PD-SH'});
set(gca,'FontSize',14, 'FontName','Calibri', 'FontWeight', 'bold')
%%%
figure; notBoxPlot(prob(1,:)',Grp,0.5,'patch',ones(length(Grp),1)); title('Pattern C')
ylabel('Rate of Pattern Occurance (Probability)', 'FontName','Calibri', 'FontSize', 14);
xticklabels({'HC','PD-NC','PD-SH'});
set(gca,'FontSize',14, 'FontName','Calibri', 'FontWeight', 'bold')
figure; notBoxPlot(prob(3,:)',Grp,0.5,'patch',ones(length(Grp),1)); title('Pattern D (Simple)')
ylabel('Rate of Pattern Occurance (Probability)', 'FontName','Calibri', 'FontSize', 14);
xticklabels({'HC','PD-NC','PD-SH'});
set(gca,'FontSize',14, 'FontName','Calibri', 'FontWeight', 'bold')
%figure; notBoxPlot((prob(1,:)./prob(3,:))',Grp,0.5,'patch',ones(length(Grp),1));
%% save results
% save('/Users/michielmeys/Desktop/patterns.mat', ...
%     'patterns', 'labels', 'prob')                                          % save results
save BrainStates_4_50__final