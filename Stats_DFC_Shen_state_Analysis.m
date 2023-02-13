% To load the data, call the DFC function and compute all DFC measures, Plot the results and to do the starts

clc;
clear;
close all;
load('C:\Research\RAJ_RP_BlueHD\NMS\ET\fMRI\Dyn\ET_HC_ShenArt.mat')
CNT = ET_HC;
ET_Pre = ET_Pre;
ET_Post = ET_Post;
%% to Compute Syncronisation and DFC matrices
[CNT_meta, CNT_sync, CNT_dFC, CNT_dFC_cos, CNT_int, CNT_seg, CNT_CDC,CNT_matrix_CDC] = dynamicMetrics_3(CNT,2.5);
[ET_Pre_meta, ET_Pre_sync, ET_Pre_dFC, ET_Pre_dFC_cos, ET_Pre_int, ET_Pre_seg, ET_Pre_CDC,ET_Pre_matrix_CDC] = dynamicMetrics_3(ET_Pre,2.5);
[ET_Post_meta, ET_Post_sync, ET_Post_dFC, ET_Post_dFC_cos, ET_Post_int, ET_Post_seg, ET_Post_CDC,ET_Post_matrix_CDC] = dynamicMetrics_3(ET_Post,2.5);
%%
[CNT_meta, CNT_sync, CNT_dFC, CNT_dFC_cos, CNT_int, CNT_seg] = dynamicMetrics_2(CNT,2.5);  
[ET_Pre_meta, ET_Pre_sync, ET_Pre_dFC, ET_Pre_dFC_cos, ET_Pre_int, ET_Pre_seg] = dynamicMetrics_2(ET_Pre,2.5);
[ET_Post_meta, ET_Post_sync, ET_Post_dFC, ET_Post_dFC_cos, ET_Post_int, ET_Post_seg] = dynamicMetrics_2(ET_Post,2.5);
%% Ploting metastability Whole Brain
Grp_con = ones(length(CNT),1); Grp_ET_Pre = 2*(ones(length(ET_Pre),1)); Grp_ET_Post = 3*(ones(length(ET_Post),1));
Grp = [Grp_con; Grp_ET_Pre; Grp_ET_Post]; %Grp = Grp';
Grp_meta = [CNT_meta'; ET_Pre_meta'; ET_Post_meta'];
figure(1); notBoxPlot(Grp_meta,Grp,0.5,'patch',ones(length(Grp_meta),1));
%% Ploting wholebrain intigration Whol Brain
CNT_wb_int = mean(CNT_int,2);ET_Pre_wb_int = mean(ET_Pre_int,2);ET_Post_wb_int = mean(ET_Post_int,2);
Grp_integ = [CNT_wb_int; ET_Pre_wb_int; ET_Post_wb_int];
figure (2); notBoxPlot(Grp_integ,Grp,0.5,'patch',ones(length(Grp_integ),1));
%% Ploting wholebrain Segrigation Whol Brain
CNT_wb_seg = mean(CNT_seg,2);ET_Pre_wb_seg = mean(ET_Pre_seg,2);ET_Post_wb_seg = mean(ET_Post_seg,2);
Grp_seg = [CNT_wb_seg; ET_Pre_wb_seg; ET_Post_wb_seg];
figure (3); notBoxPlot(Grp_seg,Grp,0.5,'patch',ones(length(Grp_seg),1));
%% %% Plot the mean connectivity
for i=1:size(CNT_dFC,2)
%CNT_dFC_all(i,:,:) = mean(abs(CNT_dFC {1,i}),3);
CNT_dFC_all(i,:,:) = mean((CNT_dFC_cos {1,i}),3);
end
for i=1:size(ET_Pre_dFC,2)
%ET_Pre_dFC_all(i,:,:) = mean(abs(ET_Pre_dFC {1,i}),3);
ET_Pre_dFC_all(i,:,:) = mean((ET_Pre_dFC_cos {1,i}),3);
end
for i=1:size(ET_Post_dFC,2)
%ET_Post_dFC_all(i,:,:) = mean(abs(ET_Post_dFC {1,i}),3);
ET_Post_dFC_all(i,:,:) = mean((ET_Post_dFC_cos {1,i}),3);
end
%%
close all
figure;
subplot(1,3,1); imagesc(squeeze(mean(CNT_dFC_all,1)));colorbar;caxis([-0.2 0.5])
title ('Mean DFC of Healthy Control')
xlabel('No of ROIs', 'FontName','Calibri', 'FontSize', 14);
ylabel('No of ROIs', 'FontName','Calibri', 'FontSize', 14);
set(gca,'FontSize',14, 'FontName','Calibri', 'FontWeight', 'bold')
%figure; 
subplot(1,3,2);imagesc(squeeze(mean(ET_Pre_dFC_all,1)));colorbar;caxis([-0.2 0.5])
title ('Mean DFC of ET patient before rTMS')
xlabel('No of ROIs', 'FontName','Calibri', 'FontSize', 14);
ylabel('No of ROIs', 'FontName','Calibri', 'FontSize', 14);
set(gca,'FontSize',14, 'FontName','Calibri', 'FontWeight', 'bold')
%figure; 
subplot(1,3,3); imagesc(squeeze(mean(ET_Post_dFC_all,1)));colorbar;caxis([-0.2 0.5])
title ('Mean DFC of ET patient after rTMS')
xlabel('No of ROIs', 'FontName','Calibri', 'FontSize', 11);
ylabel('No of ROIs', 'FontName','Calibri', 'FontSize', 11);
set(gca,'FontSize',14, 'FontName','Calibri', 'FontWeight', 'bold')
CNT_dFC_all_wb=mean(mean(CNT_dFC_all,3),2);
ET_Pre_dFC_all_wb=mean(mean(ET_Pre_dFC_all,3),2);
ET_Post_dFC_all_wb=mean(mean(ET_Post_dFC_all,3),2);
Grp_dFC_wb = [CNT_dFC_all_wb; ET_Pre_dFC_all_wb; ET_Post_dFC_all_wb];
Grp_con = ones(length(CNT_dFC_all_wb),1); Grp_ET_Pre = 2*(ones(length(ET_Pre_dFC_all_wb),1)); Grp_ET_Post = 3*(ones(length(ET_Post_dFC_all_wb),1));
Grp = [Grp_con; Grp_ET_Pre; Grp_ET_Post];
figure; notBoxPlot(Grp_dFC_wb,Grp,0.5,'patch',ones(length(Grp_dFC_wb),1)); title ('Whole brain mean dFC')
xticklabels({'HC','PD NC','PD SH'}); set(gca,'FontSize',14, 'FontName','Times New Roman', 'FontWeight', 'bold')
%%
%Matrix to plot in brain netviewer
CNT_dFC_gm=squeeze(mean(CNT_dFC_all,1));
ET_Pre_dFC_gm=squeeze(mean(ET_Pre_dFC_all,1));
ET_Post_dFC_gm=squeeze(mean(ET_Post_dFC_all,1));
%%save in at text matrix
dlmwrite('HC_dFC_gm.txt',CNT_dFC_gm)

%% whole brain mean and standard deviation of DFC
CNT_dFC_all_m=mean(mean(squeeze(mean(CNT_dFC_all,3)),2))
CNT_dFC_all_s=std(mean(squeeze(mean(CNT_dFC_all,3)),2))
%%%
ET_Pre_dFC_all_m=mean(mean(squeeze(mean(ET_Pre_dFC_all,3)),2))
ET_Pre_dFC_all_s=std(mean(squeeze(mean(ET_Pre_dFC_all,3)),2))
%%%
ET_Post_dFC_all_m=mean(mean(squeeze(mean(ET_Post_dFC_all,3)),2))
ET_Post_dFC_all_s=std(mean(squeeze(mean(ET_Post_dFC_all,3)),2))

%% whole brain mean FCD

CNT_wb_CDC = mean(CNT_CDC,2); ET_Pre_wb_CDC = mean(ET_Pre_CDC,2);ET_Post_wb_CDC = mean(ET_Post_CDC,2);
Grp_CDC = [CNT_wb_CDC; ET_Pre_wb_CDC; ET_Post_wb_CDC];
figure (2); notBoxPlot(Grp_CDC,Grp,0.5,'patch',ones(length(Grp_CDC),1));
%%
CNT_wb_CDC_m = mean(mean(CNT_CDC,2)) 
ET_Pre_wb_CDC_m = mean(mean(ET_Pre_CDC,2))
ET_Post_wb_CDC_m = mean(mean(ET_Post_CDC,2))
[h,p] = ttest2(mean(CNT_CDC,2),mean(ET_Pre_CDC,2))
[h,p] = ttest2(mean(ET_Pre_CDC,2),mean(ET_Post_CDC,2))
