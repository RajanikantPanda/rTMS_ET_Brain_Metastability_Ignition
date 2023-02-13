clc
clear

CNT = load('C:\Research\RAJ_RP_BlueHD\NMS\ET\fMRI\Dyn\Ignition\4TR_SS_Ignition_CASE1_Phases_WB.mat');
ET_Pre = load('C:\Research\RAJ_RP_BlueHD\NMS\ET\fMRI\Dyn\Ignition\4TR_SS_Ignition_CASE2_Phases_WB.mat');
ET_Post = load('C:\Research\RAJ_RP_BlueHD\NMS\ET\fMRI\Dyn\Ignition\4TR_SS_Ignition_CASE3_Phases_WB.mat');
% Bar plot of mean ignition AN
Grp_con = ones(length(CNT.mignitionAN),1); Grp_ET_Pre = 2*ones(length(ET_Pre.mignitionAN),1); Grp_ET_Post = 3*ones(length(ET_Post.mignitionAN),1);
Grp = [Grp_con; Grp_ET_Pre; Grp_ET_Post]; %Grp = Gzarp';
clear Grp_con Grp_ET_Pre Grp_ET_Post
Grp_mignitionAN = [CNT.mignitionAN'; ET_Pre.mignitionAN'; ET_Post.mignitionAN'];
figure(6); notBoxPlot(Grp_mignitionAN,Grp,0.5,'patch',ones(length(Grp_mignitionAN),1));
title ('Meaan event ignition Whole Brain')
%% Mean ignition across the events of all subjects
figure(1)
plot(CNT.mevGL, 'black')
hold on
plot(ET_Pre.mevGL, 'blue')
plot(ET_Post.mevGL, 'red')
title ('Mean ignition across the events of all subjects');
% STD ignition across the events of all subjects
figure(2)
plot(CNT.stdevGL, 'black')
hold on
plot(ET_Pre.stdevGL, 'blue')
plot(ET_Post.stdevGL, 'red')
title ('STD ignition across the events of all subjects');

%%
ci = 2.576; %99% = 2.576 / 98%	= 2.326 / 95%	= 1.96
figure(3)
a1= squeeze(median(CNT.mevokedintegSS,3));
a1s= ci*squeeze(std(CNT.mevokedintegSS,[],3))/sqrt (length(CNT.mevokedintegSS));
a2= squeeze(median(ET_Pre.mevokedintegSS,3));
a2s= ci*squeeze(std(ET_Pre.mevokedintegSS,[],3))/sqrt (length(ET_Pre.mevokedintegSS));
a3= squeeze(median(ET_Post.mevokedintegSS,3));
a3s= ci*squeeze(std(ET_Post.mevokedintegSS,[],3))/sqrt (length(ET_Post.mevokedintegSS));
errorbar (a1,a1s, 'k', 'LineWidth',2); grid on;
hold on
errorbar (a2,a2s, 'b', 'LineWidth',2); grid on;
errorbar (a3,a3s, 'r', 'LineWidth',2); grid on;
%title("Visual-2 mean evoked integration", 'FontSize', 18);
title("Whole Brain Ignition Driven Mean Integration (IDMI)", 'FontSize', 18);
legend('Control', 'ET Pre rTMS', 'ET Post rTMS', 'FontWeight', 'bold')
ylabel('Intrinsic Ignition', 'FontSize', 18);
xlabel('No of ROIs', 'FontSize', 18);
set(gca,'FontSize',18, 'FontName','Times New Roman')

create_shennifti(a1,'CNT_resion_wise_IDMI.nii.gz',0)
create_shennifti(a2,'PrerTMS_resion_wise_IDMI.nii.gz',0)
create_shennifti(a3,'PostrTMS_resion_wise_IDMI.nii.gz',0)
%%
figure(14)
a1= squeeze(mean(CNT.stdevokedintegSS,3));
a1s= ci*squeeze(std(CNT.stdevokedintegSS,[],3))/sqrt (length(CNT.stdevokedintegSS));
a2= squeeze(mean(ET_Pre.stdevokedintegSS,3));
a2s= ci*squeeze(std(ET_Pre.stdevokedintegSS,[],3))/sqrt (length(ET_Pre.stdevokedintegSS));
a3= squeeze(mean(ET_Post.stdevokedintegSS,3));
a3s= ci*squeeze(std(ET_Post.stdevokedintegSS,[],3))/sqrt (length(ET_Post.stdevokedintegSS));
errorbar (a1,a1s, 'k', 'LineWidth',2); grid on;
hold on
errorbar (a2,a2s, 'b', 'LineWidth',2); grid on;
errorbar (a3,a3s, 'r', 'LineWidth',2); grid on;
title("STD ignition across the events of all subjects", 'FontSize', 18);
legend('Control', 'ET pre rTMS', 'ET post rTMS', 'FontWeight', 'bold')
ylabel('Metastability', 'FontSize', 18);
xlabel('No of ROIs', 'FontSize', 18);
set(gca,'FontSize',14, 'FontName','Times New Roman')
%%
figure(14)
a1= squeeze(median(CNT.stdevokedintegSS,3));
a1s= ci*squeeze(std(CNT.stdevokedintegSS,[],3))/sqrt (length(CNT.stdevokedintegSS));
a2= squeeze(median(ET_Pre.stdevokedintegSS,3));
a2s= ci*squeeze(std(ET_Pre.stdevokedintegSS,[],3))/sqrt (length(ET_Pre.stdevokedintegSS));
a3= squeeze(median(ET_Post.stdevokedintegSS,3));
a3s= ci*squeeze(std(ET_Post.stdevokedintegSS,[],3))/sqrt (length(ET_Post.stdevokedintegSS));
errorbar (a1,a1s, 'k', 'LineWidth',2); grid on;
hold on
errorbar (a2,a2s, 'b', 'LineWidth',2); grid on;
errorbar (a3,a3s, 'r', 'LineWidth',2); grid on;
title("Whole Brain Ignition Driven Mean Integration (IDMI)", 'FontSize', 18);
legend('Control', 'ET pre rTMS', 'ET post rTMS', 'FontWeight', 'bold')
ylabel('Intrinsic Ignition', 'FontSize', 18);
xlabel('No of ROIs', 'FontSize', 18);
set(gca,'FontSize',14, 'FontName','Times New Roman')
%%
figure(5)
plot(CNT.mignitionAN, 'black')
hold on
plot(ET_Pre.mignitionAN, 'blue')
plot(ET_Post.mignitionAN, 'red')

%% Bar plot of mean ignition AN
Grp_con = ones(length(CNT.mignitionAN),1); Grp_ET_Pre = 2*ones(length(ET_Pre.mignitionAN),1); Grp_ET_Post = 3*ones(length(ET_Post.mignitionAN),1);
Grp = [Grp_con; Grp_ET_Pre; Grp_ET_Post]; %Grp = Gzarp';
clear Grp_con Grp_ET_Pre Grp_ET_Post
Grp_mignitionAN = [CNT.mignitionAN'; ET_Pre.mignitionAN'; ET_Post.mignitionAN'];
figure(6); notBoxPlot(Grp_mignitionAN,Grp,0.5,'patch',ones(length(Grp_mignitionAN),1));
title ('Meaan event ignition Whole Brain')
%% Bar plot of standard daviation of evoked ignition
b1= squeeze(mean(CNT.stdevokedintegSS,2));
b2= squeeze(mean(ET_Pre.stdevokedintegSS,2));
b3= squeeze(mean(ET_Post.stdevokedintegSS,2));
Grp_mevokedinteg = [b1; b2; b3];
figure(7); notBoxPlot(Grp_mevokedinteg,Grp,0.5,'patch',ones(length(Grp_mevokedinteg),1));
title ('Variance event ignition Whole Brain')
%% Bar plot of mean ignition AN
Grp_con = ones(length(CNT.mignitionAN),1); Grp_ET_Pre = 2*ones(length(ET_Pre.mignitionAN),1); Grp_ET_Post = 3*ones(length(ET_Post.mignitionAN),1);
Grp = [Grp_con; Grp_ET_Pre; Grp_ET_Post]; %Grp = Gzarp';
clear Grp_con Grp_ET_Pre Grp_ET_Post
Grp_mignitionAN = [CNT.mignitionAN'; ET_Pre.mignitionAN'; ET_Post.mignitionAN'];
figure(6); notBoxPlot(Grp_mignitionAN,Grp,0.5,'patch',ones(length(Grp_mignitionAN),1));
title ('Meaan event ignition Whole Brain')
%% Bar plot of standard daviation of evoked ignition
b1= squeeze(mean(CNT.stdevokedintegSS,2));
b2= squeeze(mean(ET_Pre.stdevokedintegSS,2));
b3= squeeze(mean(ET_Post.stdevokedintegSS,2));
Grp_mevokedinteg = [b2; b3];
figure(7); notBoxPlot(Grp_mevokedinteg,Grp,0.5,'patch',ones(length(Grp_mevokedinteg),1));
title ('Variance event ignition Whole Brain')

%% Metastability, intigration and segrigation
Grp_meta = [ET_Pre_meta'; ET_Post_meta'];
figure(8); notBoxPlot(Grp_meta,Grp,0.5,'patch',ones(length(Grp_meta),1));
title ('Whole Brain Metastability')
CNT_wb_integ = mean(CNT_integ,2);ET_Pre_wb_integ = mean(ET_Pre_integ,2);ET_Post_wb_integ = mean(ET_Post_integ,2);
Grp_integ = [ET_Pre_wb_integ; ET_Post_wb_integ];
figure (9); notBoxPlot(Grp_integ,Grp,0.5,'patch',ones(length(Grp_integ),1));
title ('Whole Brain Intigration')
CNT_wb_Q = mean(CNT_Q,2);ET_Pre_wb_Q = mean(ET_Pre_Q,2);ET_Post_wb_Q = mean(ET_Post_Q,2);
Grp_Q = [ET_Pre_wb_Q; ET_Post_wb_Q];
figure (10); notBoxPlot(Grp_Q,Grp,0.5,'patch',ones(length(Grp_Q),1));
title ('Whole Brain Segrigation')