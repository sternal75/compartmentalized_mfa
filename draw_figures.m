close all;

% Sensitivity analysis SHMT1/SHMT2
LF_left_right_BD_NEW=load('LF_left_right_BD.mat');
% LF_left_right_BD_NEW=LF_left_right_BD_NEW.LF_left_right_BD.LF_left_right_BD;
LF_left_right_BD_NEW=LF_left_right_BD_NEW.LF_left_right_BD;
HF_left_right_BD_NEW=load('HF_left_right_BD.mat');
HF_left_right_BD_NEW=HF_left_right_BD_NEW.HF_left_right_BD;

figure
subplot(1,2,1);
x=[0.01 0.2 0.4 0.7 1 2 4 7 11 16 22 27 35 40 55 65 75 85 100 120 140 170 200];
plot(x(LF_left_right_BD_NEW.exitflag_array~=-2),LF_left_right_BD_NEW.error_array(LF_left_right_BD_NEW.exitflag_array~=-2));
set(gca, 'Ylim', [0 20]);
xlabel('SHMT1/(SHMT2-Folate secretion)','FontSize',20);
ylabel('Score','FontSize',26);
title('Low Folic Acid','FontSize',28);
set(gca, 'FontSize', 13);
set(gca,'xtick',[0:1:41])
set(gca,'xticklabel',{'0';'';'';'';'4';'';'';'';'8';'';'';'';'12';'';'';'';'16';'';'';'';'20';'';'';'';'24';'';'';'';'28';'';'';'';'32';'';'';'';'36';'';'';'';'40';''});
set(gca, 'Xlim', [0 41]);
hold on;
plot(x,ones(length(x))*(min(LF_left_right_BD_NEW.error_array)+3.84),'r')
hold off;
LF_under_confidence_interval_ratios = x((LF_left_right_BD_NEW.exitflag_array~=-2) & (LF_left_right_BD_NEW.error_array<(min(LF_left_right_BD_NEW.error_array)+3.84)));
LF_best_ratio = x(find(LF_left_right_BD_NEW.error_array==min(LF_left_right_BD_NEW.error_array)));

subplot(1,2,2);
x=[0.01 0.05 0.1 0.15 0.2 0.25 0.3 0.4 0.6 0.8 1 2 4];
plot(x(HF_left_right_BD_NEW.exitflag_array~=-2),HF_left_right_BD_NEW.error_array(HF_left_right_BD_NEW.exitflag_array~=-2));
set(gca, 'Ylim', [0 35]);
xlabel('SHMT1/(SHMT2-Folate secretion)','FontSize',20);
ylabel('Score','FontSize',26);
title('High Folic Acid','FontSize',28);
set(gca, 'FontSize', 13);
set(gca,'xtick',[0:0.1:2])
set(gca,'xticklabel',{'0';'';'0.2';'';'0.4';'';'0.6';'';'0.8';'';'1';'';'1.2';'';'1.4';'';'1.6';'';'1.8';'';'2'});
set(gca, 'Xlim', [0 2]);
hold on;
plot(x,ones(length(x))*(min(HF_left_right_BD_NEW.error_array)+3.84),'r')
hold off;
HF_under_confidence_interval_ratios = x((HF_left_right_BD_NEW.exitflag_array~=-2) & (HF_left_right_BD_NEW.error_array<(min(HF_left_right_BD_NEW.error_array)+3.84)));
HF_best_ratio = x(find(HF_left_right_BD_NEW.error_array==min(HF_left_right_BD_NEW.error_array)));


% MFA results
% Low Folic Acid
load('model.mat','model');
% title('Low Folic Acid','FontSize',28);
best_LF_index=find(LF_left_right_BD_NEW.error_array==min(LF_left_right_BD_NEW.error_array));
% best_predicted_flux_LF = [196.9499;196.3298;16.6610;102.1009;101.9009;95.7743;95.5743;0.1100;0.0100;0.1100;0.0100;7.9700;0.0100;93.4101;93.6101;143.5015;142.8416;1.1400;0;1.9483;5.3583;1.9483;5.3583;1.1400;9.1211;7.5000;0.1800;6.8400;139.2315;139.2115;137.0798;136.8798;85.1538;84.9538];
% best_predicted_cy_ratio_LF = [0.0000;0.9821];
best_predicted_flux_LF = LF_left_right_BD_NEW.predicted_fluxes_array(:,best_LF_index);
best_predicted_cy_ratio_LF = LF_left_right_BD_NEW.predicted_cy_mt_ratio_array(:,best_LF_index);
% Output(best_predicted_flux_LF, best_predicted_cy_ratio_LF, 0, model, 'MFA - Low Folic Acid');
Output(best_predicted_flux_LF, best_predicted_cy_ratio_LF, 0, model, '');% best_predicted_cy_ratio_HF = [0.0000;0.364];


% best_predicted_flux_HF = [199.9754;196.8654;14.9347;80.1582;71.8655;149.0375;140.7448;4.1564;0.0100;4.1564;0.0100;0.8338;0.0100;30.8743;39.1671;34.4479;34.4479;1.1260;0;1.9038;1.7527;1.9038;1.7527;1.1260;8.9282;9.1165;1.2400;6.7505;37.5221;30.4694;126.5136;118.2208;116.4459;108.1532];
best_HF_index=find(HF_left_right_BD_NEW.error_array==min(HF_left_right_BD_NEW.error_array));
best_predicted_flux_HF = HF_left_right_BD_NEW.predicted_fluxes_array(:,best_HF_index);
best_predicted_cy_ratio_HF = HF_left_right_BD_NEW.predicted_cy_mt_ratio_array(:,best_HF_index);

% Output(best_predicted_flux_HF, best_predicted_cy_ratio_HF, 0, model, 'MFA - High Folic Acid');
Output(best_predicted_flux_HF, best_predicted_cy_ratio_HF, 0, model, '');



figure
hold on
% bar(1:2,[LF_best_ratio HF_best_ratio])
set(gca,'xtick',[1 2]);
set(gca, 'Xlim', [0.5 2.5]);
set(gca,'xticklabel',{'Low Folic Acid';'High Folic Acid'});
ylabel('SHMT1/(SHMT2-Folate secretion) (log2)','FontSize',24);
set(gca, 'FontSize', 24);
hold on
% errorbar(1:2,[LF_best_ratio;HF_best_ratio],[LF_best_ratio-LF_under_confidence_interval_ratios(1) HF_best_ratio-HF_under_confidence_interval_ratios(1)],[LF_under_confidence_interval_ratios(end)-LF_best_ratio HF_under_confidence_interval_ratios(end)-HF_best_ratio],'.','LineWidth',5)
errorbar(1:2,[log2(LF_best_ratio);log2(HF_best_ratio)],[log2(LF_best_ratio)-log2(LF_under_confidence_interval_ratios(1)) log2(HF_best_ratio)-log2(HF_under_confidence_interval_ratios(1))],[log2(LF_under_confidence_interval_ratios(end))-log2(LF_best_ratio) log2(HF_under_confidence_interval_ratios(end))-log2(HF_best_ratio)],'.','LineWidth',5)
set(gca,'yticklabel',{'-8';'-6';'-4';'-2';'0';'2';'4';'6';'8'});
get(gca,'yticklabel')
% set(gca,'YScale','log2')
hold on;
plot([0:0.1:3],zeros(1,length([0:0.1:3])),'black')
hold on;
title({'SHMT1/(SHMT2-Folate secretion)';'95% Confidence Intervals'},'FontSize',33);


SHMT1_values_For_heatmap = [0.5:0.5:5];
SHMT2_values_For_heatmap = [0.5:0.5:5];

figure
% hold on
% bar(1:2,[LF_best_ratio HF_best_ratio])
set(gca,'xtick',[1 2]);
set(gca, 'Xlim', [0.5 2.5]);
set(gca,'xticklabel',{'Low Folic Acid';'High Folic Acid'});
ylabel('SHMT1/(SHMT2-Folate secretion) (log2)','FontSize',14);
set(gca, 'FontSize', 20);
hold on
% errorbar(1:2,[LF_best_ratio;HF_best_ratio],[LF_best_ratio-LF_under_confidence_interval_ratios(1) HF_best_ratio-HF_under_confidence_interval_ratios(1)],[LF_under_confidence_interval_ratios(end)-LF_best_ratio HF_under_confidence_interval_ratios(end)-HF_best_ratio],'.','LineWidth',5)
errorbar(1:2,[log2(LF_best_ratio);log2(HF_best_ratio)],[log2(LF_best_ratio)-log2(LF_under_confidence_interval_ratios(1)) log2(HF_best_ratio)-log2(HF_under_confidence_interval_ratios(1))],[log2(LF_under_confidence_interval_ratios(end))-log2(LF_best_ratio) log2(HF_under_confidence_interval_ratios(end))-log2(HF_best_ratio)],'.','LineWidth',5)
set(gca,'yticklabel',{'-8';'-6';'-4';'-2';'0';'2';'4';'6'});
get(gca,'yticklabel')
% set(gca,'YScale','log2')
hold on;
plot([0:0.1:3],zeros(1,length([0:0.1:3])),'black')
hold on;
title('SHMT1/(SHMT2-Folate secretion) 95% Confidence Intervals','FontSize',20);

RESIZE_FACTOR = 10;
HF_error_matrix_within_confidence_intervals(HF_error_matrix_within_confidence_intervals==1)=2;
HF_error_matrix_within_confidence_intervals_resize = imresize(HF_error_matrix_within_confidence_intervals, size(HF_error_matrix_within_confidence_intervals)*RESIZE_FACTOR)
LF_error_matrix_within_confidence_intervals_resize = imresize(LF_error_matrix_within_confidence_intervals, size(LF_error_matrix_within_confidence_intervals)*RESIZE_FACTOR)
HF_error_matrix_within_confidence_intervals_resize_smooth = zeros(size(HF_error_matrix_within_confidence_intervals)*RESIZE_FACTOR);
LF_error_matrix_within_confidence_intervals_resize_smooth = zeros(size(LF_error_matrix_within_confidence_intervals)*RESIZE_FACTOR);
HF_error_matrix_within_confidence_intervals_resize_smooth(HF_error_matrix_within_confidence_intervals_resize_smooth>=1)=2;
LF_error_matrix_within_confidence_intervals_resize_smooth(LF_error_matrix_within_confidence_intervals_resize_smooth>=0.5)=1;
error_matrix_within_confidence_intervals=HF_error_matrix_within_confidence_intervals+LF_error_matrix_within_confidence_intervals;
figure;
% map=[0.7 0.7 0.7;0.15 0.44 0;0 0.44 0.73;0.7 0 0];
map=[0.7 0.7 0.7;0.15 0.44 0;0 0.44 0.73];
c=colormap(map);
h=heatmap(SHMT2_values_For_heatmap,flip(SHMT1_values_For_heatmap),flip(error_matrix_within_confidence_intervals(:,1:length(SHMT2_values_For_heatmap))),'colormap',c);
h.title({'Low & High Folic Acid';'SHMT1/(SHMT2-Folate secretion) 95% confidence intervals'});
h.xlabel('(SHMT2-Folate secretion)');
h.ylabel('SHMT1');
h.CellLabelColor = 'none';
set(gca, 'FontSize', 16);
set(gca, 'FontSize', 16);
colorbar('off')

