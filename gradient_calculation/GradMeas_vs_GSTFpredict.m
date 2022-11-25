clear all;
% Change current directory to that of this .m file
mfile_name          = mfilename('fullpath');
[pathstr,name,ext]  = fileparts(mfile_name);
cd(pathstr);

%% Set file paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax = 'x'; % Choose gradient axis: x, y, or z

if strcmp(ax,'x')
    file2load_meas = '..\trapezoid_measurements\measurement_trap_noPreemphasis_x.mat';
    file2load_gstf_1 = '..\GSTF_calculation\results\H_ref_x.mat';
    file2load_gstf_2 = '..\GSTF_calculation\results\H_fast_x.mat';
elseif strcmp(ax,'y')
    file2load_meas = '..\trapezoid_measurements\measurement_trap_noPreemphasis_y.mat';
    file2load_gstf_1 = '..\GSTF_calculation\results\H_ref_y.mat';
    file2load_gstf_2 = '..\GSTF_calculation\results\H_fast_y.mat';
elseif strcmp(ax,'z')
    file2load_meas = '..\trapezoid_measurements\measurement_trap_noPreemphasis_z.mat';
    file2load_gstf_1 = '..\GSTF_calculation\results\H_ref_z.mat';
    file2load_gstf_2 = '..\GSTF_calculation\results\H_fast_z.mat';
end
file2load_input = '..\trapezoid_measurements\input_trap_noPreemphasis.mat';

%% Load the GSTFs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GIRF_data = load(file2load_gstf_1);
H_ref = GIRF_data.H_combined;
H_refHF = GIRF_data.H_HF;
GIRF_data = load(file2load_gstf_2);
H_fast = GIRF_data.H_combined;
% Normalize GSTFs to mean value of 21 data points around omega=0
H_ref.gstf(:,2) = H_ref.gstf(:,2)./mean(abs(H_ref.gstf(ceil(end/2)-10:ceil(end/2)+10,2)),1);
H_fast.gstf(:,2) = H_fast.gstf(:,2)./mean(abs(H_fast.gstf(ceil(end/2)-10:ceil(end/2)+10,2)),1);
H_refHF.gstf(:,2) = H_refHF.gstf(:,2)./mean(abs(H_refHF.gstf(ceil(end/2)-10:ceil(end/2)+10,2)),1);

%% Set some parameters of the gradient measurement %%%%%%%%%%%%%%%%%%%%%%%%
% (For the meaning of some of these parameters, check
% ..\GSTF_calculation\main_H_fast.m or ..\GSTF_calculation\main_H_ref.m.
% The preparation-class was originally designed for the GSTF-calculation.)
prep = preparation;
prep.numRepPerGIRF = 1;         % number of acquired repetitions/averages that should be used for the gradient calculation
prep.numGIRF = 1;               % which iteration number (numIter = numRep/numRepPerGIRF) to use for the GIRF calculation (in case you acquired more averages than you want to use)
prep.numADC = 1;                % number of readouts per TR to evaluate (in case multiple readouts were acquired)
prep.numTriang = 1;             % number of test gradients per delay used in the measurement
prep.numDelays = 1;             % number of delays used in the measurement
prep.skipTriangles = cell(1,1); % Don't leave out any triangles for the forward calculation, leave this an empty cell array!
prep.singleCoil = 0;            % Set this to 0 if a weighted coil combination of the data should be used. Otherwise set this to the number of the coil element to be evaluated.
prep.skipCoils = [];            % skip certain coil elements in case of, for example, too low signal
prep.resamp_factors = [10];     % only relevant for GSTF calculation
prep.cut_time = 80;             % time in mus that is cut at the beginning and end of each measurement period (necessary to eliminate randomly scattered data points)
prep.VarPreMeas = 1;            % Was the measurement done with the variable-prephasing method?
prep.calcChannels = 2;          % number of channels that should be calculated (pay attention at appropriate number of slices)
if prep.calcChannels<2
    disp('ERROR: calcChannels must be >= 2!!!')
    return;
end

%% Preparation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare the input and output gradient data according to the parameters set above
prep.prepare_Data('',{file2load_meas},'',{file2load_input});
dts_out = prep.dts_resamp;

%% Do the forward calculation in the frequency domain %%%%%%%%%%%%%%%%%%%%%
[output_calc_1, residuals_1] = prep.forwardCalculation_fft(H_ref, 2);   % [channels, numADCs, lengthADC, triangles] (lengthADC in us)
[output_calc_2, residuals_2] = prep.forwardCalculation_fft(H_fast, 2);  % [channels, numADCs, lengthADC, triangles] (lengthADC in us)
[output_calc_3, residuals_3] = prep.forwardCalculation_fft(H_refHF, 2); % [channels, numADCs, lengthADC, triangles] (lengthADC in us)

%% Calculate root-mean-square errors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = 1e-6;
T_ADC = size(output_calc_1{1},3) * dt; % time range covered by one ADC in the output
t_axis_ADC = (0:dt:T_ADC-dt);
t_axis_output = prep.t_shift{1}(1,1)+t_axis_ADC;

% Calculate RMSE in the plateau phase of the test gradient
% disp(['t_axis_output(730) = ',num2str(t_axis_output(730)*1000)]);
% disp(['t_axis_output(1130) = ',num2str(t_axis_output(1130)*1000)]);
rmse_1_plateau = sqrt(mean((residuals_1{1}(:,:,731:1130,:)).^2, 3)); % H_ref
rmse_2_plateau = sqrt(mean((residuals_2{1}(:,:,731:1130,:)).^2, 3)); % H_fast
rmse_3_plateau = sqrt(mean((residuals_3{1}(:,:,731:1130,:)).^2, 3)); % H_refHF
disp(['RMSE_plateau(H_ref) = ',num2str(rmse_1_plateau(2,1,1))]);
disp(['RMSE_plateau(H_ref_HF) = ',num2str(rmse_3_plateau(2,1,1))]);
disp(['RMSE_plateau(H_fast) = ',num2str(rmse_2_plateau(2,1,1))]);

diff_meas_nom = prep.input_1us{1}(round(prep.t_shift{1}(1,1)*1e6):end-80,1) - squeeze(prep.output_1us{1}(2,1,:,1));
rmse_4_plateau = sqrt(mean((diff_meas_nom(731:1130,:)).^2, 1));
disp(['RMSE_plateau(nom - meas) = ',num2str(rmse_4_plateau)]);

% Calculate RMSE in the 5ms intervals marked in Figure 6a
idx = [220,5220,10220,15220,20220,25220,30220,35220,39980];
% for i=1:size(idx,2)
%     disp(['t_axis_output(idx',num2str(i),') = ',num2str(t_axis_output(idx(i))*1000)]);
% end
rmse_1 = zeros(size(residuals_1{1},1),size(residuals_1{1},2),size(idx,2)-1,size(residuals_1{1},4));
rmse_2 = zeros(size(rmse_1));
rmse_3 = zeros(size(rmse_1));
rmse_4 = zeros(size(rmse_1));
for i = 1:size(idx,2)-1
    rmse_1(:,:,i,:) = sqrt(mean((residuals_1{1}(:,:,idx(i)+1:idx(i+1),:)).^2, 3));
    rmse_2(:,:,i,:) = sqrt(mean((residuals_2{1}(:,:,idx(i)+1:idx(i+1),:)).^2, 3));
    rmse_3(:,:,i,:) = sqrt(mean((residuals_3{1}(:,:,idx(i)+1:idx(i+1),:)).^2, 3));
    rmse_4(:,:,i,:) = sqrt(mean((diff_meas_nom(idx(i)+1:idx(i+1),:)).^2, 1));
end

%% Plot Figure 6/S1/S3 of the article %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
blue = [0 0.4470 0.7410];
orange = [0.8500 0.3250 0.0980];
violet = [0.4940 0.1840 0.5560];
green = [0.4660 0.6740 0.1880];
yellow = [0.9290 0.6940 0.1250];
lblue = [0.3010 0.7450 0.9330];

dt = 1e-6;
T_ADC = size(output_calc_1{1},3) * dt; % time range covered by one ADC in the output
t_axis_ADC = (0:dt:T_ADC-dt);
t_axis_output = prep.t_shift{1}(1,1)+t_axis_ADC;
T_input = size(prep.input_1us{1},1) * dt;
t_axis_input = (0:dt:T_input-dt);

fig6 = figure('Units','centimeters', 'InnerPosition',[0 0 17.56 19]);

% Subplot (a): Whole time window
ax1 = subplot('Position',[0.07 0.88 0.91 0.08]);
plot(t_axis_input*1000, prep.input_1us{1}(:,1)*1000,'LineStyle',':','DisplayName','nominal','LineWidth',0.7);
hold on;
plot(t_axis_output*1000, squeeze(prep.output_1us{1}(2,1,:,1))*1000, 'DisplayName','measured','Color',yellow,'LineWidth',1.5);
plot(t_axis_output*1000, squeeze(output_calc_1{1}(2,1,:,1))*1000, 'DisplayName','predicted with H^r^e^f','Color',orange,'LineWidth',1.2);
plot(t_axis_output*1000, squeeze(output_calc_3{1}(2,1,:,1))*1000,'-.', 'DisplayName','predicted with H^r^e^f_H_F','LineWidth',0.8);
plot(t_axis_output*1000, squeeze(output_calc_2{1}(2,1,:,1))*1000,'--', 'DisplayName','predicted with H^f^a^s^t','Color',lblue,'LineWidth',1);
for i=1:size(idx,2)
    xline(t_axis_output(idx(i))*1000,'--');
end
labels = {'I','II','III','IV','V','VI','VII','VIII'};
for i=1:size(idx,2)-1
    text(t_axis_output(idx(i))*1000+4.5,5,labels{i},'HorizontalAlignment','right','FontName','Times');
end
labels = {'nominal','measured','predicted with H^r^e^f','predicted with H^r^e^f_H_F','predicted with H^f^a^s^t'};
leg = legend(labels,'NumColumns',5,'Position',[0.07,0.965,0.91,0.028]);
leg.ItemTokenSize = [20,5];
ylabel('Gradient (mT/m)');
ylim([-0.5 6]);
set(gca,'FontName','Times','Fontsize',8);
text(0.3,5.3,'a','FontName','Arial','Fontsize',12,'FontWeight','bold');
rectangle('Position',[8.8 -0.4 1.7 6.3],'LineStyle','-.','EdgeColor',[0.5 0.5 0.5]);
text(10.6,0.9,'b','FontName','Arial','Fontsize',10,'Color',[0.5 0.5 0.5]);

% Subplot (b): Zoom on the trapezoid
ax2 = subplot('Position',[0.07 0.725 0.27 0.12]);
plot(t_axis_input*1000, prep.input_1us{1}(:,1)*1000,'LineStyle',':','DisplayName','nominal','LineWidth',0.7);
hold on;
plot(t_axis_output*1000, squeeze(prep.output_1us{1}(2,1,:,1))*1000, 'DisplayName','measured','Color',yellow,'LineWidth',1.5);
plot(t_axis_output*1000, squeeze(output_calc_1{1}(2,1,:,1))*1000, 'DisplayName','predicted with H^r^e^f','Color',orange,'LineWidth',1.2);
plot(t_axis_output*1000, squeeze(output_calc_3{1}(2,1,:,1))*1000,'-.', 'DisplayName','predicted with H^r^e^f_H_F','LineWidth',0.8);
plot(t_axis_output*1000, squeeze(output_calc_2{1}(2,1,:,1))*1000,'--', 'DisplayName','predicted with H^f^a^s^t','Color',lblue,'LineWidth',1);
ylabel('Gradient (mT/m)');
xlim([8.8 10.5]);
ylim([-0.5 6]);
set(gca,'FontName','Times','Fontsize',8);
text(8.83,5.3,'b','FontName','Arial','Fontsize',12,'FontWeight','bold');
rectangle('Position',[9 5.2 0.4 0.4],'LineStyle','-.','EdgeColor',[0.5 0.5 0.5]);
text(9.42,5.5,'c','FontName','Arial','Fontsize',10,'Color',[0.5 0.5 0.5]);
rectangle('Position',[8.81 -0.1 1.68 0.25],'LineStyle','-.','EdgeColor',[0.5 0.5 0.5]);
text(10.4,0.7,'d','FontName','Arial','Fontsize',10,'Color',[0.5 0.5 0.5]);

% Subplot (c): Zoom on the plateau of the trapezoid
ax3 = subplot('Position',[0.38 0.725 0.27 0.12]);
plot(t_axis_input*1000, prep.input_1us{1}(:,1)*1000,'LineStyle',':','DisplayName','nominal','LineWidth',0.7);
hold on;
plot(t_axis_output*1000, squeeze(prep.output_1us{1}(2,1,:,1))*1000, 'DisplayName','measured','Color',yellow,'LineWidth',1.5);
plot(t_axis_output*1000, squeeze(output_calc_1{1}(2,1,:,1))*1000, 'DisplayName','predicted with H^r^e^f','Color',orange,'LineWidth',1.2);
plot(t_axis_output*1000, squeeze(output_calc_3{1}(2,1,:,1))*1000,'-.', 'DisplayName','predicted with H^r^e^f_H_F','LineWidth',0.8);
plot(t_axis_output*1000, squeeze(output_calc_2{1}(2,1,:,1))*1000,'--', 'DisplayName','predicted with H^f^a^s^t','Color',lblue,'LineWidth',1);
xlim([9 9.4]);
ylim([5.2 5.6]);
set(gca,'FontName','Times','Fontsize',8);
text(9.005,5.57,'c','FontName','Arial','Fontsize',12,'FontWeight','bold');

% Subplot (d): Zoom on the leading and trailing edges of the trapezoid
ax4 = subplot('Position',[0.71 0.725 0.27 0.12]);
plot(t_axis_input*1000, prep.input_1us{1}(:,1)*1000,'LineStyle',':','DisplayName','nominal','LineWidth',0.7);
hold on;
plot(t_axis_output*1000, squeeze(prep.output_1us{1}(2,1,:,1))*1000, 'DisplayName','measured','Color',yellow,'LineWidth',1.5);
plot(t_axis_output*1000, squeeze(output_calc_1{1}(2,1,:,1))*1000, 'DisplayName','predicted with H^r^e^f','Color',orange,'LineWidth',1.2);
plot(t_axis_output*1000, squeeze(output_calc_3{1}(2,1,:,1))*1000,'-.', 'DisplayName','predicted with H^r^e^f_H_F','LineWidth',0.8);
plot(t_axis_output*1000, squeeze(output_calc_2{1}(2,1,:,1))*1000,'--', 'DisplayName','predicted with H^f^a^s^t','Color',lblue,'LineWidth',1);
xlim([8.8 10.5]);
ylim([-0.02 0.08]);
set(gca,'FontName','Times','Fontsize',8);
text(8.83,0.07,'d','FontName','Arial','Fontsize',12,'FontWeight','bold');

% Subplot (e): Zoom on lingering field oscillations during first half of the readout
ax5 = subplot('Position',[0.07 0.535 0.91 0.16]);
plot(t_axis_input*1000, prep.input_1us{1}(:,1)*1000,'LineStyle',':','DisplayName','nominal','LineWidth',0.7);
hold on;
plot(t_axis_output*1000, squeeze(prep.output_1us{1}(2,1,:,1))*1000, 'DisplayName','measured','Color',yellow,'LineWidth',1.5);
plot(t_axis_output*1000, squeeze(output_calc_1{1}(2,1,:,1))*1000, 'DisplayName','predicted with H^r^e^f','Color',orange,'LineWidth',1.2);
plot(t_axis_output*1000, squeeze(output_calc_3{1}(2,1,:,1))*1000,'-.', 'DisplayName','predicted with H^r^e^f_H_F','LineWidth',0.8);
plot(t_axis_output*1000, squeeze(output_calc_2{1}(2,1,:,1))*1000,'--', 'DisplayName','predicted with H^f^a^s^t','Color',lblue,'LineWidth',1);
ylabel('Gradient (mT/m)');
xlim([7 28]);
ylim([-0.025 0.03]);
set(gca,'FontName','Times','Fontsize',8);
text(7.1,0.027,'e','FontName','Arial','Fontsize',12,'FontWeight','bold');

% Subplot (f): Zoom on lingering field oscillations during second half of the readout
ax6 = subplot('Position',[0.07 0.34 0.91 0.16]);
plot(t_axis_input*1000, prep.input_1us{1}(:,1)*1000,'LineStyle',':','DisplayName','nominal','LineWidth',0.7);
hold on;
plot(t_axis_output*1000, squeeze(prep.output_1us{1}(2,1,:,1))*1000, 'DisplayName','measured','Color',yellow,'LineWidth',1.5);
plot(t_axis_output*1000, squeeze(output_calc_1{1}(2,1,:,1))*1000, 'DisplayName','predicted with H^r^e^f','Color',orange,'LineWidth',1.2);
plot(t_axis_output*1000, squeeze(output_calc_3{1}(2,1,:,1))*1000,'-.', 'DisplayName','predicted with H^r^e^f_H_F','LineWidth',0.8);
plot(t_axis_output*1000, squeeze(output_calc_2{1}(2,1,:,1))*1000,'--', 'DisplayName','predicted with H^f^a^s^t','Color',lblue,'LineWidth',1);
xlabel('Time (ms)');
ylabel('Gradient (mT/m)');
xlim([27 48]);
ylim([-0.0051 0.006]);
set(gca,'FontName','Times','Fontsize',8);
text(27.1,0.005,'f','FontName','Arial','Fontsize',12,'FontWeight','bold');

% Subplot (g): RMSEs
ax7 = subplot('Position',[0.07 0.045 0.91 0.24]);
semilogy(squeeze(rmse_4(2,1,:,1))*1000,':x','Color',yellow,'DisplayName','nominal - measured','LineWidth',1.5);
hold on;
semilogy(squeeze(rmse_1(2,1,:,1))*1000,'-*','Color',orange,'LineWidth',1.2);
semilogy(squeeze(rmse_3(2,1,:,1))*1000,'-.o','Color',violet,'LineWidth',0.8);
semilogy(squeeze(rmse_2(2,1,:,1))*1000,'--d','Color',lblue,'LineWidth',1);
labels = {'nominal - measured','measured - predicted with H^r^e^f','measured - predicted with H^r^e^f_H_F','measured - predicted with H^f^a^s^t'};
leg = legend(labels);
leg.ItemTokenSize = [20,5];
xlim([0 9]);
ylim([3e-4 0.11]);
labels = {'','I','II','III','IV','V','VI','VII','VIII'};
xticklabels(labels);
xlabel('Time interval as marked in (a)');
ylabel('RMSE (mT/m)');
set(gca,'FontName','Times','Fontsize',8);
text(0.09,0.087,'g','FontName','Arial','Fontsize',12,'FontWeight','bold');

%% Plot Figure 7/S2/S4 of the article %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate k-space coordinates by integration
gamma = 267.513*10^6; % rad/s/T
k_input = gamma*cumsum(prep.input_1us{1}(:,1), 1)*dt;
k_output = gamma*cumsum(squeeze(prep.output_1us{1}(2,1,:,1)), 1)*dt;
k_calc_1 = gamma*cumsum(squeeze(output_calc_1{1}(2,1,:,1)), 1)*dt;
k_calc_2 = gamma*cumsum(squeeze(output_calc_2{1}(2,1,:,1)), 1)*dt;
k_calc_3 = gamma*cumsum(squeeze(output_calc_3{1}(2,1,:,1)), 1)*dt;
% Calculate RMSEs
residuals_k1 = k_output - k_calc_1;
residuals_k2 = k_output - k_calc_2;
residuals_k3 = k_output - k_calc_3;
residuals_k_nom = k_input(round(prep.t_shift{1}(1,1)*1e6):end-80,1) - k_output;
rmse_k1 = zeros(size(idx,2)-1,1);
rmse_k2 = zeros(size(rmse_k1));
rmse_k3 = zeros(size(rmse_k1));
rmse_k4 = zeros(size(rmse_k1));
for i = 1:size(idx,2)-1
    tmp = residuals_k1(idx(i)+1:idx(i+1));
    rmse_k1(i) = sqrt(mean((tmp.*tmp), 1));
    rmse_k2(i) = sqrt(mean((residuals_k2(idx(i)+1:idx(i+1))).^2, 1));
    rmse_k3(i) = sqrt(mean((residuals_k3(idx(i)+1:idx(i+1))).^2, 1));
    rmse_k4(i) = sqrt(mean((residuals_k_nom(idx(i)+1:idx(i+1))).^2, 1));
end

fig7 = figure('Units','centimeters', 'InnerPosition',[0 0 17.56 19]);

% Subplot (a): Whole time window
ax1 = subplot('Position',[0.07 0.88 0.91 0.08]);
plot(t_axis_input*1000, k_input,'LineStyle',':','DisplayName','nominal','LineWidth',0.7);
hold on;
plot(t_axis_output*1000, k_output, 'DisplayName','measured','Color',yellow,'LineWidth',1.5);
plot(t_axis_output*1000, k_calc_1, 'DisplayName','predicted with H^r^e^f','Color',orange,'LineWidth',1.2);
plot(t_axis_output*1000, k_calc_3,'-.', 'DisplayName','predicted with H^r^e^f_H_F','LineWidth',0.8);
plot(t_axis_output*1000, k_calc_2,'--', 'DisplayName','predicted with H^f^a^s^t','Color',lblue,'LineWidth',1);
for i=1:size(idx,2)
    xline(t_axis_output(idx(i))*1000,'--');
end
labels = {'I','II','III','IV','V','VI','VII','VIII'};
for i=1:size(idx,2)-1
    text(t_axis_output(idx(i))*1000+4.5,50,labels{i},'HorizontalAlignment','right','FontName','Times');
end
labels = {'nominal','measured','predicted with H^r^e^f','predicted with H^r^e^f_H_F','predicted with H^f^a^s^t'};
leg = legend(labels,'NumColumns',5,'Position',[0.07,0.965,0.91,0.028]);
leg.ItemTokenSize = [20,5];
ylabel('k (rad/m)');
ylim([-50 600]);
set(gca,'FontName','Times','Fontsize',8);
text(0.3,500,'a','FontName','Arial','Fontsize',12,'FontWeight','bold');
rectangle('Position',[8.8 -40 1.7 630],'LineStyle','-.','EdgeColor',[0.5 0.5 0.5]);
text(10.6,50,'b','FontName','Arial','Fontsize',10,'Color',[0.5 0.5 0.5]);

% Subplot (b): Zoom on linear increase in k during the trapezoid
ax2 = subplot('Position',[0.07 0.725 0.27 0.12]);
plot(t_axis_input*1000, k_input,'LineStyle',':','DisplayName','nominal','LineWidth',0.7);
hold on;
plot(t_axis_output*1000, k_output, 'DisplayName','measured','Color',yellow,'LineWidth',1.5);
plot(t_axis_output*1000, k_calc_1, 'DisplayName','predicted with H^r^e^f','Color',orange,'LineWidth',1.2);
plot(t_axis_output*1000, k_calc_3,'-.', 'DisplayName','predicted with H^r^e^f_H_F','LineWidth',0.8);
plot(t_axis_output*1000, k_calc_2,'--', 'DisplayName','predicted with H^f^a^s^t','Color',lblue,'LineWidth',1);
ylabel('k (rad/m)');
xlim([8.8 10.5]);
ylim([-50 600]);
set(gca,'FontName','Times','Fontsize',8);
rectangle('Position',[8.6 -20 0.45 40],'LineStyle','-.','EdgeColor',[0.5 0.5 0.5]);
text(9.07,10,'c','FontName','Arial','Fontsize',10,'Color',[0.5 0.5 0.5]);
rectangle('Position',[8.81 469 1.68 43],'LineStyle','-.','EdgeColor',[0.5 0.5 0.5]);
text(10.4,400,'d','FontName','Arial','Fontsize',10,'Color',[0.5 0.5 0.5]);
text(8.83,520,'b','FontName','Arial','Fontsize',12,'FontWeight','bold');

% Subplot (c): Zoom on leading edge of linear increase
ax3 = subplot('Position',[0.38 0.725 0.27 0.12]);
plot(t_axis_input*1000, k_input,'LineStyle',':','DisplayName','nominal','LineWidth',0.7);
hold on;
plot(t_axis_output*1000, k_output, 'DisplayName','measured','Color',yellow,'LineWidth',1.5);
plot(t_axis_output*1000, k_calc_1, 'DisplayName','predicted with H^r^e^f','Color',orange,'LineWidth',1.2);
plot(t_axis_output*1000, k_calc_3,'-.', 'DisplayName','predicted with H^r^e^f_H_F','LineWidth',0.8);
plot(t_axis_output*1000, k_calc_2,'--', 'DisplayName','predicted with H^f^a^s^t','Color',lblue,'LineWidth',1);
xlim([8.6 9.05]);
ylim([-0.41 0.41]);
set(gca,'FontName','Times','Fontsize',8);
text(8.607,0.32,'c','FontName','Arial','Fontsize',12,'FontWeight','bold');

% Subplot (d): Zoom on trailing edge of linear increase
ax4 = subplot('Position',[0.71 0.725 0.27 0.12]);
plot(t_axis_input*1000, k_input,'LineStyle',':','DisplayName','nominal','LineWidth',0.7);
hold on;
plot(t_axis_output*1000, k_output, 'DisplayName','measured','Color',yellow,'LineWidth',1.5);
plot(t_axis_output*1000, k_calc_1, 'DisplayName','predicted with H^r^e^f','Color',orange,'LineWidth',1.2);
plot(t_axis_output*1000, k_calc_3,'-.', 'DisplayName','predicted with H^r^e^f_H_F','LineWidth',0.8);
plot(t_axis_output*1000, k_calc_2,'--', 'DisplayName','predicted with H^f^a^s^t','Color',lblue,'LineWidth',1);
xlim([8.8 10.5]);
ylim([489 492]);
set(gca,'FontName','Times','Fontsize',8);
text(8.83,491.6,'d','FontName','Arial','Fontsize',12,'FontWeight','bold');

% Subplot (e): Zoom on lingering oscillations during first half of the readout
ax5 = subplot('Position',[0.07 0.535 0.91 0.16]);
plot(t_axis_input*1000, k_input,'LineStyle',':','DisplayName','nominal','LineWidth',0.7);
hold on;
plot(t_axis_output*1000, k_output, 'DisplayName','measured','Color',yellow,'LineWidth',1.5);
plot(t_axis_output*1000, k_calc_1, 'DisplayName','predicted with H^r^e^f','Color',orange,'LineWidth',1.2);
plot(t_axis_output*1000, k_calc_3,'-.', 'DisplayName','predicted with H^r^e^f_H_F','LineWidth',0.8);
plot(t_axis_output*1000, k_calc_2,'--', 'DisplayName','predicted with H^f^a^s^t','Color',lblue,'LineWidth',1);
ylabel('k (rad/m)');
xlim([7 28]);
ylim([488 493]);
set(gca,'FontName','Times','Fontsize',8);
text(7.1,492.7,'e','FontName','Arial','Fontsize',12,'FontWeight','bold');

% Subplot (f): Zoom on lingering oscillations during second half of the readout
ax6 = subplot('Position',[0.07 0.34 0.91 0.16]);
plot(t_axis_input*1000, k_input,'LineStyle',':','DisplayName','nominal','LineWidth',0.7);
hold on;
plot(t_axis_output*1000, k_output, 'DisplayName','measured','Color',yellow,'LineWidth',1.5);
plot(t_axis_output*1000, k_calc_1, 'DisplayName','predicted with H^r^e^f','Color',orange,'LineWidth',1.2);
plot(t_axis_output*1000, k_calc_3,'-.', 'DisplayName','predicted with H^r^e^f_H_F','LineWidth',0.8);
plot(t_axis_output*1000, k_calc_2,'--', 'DisplayName','predicted with H^f^a^s^t','Color',lblue,'LineWidth',1);
xlabel('Time (ms)');
ylabel('k (rad/m)');
xlim([27 48]);
ylim([488 493]);
set(gca,'FontName','Times','Fontsize',8);
text(27.1,492.5,'f','FontName','Arial','Fontsize',12,'FontWeight','bold');

k_max = max(k_input,[],'all');
% Subplot (g): Normalized RMSEs
ax7 = subplot('Position',[0.07 0.045 0.91 0.24]);
semilogy(rmse_k4/k_max,':x','Color',yellow,'DisplayName','nominal - measured','LineWidth',1.5);
hold on;
semilogy(rmse_k1/k_max,'-*','Color',orange,'DisplayName','predicted with H^r^e^f','LineWidth',1.2);
semilogy(rmse_k3/k_max,'-.o','Color',violet,'DisplayName','predicted with H^r^e^f_H_F','LineWidth',0.8);
semilogy(rmse_k2/k_max,'--d','Color',lblue,'DisplayName','predicted with H^f^a^s^t','LineWidth',1);
labels = {'nominal - measured','measured - predicted with H^r^e^f','measured - predicted with H^r^e^f_H_F','measured - predicted with H^f^a^s^t'};
if strcmp(ax,'z')
    leg = legend(labels,'NumColumns',2,'Location','southeast');
else
    leg = legend(labels,'NumColumns',2,'Location','northeast');
end
leg.ItemTokenSize = [20,5];
xlim([0 9]);
ylim([1e-5 2e-2]);
labels = {'','I','II','III','IV','V','VI','VII','VIII'};
xticklabels(labels);
xlabel('Time interval as marked in (a)');
ylabel('RMSE/max(k_n_o_m_i_n_a_l)');
set(gca,'FontName','Times','Fontsize',8);
text(0.07,0.015,'g','FontName','Arial','Fontsize',12,'FontWeight','bold');





























    
    