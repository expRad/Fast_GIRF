clear all;
% Change current directory to that of this .m file
mfile_name          = mfilename('fullpath');
[pathstr,name,ext]  = fileparts(mfile_name);
cd(pathstr);

%% Set file paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax = 'x'; % Choose gradient axis: x, y, or z

if strcmp(ax,'x')
    file2load_meas_noPreemph = '..\trapezoid_measurements\measurement_trap_noPreemphasis_x.mat';
    file2load_meas_preemphHfast = '..\trapezoid_measurements\measurement_trap_withPreemphasis_Hfast_x.mat';
    file2load_meas_preemphHref = '..\trapezoid_measurements\measurement_trap_withPreemphasis_Href_x.mat';
    file2load_meas_preemphHrefHF = '..\trapezoid_measurements\measurement_trap_withPreemphasis_HrefHF_x.mat';
    file2load_Href = '..\GSTF_calculation\results\H_ref_x.mat';
elseif strcmp(ax,'y')
    file2load_meas_noPreemph = '..\trapezoid_measurements\measurement_trap_noPreemphasis_y.mat';
    file2load_meas_preemphHfast = '..\trapezoid_measurements\measurement_trap_withPreemphasis_Hfast_y.mat';
    file2load_meas_preemphHref = '..\trapezoid_measurements\measurement_trap_withPreemphasis_Href_y.mat';
    file2load_meas_preemphHrefHF = '..\trapezoid_measurements\measurement_trap_withPreemphasis_HrefHF_y.mat';
    file2load_Href = '..\GSTF_calculation\results\H_ref_y.mat';
elseif strcmp(ax,'z')
    file2load_meas_noPreemph = '..\trapezoid_measurements\measurement_trap_noPreemphasis_z.mat';
    file2load_meas_preemphHfast = '..\trapezoid_measurements\measurement_trap_withPreemphasis_Hfast_z.mat';
    file2load_meas_preemphHref = '..\trapezoid_measurements\measurement_trap_withPreemphasis_Href_z.mat';
    file2load_meas_preemphHrefHF = '..\trapezoid_measurements\measurement_trap_withPreemphasis_HrefHF_z.mat';
    file2load_Href = '..\GSTF_calculation\results\H_ref_z.mat';
end
file2load_input = '..\trapezoid_measurements\input_trap_noPreemphasis.mat';

%% Load GSTF from file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H_ref_data = load(file2load_Href);
H_ref = H_ref_data.H_combined;
% Normalize GSTF to mean value of 21 data points around omega=0
H_ref_gstf = H_ref.gstf(:,2) / mean(H_ref.gstf(ceil(end/2)-10:ceil(end/2)+10,2),1);

%% Define colors for plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
blue = [0 0.4470 0.7410];
orange = [0.8500 0.3250 0.0980];
violet = [0.4940 0.1840 0.5560];
green = [0.4660 0.6740 0.1880];
yellow = [0.9290 0.6940 0.1250];
lblue = [0.3010 0.7450 0.9330];
dred = [0.6350 0.0780 0.1840];

%% Load gradient input waveform from file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input = load(file2load_input);
index_shift = input.shift;                      % start of the readout relative to start of TR
input = input.grad_input(1:10:end,:);           % gradient waveform
t_axis_input = (0:1e-5:(size(input,1)-1)*1e-5); % time axis for input gradient waveform in seconds

%% Calculate the pre-emphasized gradient waveform %%%%%%%%%%%%%%%%%%%%%%%%%
% Define Gaussian target filter
sigma = 14000;
H_tar = exp(-(H_ref.f_axis.*H_ref.f_axis) / (2*sigma^2)).';

% Determine large enough array length to avoid Gibbs ringing
N1 = 4096;
while (N1 < size(H_ref_gstf,1))
    N1 = N1 * 2;
end
while (N1 < (3 * size(input,1)))
    N1 = N1 * 2;
end
% Bring GSTF, target filter, and input waveform to the desired array length
f_axis_N1 = linspace(H_ref.f_axis(1), H_ref.f_axis(end), N1);   % frequency axis
gstf_N1 = interp1(H_ref.f_axis, H_ref_gstf, f_axis_N1).';       % interpolation
H_tar_N1 = exp(-(f_axis_N1.*f_axis_N1) / (2*sigma^2)).';        % interpolation
extraSamples = round((N1-size(input,1))/2);
if mod(size(input,1),2)                                         % zerop-filling
    input_N1 = [zeros(extraSamples,1);input;zeros(extraSamples-1,1)];
else
    input_N1 = [zeros(extraSamples,1);input;zeros(extraSamples,1)];
end
% Calculate the pre-emphasis
input_fft = fft_1D(input_N1,1);                                 % Fourier-transform input waveform
input_fft_preemph = input_fft.*H_tar_N1./gstf_N1;               % cf. Equations (8) and (10) in the article
input_preemph = real(ifft_1D(input_fft_preemph,1));             % inverse Fourier-transform
if mod(size(input,1),2)                                         % remove extra zeros
    input_preemph = input_preemph(extraSamples+1:end-extraSamples+1);
else
    input_preemph = input_preemph(extraSamples+1:end-extraSamples);
end

%% Determine measured gradient time-courses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numRepPerGrad = 1; numGrad = 1; calcChannels = 2; singleCoil = 0; skipCoils = [];
[ out_grad_noPreemph, dwelltime ] = calcMeasOutput_varPre_wComp_3D( file2load_meas_noPreemph, numRepPerGrad, numGrad, calcChannels, singleCoil, skipCoils);
[ out_grad_preemph_Hfast, ~ ] = calcMeasOutput_varPre_wComp_3D( file2load_meas_preemphHfast, numRepPerGrad, numGrad, calcChannels, singleCoil, skipCoils);
[ out_grad_preemph_Href, ~ ] = calcMeasOutput_varPre_wComp_3D( file2load_meas_preemphHref, numRepPerGrad, numGrad, calcChannels, singleCoil, skipCoils);
[ out_grad_preemph_HrefHF, ~ ] = calcMeasOutput_varPre_wComp_3D( file2load_meas_preemphHrefHF, numRepPerGrad, numGrad, calcChannels, singleCoil, skipCoils);
t_axis_out_grad = (0:dwelltime:(size(out_grad_noPreemph,2)-1)*dwelltime) + (index_shift-1)*1e-6;

%% Plot Figure 8 of the article %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig8 = figure('Units','centimeters', 'InnerPosition',[0 0 17.56 10]);
% Define Positions for subplots
xa = 0.05;
xb = 0.55;
dx = 0.43;
ya = 0.55;
yb = 0.09;
dy = 0.36;

ax1 = subplot('Position',[xa ya dx dy]);
plot(H_ref.f_axis/1000, abs(H_ref_gstf), 'DisplayName',['H^r^e^f_',ax],'Color',yellow,'LineWidth',1);
hold on;
plot(H_ref.f_axis/1000, abs(1./H_ref_gstf), 'DisplayName',['(H^r^e^f_',ax,')^-^1'],'Color',green,'LineWidth',1);
plot(H_ref.f_axis/1000, H_tar,'-.', 'DisplayName','H_t_a_r','Color',blue,'LineWidth',1);
plot(H_ref.f_axis/1000, abs(H_tar./H_ref_gstf),'--', 'DisplayName','H_p_r_e','Color',dred,'LineWidth',1);
leg = legend('Position',[xa 0.92 dx 0.07], 'NumColumns',4);
leg.ItemTokenSize = [15,5];
xlabel('Frequency (kHz)', 'FontSize',8, 'FontName','Times');
ylabel('Magnitude', 'FontSize',8, 'FontName','Times');
xlim([-50 50]);
ylim([0 2]);
set(gca,'FontName','Times','Fontsize',8);
text(-49,1.85,'a','FontName','Arial','Fontsize',12,'FontWeight','bold');

ax2 = subplot('Position',[xb ya dx dy]);
plot(H_ref.f_axis/1000, angle(H_ref_gstf), 'DisplayName',['H^r^e^f_',ax],'Color',yellow,'LineWidth',1);
hold on;
plot(H_ref.f_axis/1000, angle(1./H_ref_gstf), 'DisplayName',['(H^r^e^f_',ax,')^-^1'],'Color',green,'LineWidth',1);
plot(H_ref.f_axis/1000, zeros(size(H_tar)),'-.', 'DisplayName','H_t_a_r','Color',blue,'LineWidth',1);%,'Visible','off');
plot(H_ref.f_axis/1000, angle(H_tar./H_ref_gstf),'--', 'DisplayName','H_p_r_e','Color',dred,'LineWidth',1);
xlabel('Frequency (kHz)', 'FontSize',8, 'FontName','Times');
ylabel('Phase (rad)', 'FontSize',8, 'FontName','Times');
xlim([-50 50]);
set(gca,'FontName','Times','Fontsize',8);
text(-49,3.5,'b','FontName','Arial','Fontsize',12,'FontWeight','bold');

ax3 = subplot('Position',[xa yb 0.27 dy]);
plot(t_axis_input*1000, input, ':', 'DisplayName','nominal','LineWidth',0.7);
hold on;
plot(t_axis_out_grad*1000, out_grad_noPreemph(2,:)*1000,'--', 'DisplayName',sprintf('measured output\n(no pre-emphasis)'),'Color',yellow,'LineWidth',1);
plot(t_axis_input*1000, input_preemph, 'DisplayName',sprintf('input with\npre-emphasis'),'Color',dred,'LineWidth',1);
leg = legend('Location','northeast');
leg.ItemTokenSize = [15,5];
xlim([8.8 10.5]);
ylim([-0.5 6]);
xlabel('Time (ms)');
ylabel('Gradient (mT/m)');
set(gca,'FontName','Times','Fontsize',8);
text(8.84,5.6,'c','FontName','Arial','Fontsize',12,'FontWeight','bold');

ax4 = subplot('Position',[0.38 yb 0.27 dy]);
plot(t_axis_input*1000, input, ':', 'DisplayName','nominal','LineWidth',0.7);
hold on;
plot(t_axis_out_grad*1000, out_grad_noPreemph(2,:)*1000,'--', 'DisplayName',sprintf('measured output\n(no pre-emphasis)'),'Color',yellow,'LineWidth',1);
plot(t_axis_input*1000, input_preemph, 'DisplayName',sprintf('input with\npre-emphasis'),'Color',dred,'LineWidth',1);
xlim([9 9.4]);
ylim([5.2 5.8]);
xlabel('Time (ms)');
set(gca,'FontName','Times','Fontsize',8);
text(9.007,5.76,'d','FontName','Arial','Fontsize',12,'FontWeight','bold');

ax5 = subplot('Position',[0.71 yb 0.27 dy]);
plot(t_axis_input*1000, input, ':', 'DisplayName','nominal','LineWidth',0.7);
hold on;
plot(t_axis_out_grad*1000, out_grad_noPreemph(2,:)*1000,'--', 'DisplayName',sprintf('measured output\n(no pre-emphasis)'),'Color',yellow,'LineWidth',1);
plot(t_axis_input*1000, input_preemph, 'DisplayName',sprintf('input with\npre-emphasis'),'Color',dred,'LineWidth',1);
xlim([8.8 10.5]);
ylim([-0.1 0.08]);
xlabel('Time (ms)');
set(gca,'FontName','Times','Fontsize',8);
text(8.84,0.07,'e','FontName','Arial','Fontsize',12,'FontWeight','bold');


%% Plot Figure 9/S5/S6 of the article %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig9 = figure('Units','centimeters', 'InnerPosition',[0 0 17.56 15]);
% Define Positions for subplots
xa = 0.07;
xb = 0.55;
dx = 0.43;
ya = 0.7;
yb = 0.38;
yc = 0.06;
dy = 0.27;
dya = 0.235;

ax1 = subplot('Position',[xa ya 0.27 dya]);
plot(t_axis_input*1000, input, ':', 'DisplayName','nominal','LineWidth',0.7);
hold on;
plot(t_axis_out_grad*1000, out_grad_noPreemph(2,:)*1000, 'DisplayName','no pre-emphasis','Color',yellow,'LineWidth',1.5);
plot(t_axis_out_grad*1000, out_grad_preemph_Href(2,:)*1000, 'DisplayName','pre-emphasis with H^r^e^f','Color',orange,'LineWidth',1.2);
plot(t_axis_out_grad*1000, out_grad_preemph_HrefHF(2,:)*1000,'-.', 'DisplayName','pre-emphasis with H^r^e^f_H_F','Color',violet,'LineWidth',0.8);
plot(t_axis_out_grad*1000, out_grad_preemph_Hfast(2,:)*1000,'--', 'DisplayName','pre-emphasis with H^f^a^s^t','Color',lblue,'LineWidth',1);
leg = legend('NumColumns',5,'Position',[xa 0.95 0.91 0.04]);
leg.ItemTokenSize = [20,5];
xlim([8.8 10.5]);
ylim([-0.5 6]);
ylabel('Gradient (mT/m)');
set(gca,'FontName','Times','Fontsize',8);
text(8.84,5.6,'a','FontName','Arial','Fontsize',12,'FontWeight','bold');

ax2 = subplot('Position',[0.38 ya 0.27 dya]);
plot(t_axis_input*1000, input, ':', 'DisplayName','nominal','LineWidth',0.7);
hold on;
plot(t_axis_out_grad*1000, out_grad_noPreemph(2,:)*1000, 'DisplayName','no pre-emphasis','Color',yellow,'LineWidth',1.5);
plot(t_axis_out_grad*1000, out_grad_preemph_Href(2,:)*1000, 'DisplayName','pre-emphasis with H^r^e^f','Color',orange,'LineWidth',1.2);
plot(t_axis_out_grad*1000, out_grad_preemph_HrefHF(2,:)*1000,'-.', 'DisplayName','pre-emphasis with H^r^e^f_H_F','Color',violet,'LineWidth',0.8);
plot(t_axis_out_grad*1000, out_grad_preemph_Hfast(2,:)*1000,'--', 'DisplayName','pre-emphasis with H^f^a^s^t','Color',lblue,'LineWidth',1);
xlim([9 9.4]);
ylim([5.2 5.6]);
set(gca,'FontName','Times','Fontsize',8);
text(9.01,5.57,'b','FontName','Arial','Fontsize',12,'FontWeight','bold');
rectangle('Position',[9.05 5.38 0.29 0.04],'LineStyle','-.','EdgeColor',[0.5 0.5 0.5]);
line([9.05;9.1],[5.42;5.57],'LineStyle','-.','Color',[0.5 0.5 0.5]);
line([9.34;9.383],[5.38;5.47],'LineStyle','-.','Color',[0.5 0.5 0.5]);
ax2_inset = axes('Position',[0.45 0.86 0.19 0.06]);
plot(t_axis_input*1000, input, ':', 'DisplayName','nominal','LineWidth',0.7);
hold on;
plot(t_axis_out_grad*1000, out_grad_noPreemph(2,:)*1000, 'DisplayName','no pre-emphasis','Color',yellow,'LineWidth',1.5);
plot(t_axis_out_grad*1000, out_grad_preemph_Href(2,:)*1000, 'DisplayName','pre-emphasis with H^r^e^f','Color',orange,'LineWidth',1.2);
plot(t_axis_out_grad*1000, out_grad_preemph_HrefHF(2,:)*1000,'-.', 'DisplayName','pre-emphasis with H^r^e^f_H_F','Color',violet,'LineWidth',0.8);
plot(t_axis_out_grad*1000, out_grad_preemph_Hfast(2,:)*1000,'--', 'DisplayName','pre-emphasis with H^f^a^s^t','Color',lblue,'LineWidth',1);
xlim([9.05 9.34]);
ylim([5.38 5.42]);
xticklabels('');
yticklabels('');

ax3 = subplot('Position',[0.71 ya 0.27 dya]);
plot(t_axis_input*1000, input, ':', 'DisplayName','nominal','LineWidth',0.7);
hold on;
plot(t_axis_out_grad*1000, out_grad_noPreemph(2,:)*1000, 'DisplayName','no pre-emphasis','Color',yellow,'LineWidth',1.5);
plot(t_axis_out_grad*1000, out_grad_preemph_Href(2,:)*1000, 'DisplayName','pre-emphasis with H^r^e^f','Color',orange,'LineWidth',1.2);
plot(t_axis_out_grad*1000, out_grad_preemph_HrefHF(2,:)*1000,'-.', 'DisplayName','pre-emphasis with H^r^e^f_H_F','Color',violet,'LineWidth',0.8);
plot(t_axis_out_grad*1000, out_grad_preemph_Hfast(2,:)*1000,'--', 'DisplayName','pre-emphasis with H^f^a^s^t','Color',lblue,'LineWidth',1);
xlim([8.8 10.5]);
ylim([-0.02 0.08]);
set(gca,'FontName','Times','Fontsize',8);
text(8.84,0.075,'c','FontName','Arial','Fontsize',12,'FontWeight','bold');
rectangle('Position',[9.4 -0.01 0.3 0.03],'LineStyle','-.','EdgeColor',[0.5 0.5 0.5]);
line([9.4;9.75],[0.02;0.076],'LineStyle','-.','Color',[0.5 0.5 0.5]);
line([9.7;10.4],[-0.01;0.026],'LineStyle','-.','Color',[0.5 0.5 0.5]);
ax2_inset = axes('Position',[0.86 0.815 0.11 0.11]);
plot(t_axis_input*1000, input, ':', 'DisplayName','nominal','LineWidth',0.7);
hold on;
plot(t_axis_out_grad*1000, out_grad_noPreemph(2,:)*1000, 'DisplayName','no pre-emphasis','Color',yellow,'LineWidth',1.5);
plot(t_axis_out_grad*1000, out_grad_preemph_Href(2,:)*1000, 'DisplayName','pre-emphasis with H^r^e^f','Color',orange,'LineWidth',1.2);
plot(t_axis_out_grad*1000, out_grad_preemph_HrefHF(2,:)*1000,'-.', 'DisplayName','pre-emphasis with H^r^e^f_H_F','Color',violet,'LineWidth',0.8);
plot(t_axis_out_grad*1000, out_grad_preemph_Hfast(2,:)*1000,'--', 'DisplayName','pre-emphasis with H^f^a^s^t','Color',lblue,'LineWidth',1);
xlim([9.4 9.7]);
ylim([-0.01 0.02]);
xticklabels('');
yticklabels('');

ax4 = subplot('Position',[xa yb 0.91 dy]);
plot(t_axis_input*1000, input, ':', 'DisplayName','nominal','LineWidth',0.7);
hold on;
plot(t_axis_out_grad*1000, out_grad_noPreemph(2,:)*1000, 'DisplayName','no pre-emphasis','Color',yellow,'LineWidth',1.5);
plot(t_axis_out_grad*1000, out_grad_preemph_Href(2,:)*1000, 'DisplayName','pre-emphasis with H^r^e^f','Color',orange,'LineWidth',1.2);
plot(t_axis_out_grad*1000, out_grad_preemph_HrefHF(2,:)*1000,'-.', 'DisplayName','pre-emphasis with H^r^e^f_H_F','Color',violet,'LineWidth',0.8);
plot(t_axis_out_grad*1000, out_grad_preemph_Hfast(2,:)*1000,'--', 'DisplayName','pre-emphasis with H^f^a^s^t','Color',lblue,'LineWidth',1);
xlim([7 28]);
ylim([-0.02 0.03]);
ylabel('Gradient (mT/m)');
set(gca,'FontName','Times','Fontsize',8);
text(7.2,0.0265,'d','FontName','Arial','Fontsize',12,'FontWeight','bold');

ax5 = subplot('Position',[xa yc 0.91 dy]);
plot(t_axis_input*1000, input, ':', 'DisplayName','nominal','LineWidth',0.7);
hold on;
plot(t_axis_out_grad*1000, out_grad_noPreemph(2,:)*1000, 'DisplayName','no pre-emphasis','Color',yellow,'LineWidth',1.5);
plot(t_axis_out_grad*1000, out_grad_preemph_Href(2,:)*1000, 'DisplayName','pre-emphasis with H^r^e^f','Color',orange,'LineWidth',1.2);
plot(t_axis_out_grad*1000, out_grad_preemph_HrefHF(2,:)*1000,'-.', 'DisplayName','pre-emphasis with H^r^e^f_H_F','Color',violet,'LineWidth',0.8);
plot(t_axis_out_grad*1000, out_grad_preemph_Hfast(2,:)*1000,'--', 'DisplayName','pre-emphasis with H^f^a^s^t','Color',lblue,'LineWidth',1);
xlim([27 48]);
ylim([-0.005 0.005]);
xlabel('Time (ms)');
ylabel('Gradient (mT/m)');
set(gca,'FontName','Times','Fontsize',8);
text(27.2,0.0045,'e','FontName','Arial','Fontsize',12,'FontWeight','bold');














