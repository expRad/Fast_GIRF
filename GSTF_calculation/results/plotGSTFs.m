clear all;
% Change current directory to that of this .m file
mfile_name          = mfilename('fullpath');
[pathstr,name,ext]  = fileparts(mfile_name);
cd(pathstr);

%% Load GSTFs from files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H_ref_x_data = load('H_ref_x.mat');
H_ref_y_data = load('H_ref_y.mat');
H_ref_z_data = load('H_ref_z.mat');

H_fast_x_data = load('H_fast_x.mat');
H_fast_y_data = load('H_fast_y.mat');
H_fast_z_data = load('H_fast_z.mat');

H_ref_x = H_ref_x_data.H_combined;
H_ref_y = H_ref_y_data.H_combined;
H_ref_z = H_ref_z_data.H_combined;

H_ref_x_HF = H_ref_x_data.H_HF;
H_ref_y_HF = H_ref_y_data.H_HF;
H_ref_z_HF = H_ref_z_data.H_HF;

H_fast_x = H_fast_x_data.H_combined;
H_fast_y = H_fast_y_data.H_combined;
H_fast_z = H_fast_z_data.H_combined;

%% Define colors for plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
blue = [0 0.4470 0.7410];
orange = [0.8500 0.3250 0.0980];
violet = [0.4940 0.1840 0.5560];
green = [0.4660 0.6740 0.1880];
yellow = [0.9290 0.6940 0.1250];
lblue = [0.3010 0.7450 0.9330];
dred = [0.6350 0.0780 0.1840];

%% Plot Figure 4 of the article %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H_ref_x_FFT = H_ref_x_data.H_fft;
H_ref_x_LF = H_ref_x_data.H_LF;

H_fast_x_FFT = H_fast_x_data.H_fft;
H_fast_x_LF = H_fast_x_data.H_LF;
H_fast_x_HF = H_fast_x_data.H_HF;

fig4 = figure('Units','centimeters', 'InnerPosition',[0.5 0.5 17.56 10]);
% Define Positions for subplots
xa = 0.05;
xb = 0.55;
dx = 0.43;
ya = 0.59;
yb = 0.09;
dy = 0.39;

% Subplot (a): Magnitude of H^ref_FFT,x and H^ref_HF,x
ax1 = subplot('Position',[xa ya dx dy]);
plot(H_ref_x_FFT.f_axis/1000, abs(H_ref_x_FFT.gstf(:,2)), 'DisplayName','H^r^e^f_F_F_T_,_x','Color',green,'LineWidth',1);
hold on;
plot(H_ref_x_HF.f_axis/1000, abs(H_ref_x_HF.gstf(:,2)), 'DisplayName','H^r^e^f_H_F_,_x','Color',violet,'LineWidth',1);
leg = legend('Location','north','NumColumns',2);
leg.ItemTokenSize = [15,5];
xlabel('Frequency (kHz)');
ylabel('Magnitude');
xlim([-51 51]);
ylim([0 1.4]);
xticks([-50 -40 -30 -20 -10 0 10 20 30 40 50]);
set(gca,'FontName','Times','Fontsize',8);
text(-49,1.32,'a','FontName','Arial','Fontsize',12,'FontWeight','bold');

% Subplot (b): Phase of H^ref_FFT,x and H^ref_HF,x
ax2 = subplot('Position',[xb ya dx dy]);
plot(H_ref_x_FFT.f_axis/1000, angle(H_ref_x_FFT.gstf(:,2)), 'DisplayName','H^r^e^f_F_F_T_,_x','Color',green,'LineWidth',1);
hold on;
plot(H_ref_x_HF.f_axis/1000, angle(H_ref_x_HF.gstf(:,2)), 'DisplayName','H^r^e^f_H_F_,_x','Color',violet,'LineWidth',1);
leg = legend('Location','north','NumColumns',2);
leg.ItemTokenSize = [15,5];
xlabel('Frequency (kHz)');
ylabel('Phase (rad)');
xlim([-51 51]);
ylim([-4.1 4.5]);
xticks([-50 -40 -30 -20 -10 0 10 20 30 40 50]);
set(gca,'FontName','Times','Fontsize',8);
text(-50,4,'b','FontName','Arial','Fontsize',12,'FontWeight','bold');

% Subplot (c): Magnitude of H^fast_FFT,x and H^fast_HF,x
ax3 = subplot('Position',[xa yb dx dy]);
plot(H_fast_x_FFT.f_axis/1000, abs(H_fast_x_FFT.gstf(:,2)), 'DisplayName','H^f^a^s^t_F_F_T_,_x','Color',blue,'LineWidth',1);
hold on;
plot(H_fast_x_HF.f_axis/1000, abs(H_fast_x_HF.gstf(:,2)), 'DisplayName','H^f^a^s^t_H_F_,_x','Color',dred,'LineWidth',1);
leg = legend('Location','north','NumColumns',2);
leg.ItemTokenSize = [15,5];
xlabel('Frequency (kHz)');
ylabel('Magnitude');
xlim([-51 51]);
ylim([0 1.4]);
xticks([-50 -40 -30 -20 -10 0 10 20 30 40 50]);
set(gca,'FontName','Times','Fontsize',8);
text(-49,1.32,'c','FontName','Arial','Fontsize',12,'FontWeight','bold');

% Subplot (d): Magnitude of H^fast_LF,x and H^fast_HF,x
ax4 = subplot('Position',[xb yb dx dy]);
plot(H_fast_x_LF.f_axis/1000, abs(H_fast_x_LF.gstf(:,2)), 'DisplayName','H^f^a^s^t_L_F_,_x','Color',yellow,'LineWidth',1);
hold on;
plot(H_fast_x_HF.f_axis/1000, abs(H_fast_x_HF.gstf(:,2)), 'DisplayName','H^f^a^s^t_H_F_,_x','Color',dred,'LineWidth',1);
leg = legend('Location','north','NumColumns',2);
leg.ItemTokenSize = [15,5];
xlabel('Frequency (kHz)');
ylabel('Magnitude');
xlim([-51 51]);
ylim([0 1.4]);
xticks([-50 -40 -30 -20 -10 0 10 20 30 40 50]);
set(gca,'FontName','Times','Fontsize',8);
text(-50,1.32,'d','FontName','Arial','Fontsize',12,'FontWeight','bold');

%% Normalize GSTFs to mean value of 21 data points around omega=0 %%%%%%%%%
H_10T_x_gstf = H_ref_x.gstf(:,2) / mean(H_ref_x.gstf(ceil(end/2)-10:ceil(end/2)+10,2),1);
H_10T_y_gstf = H_ref_y.gstf(:,2) / mean(H_ref_y.gstf(ceil(end/2)-10:ceil(end/2)+10,2),1);
H_10T_z_gstf = H_ref_z.gstf(:,2) / mean(H_ref_z.gstf(ceil(end/2)-10:ceil(end/2)+10,2),1);
H_10T_x_HF_gstf = H_ref_x_HF.gstf(:,2) / mean(H_ref_x_HF.gstf(ceil(end/2)-10:ceil(end/2)+10,2),1);
H_10T_y_HF_gstf = H_ref_y_HF.gstf(:,2) / mean(H_ref_y_HF.gstf(ceil(end/2)-10:ceil(end/2)+10,2),1);
H_10T_z_HF_gstf = H_ref_z_HF.gstf(:,2) / mean(H_ref_z_HF.gstf(ceil(end/2)-10:ceil(end/2)+10,2),1);
H_1T_x_gstf = H_fast_x.gstf(:,2) / mean(H_fast_x.gstf(ceil(end/2)-10:ceil(end/2)+10,2),1);
H_1T_y_gstf = H_fast_y.gstf(:,2) / mean(H_fast_y.gstf(ceil(end/2)-10:ceil(end/2)+10,2),1);
H_1T_z_gstf = H_fast_z.gstf(:,2) / mean(H_fast_z.gstf(ceil(end/2)-10:ceil(end/2)+10,2),1);

%% Plot Figure 5 of the article %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig5 = figure('Units','centimeters', 'InnerPosition',[0.5 0.5 17.56 10]);
% Define positions for subplots
xa = 0.05;
dxa = 0.31;
xb = 0.41;
dxb = 0.33;
xc = 0.79;
dxc = 0.2;
y1 = 0.73;
dy1 = 0.26;
y2 = 0.41;
dy2 = 0.26;
y3 = 0.09;
dy3 = 0.26;

% Subplot (a)-(c): Magnitude of H^ref_x, H^ref_HF,x, H^fast_x
% Subplot (a): Full frequency range
ax1a = subplot('Position',[xa y1 dxa dy1]);
plot(H_ref_x.f_axis/1000, abs(H_10T_x_gstf)+0.2, 'DisplayName','H^r^e^f_x', 'LineWidth',1.2, 'Color',orange);
hold on;
plot(H_ref_x_HF.f_axis/1000, abs(H_10T_x_HF_gstf)+0.1,'-.', 'DisplayName','H^r^e^f_H_F_,_x', 'LineWidth',0.9,'Color',violet);
plot(H_fast_x.f_axis/1000, abs(H_1T_x_gstf),'--', 'DisplayName','H^f^a^s^t_x', 'LineWidth',1,'Color',lblue);
leg = legend('Location','north', 'FontSize',8, 'FontName','Times','NumColumns',3);
leg.ItemTokenSize = [15,5];
ylabel('Magnitude', 'FontSize',8, 'FontName','Times');
xlim([-51 51]);
ylim([0 1.85]);
xticks([-50 -40 -30 -20 -10 0 10 20 30 40 50]);
set(gca,'FontName','Times','Fontsize',8);
text(-49,1.7,'a','FontName','Arial','Fontsize',12,'FontWeight','bold');
rectangle('Position',[-2 0.95 4 0.065],'LineStyle','-.','EdgeColor',[0.5 0.5 0.5]);
text(0,0.84,'b','FontName','Arial','Fontsize',10,'Color',[0.5 0.5 0.5],'HorizontalAlignment','center');

% Subplot (b): Zoom into frequency range between -2 and 2 kHz
ax1b = subplot('Position',[xb y1 dxb dy1]);
plot(H_ref_x.f_axis/1000, abs(H_10T_x_gstf), 'DisplayName','H^r^e^f_x', 'LineWidth',1.2, 'Color',orange);
hold on;
plot(H_ref_x_HF.f_axis/1000, abs(H_10T_x_HF_gstf),'-.', 'DisplayName','H^r^e^f_H_F_,_x', 'LineWidth',0.9,'Color',violet);
plot(H_fast_x.f_axis/1000, abs(H_1T_x_gstf),'--', 'DisplayName','H^f^a^s^t_x', 'LineWidth',1,'Color',lblue);
xlim([-2.1 2.1]);
ylim([0.955 1.01]);
set(gca,'FontName','Times','Fontsize',8);
text(-2,1.005,'b','FontName','Arial','Fontsize',12,'FontWeight','bold');
rectangle('Position',[0.8 0.956 0.7 0.053],'LineStyle','-.','EdgeColor',[0.5 0.5 0.5]);
text(1.55,0.96,'c','FontName','Arial','Fontsize',10,'Color',[0.5 0.5 0.5]);

% Subplot (c): Zoom on mechanical resonances
ax1c = subplot('Position',[xc y1 dxc dy1]);
plot(H_ref_x.f_axis/1000, abs(H_10T_x_gstf), 'DisplayName','H^r^e^f_x', 'LineWidth',1.2, 'Color',orange);
hold on;
plot(H_ref_x_HF.f_axis/1000, abs(H_10T_x_HF_gstf),'-.', 'DisplayName','H^r^e^f_H_F_,_x', 'LineWidth',0.9,'Color',violet);
plot(H_fast_x.f_axis/1000, abs(H_1T_x_gstf),'--', 'DisplayName','H^f^a^s^t_x', 'LineWidth',1,'Color',lblue);
xlim([0.8 1.5]);
ylim([0.955 1.01]);
set(gca,'FontName','Times','Fontsize',8);
text(0.82,1.005,'c','FontName','Arial','Fontsize',12,'FontWeight','bold');

% Subplot (d)-(f): Magnitude of H^ref_y, H^ref_HF,y, H^fast_y
% Subplot (d): Full frequency range
ax2a = subplot('Position',[xa y2 dxa dy2]);
plot(H_ref_y.f_axis/1000, abs(H_10T_y_gstf)+0.2, 'DisplayName','H^r^e^f_y', 'LineWidth',1.2, 'Color',orange);
hold on;
plot(H_ref_y_HF.f_axis/1000, abs(H_10T_y_HF_gstf)+0.1,'-.', 'DisplayName','H^r^e^f_H_F_,_y', 'LineWidth',0.9,'Color',violet);
plot(H_fast_y.f_axis/1000, abs(H_1T_y_gstf),'--', 'DisplayName','H^f^a^s^t_y', 'LineWidth',1,'Color',lblue);
leg = legend('Location','north', 'FontSize',8, 'FontName','Times','NumColumns',3);
leg.ItemTokenSize = [15,5];
ylabel('Magnitude', 'FontSize',8, 'FontName','Times');
xlim([-51 51]);
ylim([0 1.85]);
xticks([-50 -40 -30 -20 -10 0 10 20 30 40 50]);
set(gca,'FontName','Times','Fontsize',8);
text(-49,1.7,'d','FontName','Arial','Fontsize',12,'FontWeight','bold');
rectangle('Position',[-2 0.96 4 0.055],'LineStyle','-.','EdgeColor',[0.5 0.5 0.5]);
text(0,0.88,'e','FontName','Arial','Fontsize',10,'Color',[0.5 0.5 0.5],'HorizontalAlignment','center');

% Subplot (e): Zoom into frequency range between -2 and 2 kHz
ax2b = subplot('Position',[xb y2 dxb dy2]);
plot(H_ref_y.f_axis/1000, abs(H_10T_y_gstf), 'DisplayName','H^r^e^f_y', 'LineWidth',1.2, 'Color',orange);
hold on;
plot(H_ref_y_HF.f_axis/1000, abs(H_10T_y_HF_gstf),'-.', 'DisplayName','H^r^e^f_H_F_,_y', 'LineWidth',0.9,'Color',violet);
plot(H_fast_y.f_axis/1000, abs(H_1T_y_gstf),'--', 'DisplayName','H^f^a^s^t_y', 'LineWidth',1,'Color',lblue);
xlim([-2.1 2.1]);
ylim([0.965 1.01]);
set(gca,'FontName','Times','Fontsize',8);
text(-2,1.005,'e','FontName','Arial','Fontsize',12,'FontWeight','bold');
rectangle('Position',[0.8 0.966 0.7 0.043],'LineStyle','-.','EdgeColor',[0.5 0.5 0.5]);
text(1.55,0.969,'f','FontName','Arial','Fontsize',10,'Color',[0.5 0.5 0.5]);

% Subplot (f): Zoom on mechanical resonances
ax2c = subplot('Position',[xc y2 dxc dy2]);
plot(H_ref_y.f_axis/1000, abs(H_10T_y_gstf), 'DisplayName','H^r^e^f_y', 'LineWidth',1.2, 'Color',orange);
hold on;
plot(H_ref_y_HF.f_axis/1000, abs(H_10T_y_HF_gstf),'-.', 'DisplayName','H^r^e^f_H_F_,_y', 'LineWidth',0.9,'Color',violet);
plot(H_fast_y.f_axis/1000, abs(H_1T_y_gstf),'--', 'DisplayName','H^f^a^s^t_y', 'LineWidth',1,'Color',lblue);
xlim([0.8 1.5]);
ylim([0.965 1.01]);
set(gca,'FontName','Times','Fontsize',8);
text(0.82,1.005,'f','FontName','Arial','Fontsize',12,'FontWeight','bold');

% Subplot (g)-(i): Magnitude of H^ref_z, H^ref_HF,z, H^fast_z
% Subplot (g): Full frequency range
ax3a = subplot('Position',[xa y3 dxa dy3]);
plot(H_ref_z.f_axis/1000, abs(H_10T_z_gstf)+0.2, 'DisplayName','H^r^e^f_z', 'LineWidth',1.2, 'Color',orange);
hold on;
plot(H_ref_z_HF.f_axis/1000, abs(H_10T_z_HF_gstf)+0.1,'-.', 'DisplayName','H^r^e^f_H_F_,_z', 'LineWidth',0.9,'Color',violet);
plot(H_fast_z.f_axis/1000, abs(H_1T_z_gstf),'--', 'DisplayName','H^f^a^s^t_z', 'LineWidth',1,'Color',lblue);
leg = legend('Location','north', 'FontSize',8, 'FontName','Times','NumColumns',3);
leg.ItemTokenSize = [15,5];
xlabel('Frequency (kHz)', 'FontSize',8, 'FontName','Times');
ylabel('Magnitude', 'FontSize',8, 'FontName','Times');
xlim([-51 51]);
ylim([0 1.85]);
xticks([-50 -40 -30 -20 -10 0 10 20 30 40 50]);
set(gca,'FontName','Times','Fontsize',8);
text(-49,1.7,'g','FontName','Arial','Fontsize',12,'FontWeight','bold');
rectangle('Position',[-2 0.955 4 0.08],'LineStyle','-.','EdgeColor',[0.5 0.5 0.5]);
text(0,0.84,'h','FontName','Arial','Fontsize',10,'Color',[0.5 0.5 0.5],'HorizontalAlignment','center');

% Subplot (h): Zoom into frequency range between -2 and 2 kHz
ax3b = subplot('Position',[xb y3 dxb dy3]);
plot(H_ref_z.f_axis/1000, abs(H_10T_z_gstf), 'DisplayName','H^r^e^f_z', 'LineWidth',1.2, 'Color',orange);
hold on;
plot(H_ref_z_HF.f_axis/1000, abs(H_10T_z_HF_gstf),'-.', 'DisplayName','H^r^e^f_H_F_,_z', 'LineWidth',0.9,'Color',violet);
plot(H_fast_z.f_axis/1000, abs(H_1T_z_gstf),'--', 'DisplayName','H^f^a^s^t_z', 'LineWidth',1,'Color',lblue);
xlabel('Frequency (kHz)', 'FontSize',8, 'FontName','Times');
xlim([-2.1 2.1]);
ylim([0.96 1.03]);
set(gca,'FontName','Times','Fontsize',8);
text(-2,1.021,'h','FontName','Arial','Fontsize',12,'FontWeight','bold');
rectangle('Position',[1.1 0.961 0.98 0.068],'LineStyle','-.','EdgeColor',[0.5 0.5 0.5]);
text(0.99,0.967,'i','FontName','Arial','Fontsize',10,'Color',[0.5 0.5 0.5]);

% Subplot (i): Zoom on mechanical resonances
ax3c = subplot('Position',[xc y3 dxc dy3]);
plot(H_ref_z.f_axis/1000, abs(H_10T_z_gstf), 'DisplayName','H^r^e^f_z', 'LineWidth',1.2, 'Color',orange);
hold on;
plot(H_ref_z_HF.f_axis/1000, abs(H_10T_z_HF_gstf),'-.', 'DisplayName','H^r^e^f_H_F_,_z', 'LineWidth',0.9,'Color',violet);
plot(H_fast_z.f_axis/1000, abs(H_1T_z_gstf),'--', 'DisplayName','H^f^a^s^t_z', 'LineWidth',1,'Color',lblue);
xlabel('Frequency (kHz)', 'FontSize',8, 'FontName','Times');
xlim([1.1 2.1]);
ylim([0.96 1.03]);
set(gca,'FontName','Times','Fontsize',8);
text(1.12,1.021,'i','FontName','Arial','Fontsize',12,'FontWeight','bold');













