
clear all;
% Change current directory to that of this .m file
mfile_name          = mfilename('fullpath');
[pathstr,name,ext]  = fileparts(mfile_name);
cd(pathstr);

t_tic_1 = tic; % For measuring the time

%% Select data files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path = '../triangle_data_demo/';     % Specifiy the path to the folder where the measurement data are stored

meas_name = 'measurement_H_demo.mat';  % file name of the measurement data
files = dir([path,meas_name]); 
measfiles = {files.name};

ax_name = 'x';  % Specify the axis corresponding to the measurement specified above (for the titles of the plots later on)
% ax_name = 'y';
% ax_name = 'z';

input_path = path;
input_name = 'input_H_demo.mat';         % file name for the input data
files = dir([input_path,input_name]); 
inputfiles = {files.name};

% IMPORTANT: Check that the measurement files and the simulation files correspond to each other!!!
for i=1:1:length(measfiles)
    disp(['measurement file #',num2str(i),': ',measfiles{i}])
    disp(['      input file #',num2str(i),': ',inputfiles{i}])
    disp(' ')
end

doSaveGIRFs = 0;                % Set to 1 if you want to save the result of the calculation, otherwise set to 0
path2save = './results/';       % Specify the folder in which to save the results
name2save = meas_name(13:end);  % Specify the file name

term2plot = 2;                  % Which GSTF term to plot: 1=B0-cross-term, 2=self-term (for higher-order calculation: 2=X, 3=Y, 4=Z, ...)

%% Set parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prep = preparation;
% The preparation object will hold the input and measured output data and
% associated parameters necessary for preparing the data for the GIRF calculation.
prep.numRepPerGIRF = 1;         % number of acquired repetitions/averages that should be used for the GIRF calculation
prep.numGIRF = 1;               % which iteration number (numIter = numRep/numRepPerGIRF) to use for the GIRF calculation (in case you acquired more averages than you want to use)
prep.numADC = 1;                % number of readouts per TR to evaluate (in case multiple readouts were acquired)
prep.numTriang = 10;            % number of triangles per delay used in the measurement
prep.numTriangs4fft = 10;       % number of triangles to use for the calculation of H_FFT
prep.numDelays = 1;             % number of delays used in the measurement
prep.skipTriangles = cell(1,length(measfiles)); % prepare cells to specify which triangles to leave out of the calculation (if needed)
prep.singleCoil = 0;            % Set this to 0 if a weighted coil combination of the data should be used. Otherwise set this to the number of the coil element to be evaluated.
prep.skipCoils = [];            % skip certain coil elements in case of, for example, too low signal
prep.resamp_factors = [10];     % time resolution of the data to be used for the matrix calculation in us (10 us = GRT) (array needs as many entries as measurement files)
prep.cut_time = 80;             % time in mus that is cut at the beginning and end of each measurement period (necessary to eliminate randomly scattered data points)
prep.VarPreMeas = 0;            % Was the measurement done with the variable-prephasing method?
prep.calcChannels = 2;          % number of channels that should be calculated, e.g. 2 = B0-cross- & self-term (pay attention at appropriate number of slices)
if prep.calcChannels<2
    disp('ERROR: calcChannels must be >= 2!!!')
    return;
end

% Some more parameters needed for the GIRF-calculation
corr_delay = 0.95e-4;       % delay to compensate for the phase error 
                            % made by looking ahead in the input matrix (c.f. preparation.calcInputMatrix()) 
                            % and resampling to the GRT time grid (c.f. preparation.prepareData())
lambda = sqrt(1e-7);        % weighting factor for Tikhonov regularization
% Parameters to determine the frequency-dependent growth rate alpha in the
% weighted Fourier operator in the physics-informed regularization
A = 600;
omega_0 = 6000;

%% Preparation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Prepare input and output data...')
% Prepare the input and output gradient data according to the parameters set above
prep.prepare_Data(path,measfiles,input_path,inputfiles);
dts_out = prep.dts_resamp;
prep.findMax_tshift();
% Dimensions:
% prep.input{k}: [numGRTTimePoints, numTriang]
% prep.output{k}: [calcChannels, numADC, numGRTTimePoints, numTriang]
% prep.input4fft: [numADCTimePoints, numTriang]
% prep.output4fft: [calcChannels, numADCTimePoints, numTriang]

disp('    Data prepared.')

%% Calculate H_FFT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Calculate H_FFT...')
H_fft = GIRF_fft(1e-6);
H_fft.calcH_fft(prep.input4fft, prep.output4fft, prep.skipTriangles{1})
% H_fft.gstf has dimensions [lengthADC, calcChannels]
disp('    H_FFT calculated.')

%% Put input and output data into matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Prepare matrices for time domain calculation...')
t_tic_2 = tic;

numInputChannels = 1;       % number of input channels
lengthADC = size(prep.output{1},3);
dt_in = dts_out(1);
% Determine length of the GIRF based on the measured time windows
lengthH = floor( (prep.findMax_tshift() - prep.t_shift{1}(1,1))/dt_in + lengthADC -990/prep.resamp_factors(1) );
% Make sure lengthH is odd (such that the central point of the GSTF is at f=0)
if ~mod(lengthH,2)
    lengthH = lengthH - 1;
end

% Put the input and output gradient data in matrix form
% Matrices for the calculation of H_LF (contain all data)
inputMatrix_LF = prep.calcInputMatrix(numInputChannels, lengthH,[],[]);
% Dimension inputMatrix: [time points, lengthH+1]
outputMatrix_LF = prep.calcOutputMatrix([],[]);
% Dimension outputMatrix: [time points, channels]

% Matrices for the calculation of H_HF (contain only data from the measurements with delay = 0)
lengthH_HF = floor( lengthADC -990/prep.resamp_factors(1) ) -200;
if ~mod(lengthH_HF,2)
    lengthH_HF = lengthH_HF - 1;
end
inputMatrix_HF = prep.calcInputMatrix_firstMeas(numInputChannels, lengthH_HF, 0, prep.numTriangs4fft/prep.numTriang);
outputMatrix_HF = prep.calcOutputMatrix_firstMeas(0, prep.numTriangs4fft/prep.numTriang);

disp(['    size(inputMatrix_LF) = ',num2str(size(inputMatrix_LF))])
disp(['    size(outputMatrix_LF) = ',num2str(size(outputMatrix_LF))])
disp(['    size(inputMatrix_HF) = ',num2str(size(inputMatrix_HF))])
disp(['    size(outputMatrix_HF) = ',num2str(size(outputMatrix_HF))])
disp('    Matrices prepared.')

%% Calculate GIRFs in the time domain with the matrix inversion method %%%%
disp('Calculate H_LF and H_HF in time domain...')
H_LF = GIRF_matrix(dt_in,lengthH);                                      % GIRF with high frequency resolution, used for low frequencies (LF) < omega_0
H_HF = GIRF_matrix(dt_in,lengthH_HF);                            % GIRF with low frequency resolution, used for high frequencies (HF)

% Calculate H_LF with generic Tikhonov regularization
H_LF.calcH_matrix_Tikhonov(inputMatrix_LF, outputMatrix_LF, lambda, 0);
% Calculate H_HF with physics-informed regularization
alpha_array = alpha_func(H_HF.f_axis, omega_0, H_HF.f_axis(end), A);    % alpha_func is defined at the end of the script
H_HF.calcH_matrix_Tikhonov_freqWeight(inputMatrix_HF, outputMatrix_HF, lambda, alpha_array);

%% Phase correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Correct the phase of the GSTF to account for the timing error introduced
% by "looking ahead" in the input matrix (corr_delay is defined in the preparation part above)
phaseAtZero_fft = angle(H_fft.gstf(floor(size(H_fft.gstf,1)/2)+3,2));
H_LF.correct_GSTF_phase(phaseAtZero_fft, corr_delay);
H_HF.correct_GSTF_phase(phaseAtZero_fft, corr_delay);

disp('    H in time domain calculated.')

%% Dwell time compensation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Perform dwelltime compensation...')

H_fft.dwelltime_compensation(prep.dwelltime);
H_LF.dwelltime_compensation(prep.dwelltime);
H_HF.dwelltime_compensation(prep.dwelltime);

disp('    Dwell time compensation done.');

%% Combine H_LF and H_HF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Combine GSTFs...')
H_combined = GIRF_combined();
H_combined.combineGSTFs_cutoffFreq(H_LF.gstf, H_LF.f_axis, H_HF.gstf, H_HF.f_axis, omega_0, H_LF.fieldOffsets, corr_delay);

disp('    GSTFs combined.')
t_toc_2 = toc(t_tic_2);

%% Save GIRFs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if doSaveGIRFs
    save([path2save,name2save],'H_combined','lambda','alpha_array','H_LF','H_HF','H_fft','dts_out','corr_delay','omega_0');
    disp('Saved results.');
else
    disp('Results were not saved.');
end
t_toc_1 = toc(t_tic_1);

%% Plot and compare results of H_fft and H_matrix %%%%%%%%%%%%%%%%%%%%%%%%%
disp('Plotting...')
plotter = GIRF_plotter();
plotter.plot_GSTFs(H_fft,['H^d^e^m^o_F_F_T_,_',ax_name], H_LF,['H^d^e^m^o_L_F_,_',ax_name], H_HF,['H^d^e^m^o_H_F_,_',ax_name], H_combined,['H^d^e^m^o_',ax_name], 4, term2plot, ax_name);
disp('    main_H_demo.m finished.')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function alpha = alpha_func(omega, omega_0, omega_max, A)
    alpha = zeros(size(omega));
    for i=1:size(omega,2)
        if abs(omega(i))<omega_0
            alpha(i) = 0;
        else
            alpha(i) = A*sqrt((abs(omega(i))-omega_0)/(omega_max-omega_0));
        end
    end
end






