function [ out_grad, dwelltime ] = calcMeasOutput_varPre_wComp_3D( file, numRepPerGrad, numGrad, calcChannels, singleCoil, skipCoils)
% This function calculates the measured output gradient time courses
% from the phases of the recorded FID signals
% of the measurements with our modified variable-prephasing sequence
% Arguments:
%   file:          path to the .mat file containing the raw data
%   numRepPerGrad: number of measurement repetitions to be averaged
%   numGrad:       if numRepetitions/numRepPerGrad > 1, which averaged dataset should be evaluated
%   calcChannels:  number of basis functions to be used for the field expansion
%   singleCoil:    number of the coil element to be evaluated if a combination is not wanted
%   skipCoils:     which coil elements to leave out of the calculation (if they have too low signal, for example)
%
% Returns:
%   out_grad:  measured gradient waveform
%   dwelltime: dwell time in seconds

% Copyright (c) 2022 Hannah Scholten

%% Read raw data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rawdata = load(file);

dwelltime = rawdata.dwelltime;      % dwelltime in seconds
FOV = rawdata.FOV;                  % FOV width and height in mm
orientation = rawdata.orientation;  % slice orientation
numSlices = rawdata.numSlices;      % number of slices
PosSlices = rawdata.PosSlices;      % slice positions (order should match the raw data)

kspace = rawdata.kspace;            % FID raw data
% dimensions: [read-out-points, coils, PE-steps(PE-dir), PE-steps(RO-dir), slices, measurement-steps, Cardiac-Phases(1), ADC readouts, repetitions(if>1)]
% If phase encoding was applied in the two directions perpendicular to the
% slice selection direction to enable the determination of higher order
% terms, the respective data should be stored along the dimensions
% "PE-steps(PE-dir)" and "PE-steps(RO-dir)" in the raw data array.
% The following code assumes that the number of PE-steps in both directions
% is equal.
% Since every variable-prephasing (VP) step is played out four times (twice 
% with and twice without the test gradient, and every second time with  
% inverted sign), size(kspace,6) = 4*number of VP steps.

weigh_equal = 0;                    % parameter for the coil weighting

%% Prepare raw data array %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (length(size(kspace))==6)
    % In case of just 1 ADC and 1 measurement: add extra dimensions, so the dimensions in the following are still correct
    kspace=repmat(kspace,1,1,1,1,1,1,1,1,numRepPerGrad);
end
clearvars rawdata;

%% Determine some parameters for easier array reshaping later on %%%%%%%%%%
% Determine number of variable-prephasing steps
% For every VP step, we do 4 measurements:
%   1) with prephasing & test gradient
%   2) like 1) but with inverted signs
%   3) with prephasing gradient and shifted slice selection
%   4) like 3) but with the sign of the prephaser inverted
% (cf Figure 3 in the article)
numVP = size(kspace,6)/4;

numRep = size(kspace,9);                % number of repetitions
numIter = floor(numRep/numRepPerGrad);  % number of loop-counts
numROP = size(kspace,1);                % number of Read Out Points
numPE = size(kspace,3);                 % number of phase encode steps
numADC = size(kspace,8);                % number of ADCs per TR
numMeas = size(kspace,6);               % number of measurements = numVP*4

%% Fourier-transform k-space data along the two phase-encoding directions %
kspace = fft_1D(kspace, 3);
kspace = fft_1D(kspace, 4);

%% Sort out undesired coil elements %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for coil=size(kspace,2):-1:1
    if ismember(coil, skipCoils)
        kspace(:,coil,:,:,:,:,:,:,:) = [];
    end
end

%% Get coil-combined magnitude and phase data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[diff_phase, magnitude] = combineCoils(kspace, dwelltime, weigh_equal, singleCoil, numRepPerGrad);
% dimensions: [numROP, numPE, numPE, numSlices, numMeas, numADC, numRep]
clearvars kspace;

%% Average phase data over numRepPerGrad repetitions %%%%%%%%%%%%%%%%%%%%%%
diff_phase_avg = zeros(numROP, numPE, numPE, numSlices, numMeas, numADC, numIter);
for i=1:1:numIter
    diff_phase_avg(:,:,:,:,:,:,i) = mean(diff_phase(:,:,:,:,:,:,(i-1)*numRepPerGrad+1:i*numRepPerGrad),7);
end
clearvars diff_phase;
% Proceed with the specified iteration
diff_phase_avg = diff_phase_avg(:,:,:,:,:,:,numGrad); % [numROP, numPE, numPE, numSlices, numMeas, numADC, 1]

%% Reorder phase data for the following calculation %%%%%%%%%%%%%%%%%%%%%%%
diff_phase = permute(diff_phase_avg, [2,3,4,5,1,6,7]); % [numPE1, numPE2, numSlices, numMeas, numROP, numADC, 1]
diff_phase = reshape(diff_phase, [numPE, numPE, numSlices, numMeas, numROP*numADC]);
numPositions = numSlices*numPE*numPE;
diff_phase = reshape(diff_phase, [numPositions, numMeas, numROP*numADC]);
% Separate test- and reference measurements
diff_phase_1 = diff_phase(:,1:4:end,:); % [numPositions, numVP, numTimePoints]
diff_phase_2 = diff_phase(:,2:4:end,:); % [numPositions, numVP, numTimePoints]
diff_phase_3 = diff_phase(:,3:4:end,:); % [numPositions, numVP, numTimePoints]
diff_phase_4 = diff_phase(:,4:4:end,:); % [numPositions, numVP, numTimePoints]
clearvars diff_phase

%% Average magnitude data over numRepPerGrad repetitions %%%%%%%%%%%%%%%%%%
% The magnitude is needed to sort out unusable voxels further down
magnitude_avg = zeros(numROP, numPE, numPE, numSlices, numMeas, numADC, numIter);
for i=1:1:numIter
    magnitude_avg(:,:,:,:,:,:,i) = mean(magnitude(:,:,:,:,:,:,(i-1)*numRepPerGrad+1:i*numRepPerGrad),7);
end
magnitude_avg = magnitude_avg(:,:,:,:,:,:,numGrad);
clearvars magnitude;
% Reorder
magnitude = permute(magnitude_avg, [2,3,4,5,1,6,7]); % [numPE, numPE, numSlices, numMeas, numROP, numADC, 1]
magnitude = reshape(magnitude, [numPE, numPE, numSlices, numMeas, numROP*numADC]);
magnitude = reshape(magnitude, [numPositions, numMeas, numROP*numADC]);
magnitude_thres = squeeze(mean(magnitude(:,2,10:50),3));
% Separate test- and reference measurements
magnitude_1 = magnitude(:,1:4:end,:); % [numPositions, numVP, numTimePoints]
magnitude_2 = magnitude(:,2:4:end,:); % [numPositions, numVP, numTimePoints]
magnitude_3 = magnitude(:,3:4:end,:); % [numPositions, numVP, numTimePoints]
magnitude_4 = magnitude(:,4:4:end,:); % [numPositions, numVP, numTimePoints]

%% Get the positions of the measured voxels %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
positions = createPositionArray(orientation, numPE, numSlices, FOV, PosSlices); % [numPE, numPE, numSlices, 3]
positions = reshape(positions, [numPE*numPE*numSlices, 3]); % [numVoxels, 3]

%% Sort out unusable voxels %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by thresholding the signal magnitude
max_mag = max(magnitude_thres);

validVoxels = zeros(size(positions,1),1) + 1;

if numPE>1
    for voxel=numPE*numPE*numSlices:-1:1
        if magnitude_thres(voxel) < max_mag*0.6
            positions(voxel,:) = [];
            diff_phase_1(voxel,:,:) = [];
            diff_phase_2(voxel,:,:) = [];
            diff_phase_3(voxel,:,:) = [];
            diff_phase_4(voxel,:,:) = [];
            validVoxels(voxel) = 0;
        end
    end
end
clearvars magnitude;

%% Build the matrices of the linear system to solve for the field terms %%%
% Reshape arrays
diff_phase_1 = reshape(diff_phase_1, [],numROP*numADC); % [numPositions*numVP, numTimePoints]
diff_phase_2 = reshape(diff_phase_2, [],numROP*numADC); % [numPositions*numVP, numTimePoints]
diff_phase_3 = reshape(diff_phase_3, [],numROP*numADC); % [numPositions*numVP, numTimePoints]
diff_phase_4 = reshape(diff_phase_4, [],numROP*numADC); % [numPositions*numVP, numTimePoints]

magnitude_1 = reshape(magnitude_1, [],numROP*numADC); % [numPositions*numVP, numTimePoints]
magnitude_2 = reshape(magnitude_2, [],numROP*numADC); % [numPositions*numVP, numTimePoints]
magnitude_3 = reshape(magnitude_3, [],numROP*numADC); % [numPositions*numVP, numTimePoints]
magnitude_4 = reshape(magnitude_4, [],numROP*numADC); % [numPositions*numVP, numTimePoints]

% Stack arrays to get the left handside of the matrix equation
diff_phase = cat(1, diff_phase_1, diff_phase_2, diff_phase_3, diff_phase_4); % [4*numPositions*numVP, numTimePoints]
magnitude = cat(1, magnitude_1, magnitude_2, magnitude_3, magnitude_4); % [4*numPositions*numVP, numTimePoints]
clearvars diff_phase_1 diff_phase_2 diff_phase_3 diff_phase_4 magnitude_1 magnitude_2 magnitude_3 magnitude_4

%% Get the probing matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% to translate voxel positions into sperical harmonics
if numPE==1
    probingMatrix = zeros(size(positions,1), calcChannels);
    for slice=1:1:size(positions,1)
        probingMatrix(slice,1) = 1;
        if calcChannels>1
            if strcmp(orientation,'dTra')
                slicePosition = positions(slice,3);
            elseif strcmp(orientation,'dSag')
                slicePosition = positions(slice,1);
            elseif strcmp(orientation,'dCor')
                slicePosition = positions(slice,2);
            end
            probingMatrix(slice,2) = slicePosition;
            if calcChannels>2
                probingMatrix(slice,3) = 2*slicePosition*slicePosition;
            end
        end
    end
else
    probingMatrix = createProbingMatrix(positions,calcChannels); % [numValidVoxels, 4/9/16], depending on the maximum expansion order
    numPositions = size(probingMatrix,1);
end

%% Create the coefficient matrix for the right handside %%%%%%%%%%%%%%%%%%%
%  of the matrix equation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_p = size(probingMatrix,1); % number of positions
N_L = size(probingMatrix,2); % number of basis functions / channels
coeff_mat = zeros(size(diff_phase,1), (1+numVP)*N_L+2*numVP*N_p); % [4*numVP*numPositions, (1+numVP)*numChannels+2*numVP*numPositions]

tmp = repmat(probingMatrix, [numVP,1]);
coeff_mat(1:numVP*N_p,1:N_L) = tmp;
coeff_mat(numVP*N_p+1:2*numVP*N_p,1:N_L) = -tmp;

tmp = zeros(numVP*N_p, numVP*N_L);
for i=0:numVP-1
    tmp(i*N_p+1:(i+1)*N_p, i*N_L+1:(i+1)*N_L) = probingMatrix;
end
coeff_mat(:,N_L+1:(1+numVP)*N_L) = cat(1,tmp,-tmp,tmp,-tmp);

coeff_mat(:,(1+numVP)*N_L+1:(1+numVP)*N_L+numVP*N_p) = cat(1,eye(numVP*N_p),eye(numVP*N_p),zeros(numVP*N_p),zeros(numVP*N_p));
coeff_mat(:,(1+numVP)*N_L+numVP*N_p+1:end) = cat(1,zeros(numVP*N_p),zeros(numVP*N_p),eye(numVP*N_p),eye(numVP*N_p));

clearvars tmp;

%% Calculate the solution matrix... %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ... for each time point separately
sol_mat = zeros(size(coeff_mat,2), size(diff_phase,2)); % [(1+numVP)*numChannels+2*numVP*numPositions, numTimePoints]
for t = 1:size(diff_phase,2)
    f = diff_phase(:,t); % [4*numVP*numPositions, 1]
    m = magnitude(:,t); % [4*numVP*numPositions, 1]
    invcov = zeros(size(f,1));
    for i = 1:size(f,1)
        invcov(i,i) = m(i)*m(i);
    end
    x = (transpose(coeff_mat)*invcov*coeff_mat) \ (transpose(coeff_mat)*invcov*f);
    sol_mat(:,t) = x;
end

%% Extract the output signal from the solution matrix %%%%%%%%%%%%%%%%%%%%%
gamma = 267.513*10^6; % rad/s/T
out_grad = 1/gamma * sol_mat(1:N_L,:); % [numChannels, numTimePoints]

end
