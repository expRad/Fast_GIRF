classdef preparation < handle

    % Copyright (c) 2022 Hannah Scholten
    
    properties
        % properties that are set in the main program
        resamp_factors
        cut_time
        numRepPerGIRF
        numGIRF
        calcChannels
        numADC
        numTriang
        numDelays
        singleCoil
        skipCoils
        model
        skipTriangles
        numTriangs4fft
        VarPreMeas
        
        % properties that will be determined in prepareData() below
        input
        input_1us
        output
        output_1us
        input4fft
        output4fft
        t_shift
        i_shift
        discardData
        dwelltime
        dts_resamp
        FIDs
        lengthADC_us
        filter
        
    end
    
    methods
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = preparation()
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % constructor
        end
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function prepare_Data(obj, path, names, simpath, simfiles)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.input = cell(1,length(names));
            obj.input_1us = cell(1,length(names));
            obj.output = cell(1,length(names));
            obj.output_1us = cell(1,length(names));
            obj.t_shift = cell(1,length(names));
            obj.i_shift = cell(1,length(names));
            obj.discardData = cell(1,length(names));
            obj.dts_resamp = zeros(1,length(names));
            
            % Iterate through the raw data files of the measurements
            for k=1:1:length(names)
                disp('    Current file: ')
                disp(names{k})

                obj.dts_resamp(k) = obj.resamp_factors(k) * 1e-6;

                % Load the input data
                in = load([simpath,simfiles{k}]);
                in_signal_sim = in.grad_input/1000;  % input gradient, time resolution = 1us
                length_ADC = in.lengthADC;           % length of the readout in us
                obj.lengthADC_us = length_ADC;
                index_shift = in.shift;              % time shifts of the readouts with respect to start of TR
                if size(in_signal_sim,2)<obj.numTriang*obj.numDelays
                    % In case of less test pulses in one measurement, append zeros, so the dimensions match with the rest
                    in_signal_sim = cat(2,in_signal_sim,zeros(size(in_signal_sim,1),obj.numTriang*obj.numDelays-size(in_signal_sim,2)));
                    index_shift = cat(2,index_shift,(zeros(size(index_shift,1),obj.numTriang*obj.numDelays-size(index_shift,2))+1));
                end
                obj.input_1us{k} = in_signal_sim;

                in_signal = in_signal_sim(1:obj.resamp_factors(1):end,:); % Take only the values on the Gradient-Raster-Time-Grid (GRT) (e.g. resamp_factor=10)

                % Determine the output gradient from the phase of the measured signal
                if obj.VarPreMeas % gradient measurement with adapted variable-prephasing method
                    [out_signals, obj.dwelltime ] = calcMeasOutput_varPre_wComp_3D([path,names{k}], obj.numRepPerGIRF, obj.numGIRF, obj.calcChannels, obj.singleCoil, obj.skipCoils);
                    out_signals = reshape(out_signals, [size(out_signals,1),1, size(out_signals,2),1]);
                else % gradient measurement with thin-slice method
                    [out_signals, obj.dwelltime] = calculateOutputGradient([path,names{k}], obj.numRepPerGIRF, obj.numGIRF, obj.numTriang, obj.numDelays, obj.calcChannels, obj.singleCoil, obj.skipCoils);
                end
                % dimensions of out_signals: [channels, numADCs, numDataPoints, triangles]
                
                % Determine time and index shifts of the readouts, needed for the input matrix later on.
                obj.t_shift{k} = (index_shift+obj.cut_time)*1e-6; % because first and last cut_points data points will be discarded;
                obj.i_shift{k} = round(obj.t_shift{k}*1e6/obj.resamp_factors(k));
                
                % Interpolate the output data from the dwelltime-grid to 1us-time-grid, like nominal input data
                t_axis_ADC = (0:obj.dwelltime:(size(out_signals,3)-1)*obj.dwelltime);
                t_axis_1us = (0:1e-6:(length_ADC-1)*1e-6);
                out_signals_interp = interpolate_output_meas(out_signals, t_axis_ADC, t_axis_1us); % [channels, numADCs, lengthADC, triangles]
                clearvars out_signals;
                out_signals_interp = out_signals_interp(:,:,obj.cut_time+1:end-obj.cut_time,:); % discard first and last cut_time data points
                obj.output_1us{k} = out_signals_interp; % [channels, numADCs, lengthADC, triangles]
                
                if k==1
                    % Select data for the GSTF calculation with the FFT-method
                    obj.input4fft = zeros(length_ADC,size(in_signal_sim,2)); % for calculation in frequency domain
                    for i=1:1:size(in_signal_sim,2)
                        % extract nominal input during first ADC for FFT calculation
                        obj.input4fft(:,i) = in_signal_sim(index_shift(1,i):index_shift(1,i)+length_ADC-1,i);
                    end
                    obj.input4fft = obj.input4fft(obj.cut_time+1:end-obj.cut_time,1:obj.numTriangs4fft);
                    % From the output, only take the first ADC and the specified number of triangles for the FFT calculation
                    obj.output4fft = squeeze(out_signals_interp(:,1,:,1:obj.numTriangs4fft));
                    
                    if ~mod(size(obj.input4fft,1),2)
                        % make sure the number of time points is odd so the spectrum and the GSTF are symmetric about 0
                        obj.input4fft = obj.input4fft(2:end,:);
                        obj.output4fft = obj.output4fft(:,2:end,:);
                    end
                end

                % Now resample all output data to the GRT-grid
                out_signal = zeros(size(out_signals_interp,1),size(out_signals_interp,2),ceil(size(out_signals_interp,3)/obj.resamp_factors(k)),size(out_signals_interp,4));
                for channel=1:1:size(out_signal,1)
                    for adc=1:1:size(out_signal,2)
                        [out_signal(channel,adc,:,:), obj.filter] = resample(squeeze(out_signals_interp(channel,adc,:,:)),1,obj.resamp_factors(k));
                    end
                end
                clearvars out_signals_interp;

                obj.output{k} = out_signal(:,1:obj.numADC,:,:);         % only take specified first ADCs (only relevant if more than 1 ADC was measured)
                obj.discardData{k} = zeros(size(obj.output{k})) + 1;    % needed for discarding unusable data points in calcInputMatrix() further down
                obj.input{k} = in_signal;
            end
        end % prepare_Data
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [ max_tshift ] = findMax_tshift(obj)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Find the readout with the maximum time shift.
            % Needed for the determination of the GIRF-time-window in the main program.
            maxs_tshift = zeros(1,length(obj.t_shift));
            for meas=1:1:length(obj.t_shift)
                maxs_tshift(meas) = max(obj.t_shift{meas}(obj.numADC,:),[],'all');
            end
            max_tshift = max(maxs_tshift);
        end % findMax_tshift
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [ inputMatrix ] = calcInputMatrix(obj,numInputChannels,lengthH,weightedMeas_s,weightFactors)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Prepare the input matrix M_in, where each row represents a shifted instance of the input waveform
            Nrows = 0;
            % Determine the number of rows, depending on the number of triangles, 
            % the number of readouts, the number of measurements, and the number of data points per readout
            for meas = 1:1:length(obj.input)
                Nrows = Nrows + obj.numTriang*obj.numDelays * size(obj.output{meas},2) * size(obj.output{meas},3);
            end
            inputMatrix = zeros(Nrows, numInputChannels*lengthH+1); % initialize input matrix
            
            row = 1;
            for meas=1:1:length(obj.input) % iterate through measurements
                in_signal = obj.input{meas};
                time_shift = obj.t_shift{meas};
                dt_out = obj.dts_resamp(meas);
                dt_in = obj.dts_resamp(1);
                numROP = size(obj.output{meas},3);
                for triang=1:1:obj.numTriang*obj.numDelays                      % iterate through triangles
                    for adc=1:1:size(obj.output{meas},2)                        % iterate through ADCs
                        if (~ismember(triang,obj.skipTriangles{meas}))          % check if current triangle should be skipped
                            for rop=1:1:numROP                                  % iterate through Readout points
                                t_out = time_shift(adc,triang) + rop*dt_out;    % absolute time of the current output data point
                                idx_in = round((t_out+1e-4)/dt_in);             % finds the index in the input array corresponding to 100 us ahead
                                                                                % (this improves the numerical stability of the matrix inversion)
                                % Determine start and end indices of the current shifted waveform instance 
                                % in the input array and in the matrix
                                input_start = max(0, (idx_in-lengthH)) +1;
                                input_end = min(idx_in, size(in_signal,1));
                                row_start = max(0, (idx_in-size(in_signal,1))) +1;
                                row_end = min(idx_in, lengthH);
                                
                                for chan=1:1:numInputChannels                   % iterate through number of input channels (in our case 1) (cf. https://doi.org/10.1109/TMI.2019.2936107)
                                    chan_idx = (chan-1)*lengthH;
                                    % Write shifted input waveform into the corresponding entries of the matrix
                                    inputMatrix((row-1)+rop, chan_idx+row_start:chan_idx+row_end) = flip(squeeze(in_signal(input_start:input_end, triang)));
                                end
                                inputMatrix((row-1)+rop, numInputChannels*lengthH+1) = 1; % the last column of the input matrix contains only ones to compute the additive field offsets
                                
                                % Set rows that should not be taken into account to zero
                                inputMatrix((row-1)+rop,:) = inputMatrix((row-1)+rop,:).*obj.discardData{meas}(1,adc,rop,triang);
                                
                                for wMeas = 1:size(weightedMeas_s,2)            % if some measurements should be weighted, apply the weighting factor
                                    if meas==weightedMeas_s(wMeas)
                                        inputMatrix((row-1)+rop,:) = inputMatrix((row-1)+rop,:).*weightFactors(wMeas);
                                    end
                                end
                            end % rop
                        end % skipTriangles
                        row = row + numROP;
                    end % adc
                end % triang
            end % meas
        end % calcInputMatrix
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [ outputMatrix ] = calcOutputMatrix(obj,weightedMeas_s,weightFactors)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Assemble the output matrix 
            % (each column is a vector holding the measured samples of the respective output channel)
            Nrows = 0;
            % As for the input matrix, determine the number of rows first
            for meas = 1:1:length(obj.output)
                Nrows = Nrows + obj.numTriang*obj.numDelays * size(obj.output{meas},2) * size(obj.output{meas},3);
            end
            outputMatrix = zeros(Nrows, obj.calcChannels);
            
            for channel=1:1:obj.calcChannels                                % iterate through output channels (e.g. B0-cross-term & self-term)
                row = 1;
                for meas=1:1:length(obj.output)                             % length(output) is the number of measurements
                    out_signal = obj.output{meas};
                    numROP = size(out_signal,3);
                    for triang=1:1:obj.numTriang*obj.numDelays              % iterate through triangles
                        for adc=1:1:size(obj.output{meas},2)                % iterate through ADCs
                            if (~ismember(triang,obj.skipTriangles{meas}))  % check if current triangle should be skipped
                                outputMatrix(row:row+numROP-1,channel) = squeeze(out_signal(channel,adc,:,triang).*obj.discardData{meas}(channel,adc,:,triang));
                                % By multiplication with obj.discardData, rows that should be ignored in the GIRF calculation are set to zero
                            end
                            for wMeas = 1:size(weightedMeas_s,2)            % if some measurements should be weighted, apply the weighting factor
                                if meas==weightedMeas_s(wMeas)
                                    outputMatrix(row:row+numROP-1,channel) = outputMatrix(row:row+numROP-1,channel).*weightFactors(wMeas);
                                end
                            end
                            row = row + numROP;
                        end % adc
                    end % triang
                end % meas
            end % channel
        end % calcOutputMatrix
        
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [ inputMatrix ] = calcInputMatrix_firstMeas(obj,numInputChannels,lengthH,discard_Data,numFirstMeas)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Create the input matrix for the calculation of H^HF.
            % Only use the first ADC of the first measurement for in this case.
            meas = 1;
            adc = 1;
            Nrows = obj.numTriang * numFirstMeas * size(obj.output{meas},3);
            inputMatrix = zeros(Nrows, numInputChannels*lengthH+1);   % initialize input matrix
            
            row = 1;
            in_signal = obj.input{meas};
            time_shift = obj.t_shift{meas};
            dt_out = obj.dts_resamp(meas);
            dt_in = obj.dts_resamp(1);
            numROP = size(obj.output{meas},3);
            for triang=1:1:obj.numTriang*numFirstMeas
                if ~ismember(triang,obj.skipTriangles{meas})
                    for rop=1:1:numROP
                        t_out = time_shift(adc,triang) + rop*dt_out;    % absolute time of the current output data point
                        idx_in = round((t_out+1e-4)/dt_in);             % finds the index in the input array corresponding to 100 mus ahead

                        input_start = max(0, (idx_in-lengthH)) +1;
                        input_end = min(idx_in, size(in_signal,1));
                        row_start = max(0, (idx_in-size(in_signal,1))) +1;
                        row_end = min(idx_in, lengthH);

                        for chan=1:1:numInputChannels
                            chan_idx = (chan-1)*lengthH;
                            % Write shifted input waveform into the corresponding entries of the matrix
                            inputMatrix((row-1)+rop, chan_idx+row_start:chan_idx+row_end) = flip(squeeze(in_signal(input_start:input_end, triang)));
                        end
                        
                        if discard_Data
                            % Set rows that should not be taken into account to zero
                            inputMatrix((row-1)+rop,:) = inputMatrix((row-1)+rop,:).*obj.discardData{meas}(1,adc,rop,triang);
                        end
                    end % rop
                end % skipTriangles
                row = row + numROP;
            end % triang
        end % calcInputMatrix_firstMeas
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [ outputMatrix ] = calcOutputMatrix_firstMeas(obj,discard_Data,numFirstMeas)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Assemble the output matrix for the calculation of H^HF.
            meas = 1;
            adc = 1;
            Nrows = obj.numTriang * numFirstMeas * size(obj.output{meas},3);
            outputMatrix = zeros(Nrows, obj.calcChannels);
            
            for channel=1:1:obj.calcChannels
                row = 1;
                out_signal = obj.output{meas};
                numROP = size(out_signal,3);
                for triang=1:1:obj.numTriang*numFirstMeas
                    if ~ismember(triang,obj.skipTriangles{meas})
                        if discard_Data
                            outputMatrix(row:row+numROP-1,channel) = squeeze(out_signal(channel,adc,:,triang).*obj.discardData{meas}(channel,adc,:,triang));
                        else
                            outputMatrix(row:row+numROP-1,channel) = squeeze(out_signal(channel,adc,:,triang));
                        end
                    end
                    row = row + numROP;
                end % triang
            end % channel
        end % calcOutputMatrix_firstMeas
        
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [output_calc, residuals] = forwardCalculation_fft(obj, GIRF, chanGIRF)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % This function calculates the GIRF-predicted gradient progression
            % and the difference between measurement and prediction.
            % The calculation is performed in the frequency domain.
            % GIRF: GIRF-object (either GIRF_fft or GIRF_matrix or GIRF_combined)
            % chanGIRF: Which channel of the GIRF to use for the self-term calculation
            %           (This is needed in case the measured data only has one or two channels, 
            %           e.g. B0-cross- and self-term, and the GIRF has more channels.)
            
            % Dimensions:
            % input{i}: [timepoints, numTriang] (timepoints from start of TR until end of last ADC of the ith measurement, in us)
            % output{i}: [channels, ADCs, lengthADC, triangles] (lengthADC in us)
            lengthADC = size(obj.output_1us{1},3);
            num_ADC = size(obj.output_1us{1},2);
            numMeas = length(obj.output_1us);
            numChannelsMeas = size(obj.output_1us{1},1);
            numChannelsGIRF = size(GIRF.gstf, 2);
            
            output_calc = cell(1,numMeas);
            residuals = cell(1,numMeas);
            
            for meas=1:1:numMeas
                input_meas = obj.input_1us{meas};
                output_calc_meas = zeros(numChannelsGIRF, size(input_meas,1), size(input_meas,2));
                
                % Append zeros to avoid side effects by the fft calculation
                nExtra = round((1e6-size(input_meas,1))/2);
                input_meas = [zeros(nExtra,size(input_meas,2)); input_meas; zeros(nExtra,size(input_meas,2))];
                input_meas = fft_1D(input_meas,1);
                for channel=1:numChannelsGIRF
                    % Interpolate the GSTF to the frequency grid corresponding to the zero-filled input array.
                    gstf_interp = interp1(GIRF.f_axis, squeeze(GIRF.gstf(:,channel)), linspace(-500000,500000,size(input_meas,1)),'makima',0);
                    calc_meas = input_meas.*gstf_interp.';                              % multiply input spectrum with GSTF
                    calc_meas = real(ifft_1D(calc_meas,1));                             % transform back to time domain
                    output_calc_meas(channel,:,:) = calc_meas(nExtra+1:end-nExtra,:);   % remove extra zeros
                end
                clearvars input_meas calc_meas gstf_interp;
                
                calc_meas = zeros(size(obj.output_1us{meas}));                          % [channels, ADCs, lengthADC, triangles]
                for channel=1:1:numChannelsMeas
                    for triang=1:1:obj.numTriang*obj.numDelays
                        for adc=1:1:num_ADC
                            idx = round(obj.t_shift{meas}(adc,triang)*1e6);
                            if channel==numChannelsMeas % e.g. 2 -> self-term
                                calc_meas(channel,adc,:,triang) = output_calc_meas(chanGIRF,idx:idx+lengthADC-1,triang) + GIRF.fieldOffsets(chanGIRF);
                            else % e.g. 1 -> B0-cross-term
                                calc_meas(channel,adc,:,triang) = output_calc_meas(channel,idx:idx+lengthADC-1,triang) + GIRF.fieldOffsets(channel);
                            end
                        end % adc
                    end % triang
                end % channel
                output_calc{meas} = calc_meas;
                % Calculate residuals between measured and calculated output signal
                residuals{meas} = obj.output_1us{meas} - calc_meas;
            end % meas
            clearvars output_calc_meas calc_meas;
            
        end % forwardCalculation_fft
        
    end % methods
    
end % classdef







