classdef GIRF_matrix < handle
% This class defines a GIRF-object, in which the GIRF is calculated by matrix inversion in the time domain.
    properties
        girf
        gstf
        fieldOffsets
        t_axis
        dt
        f_axis
        df
        dwelltime_compensated
        lengthH
    end % properties
    
    methods
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = GIRF_matrix(dt, lengthH)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % constructor
            obj.dt = dt;
            obj.dwelltime_compensated = 0;
            obj.lengthH = lengthH;
            % Calculate the time and frequency axes
            F = 1/obj.dt;
            obj.df = F/(lengthH-1);
            obj.f_axis = (-F/2:obj.df:F/2);
            T = 1/obj.df;
            obj.t_axis = (0:dt:T);
        end
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function calcH_matrix(obj, inputMatrix, outputMatrix)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Calculate the GIRF by simple matrix inversion, without regularization
            % matrix dimensions:
            % inputMatrix: [timepoints, lengthH + 1]
            % output: [timepoints, channels]
            
            % GIRF calculation without regularization
            H = sparse(inputMatrix)\outputMatrix; % [lengthH + 1, channels]
            H(end,:)
            obj.girf = H(1:end-1,:); % GIRF is contained in all but the last row
            obj.fieldOffsets = H(end,:); % last row contains the field offsets
            
            obj.gstf = fft_1D(obj.girf,1); % [lengthH, channels]
        end % calcH_matrix
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function calcH_matrix_Tikhonov(obj, inputMatrix, outputMatrix, lambda, alpha)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Calculate the GIRF by matrix inversion with a generic
            % Tikhonov regularization
            
            % matrix dimensions:
            % inputMatrix: [timepoints, lengthH + 1]
            % output: [timepoints, channels]
            obj.dwelltime_compensated = 0;
            disp('calcH_matrix_Tikhonov...');
            
            % Prepare regularization matrix
            R = eye(size(inputMatrix,2)); % identity matrix
            % Using the identity matrix for the Tikhonov regularization favors the solution with the smalles norm.
            % Optionally, setting alpha != 0 introduces an exponential weighting into the regularization (not frequency-dependent).
            R(1:obj.lengthH,1:obj.lengthH) = exp(alpha*obj.t_axis.') .* R(1:obj.lengthH,1:obj.lengthH);
            b = zeros(size(inputMatrix,2), size(outputMatrix,2));
            H = sparse([inputMatrix; lambda*R])\[outputMatrix; b]; % [lengthH + 1, channels]
            % https://www.math.uni-frankfurt.de/~harrach/talks/2014Bangalore_Lecture2.pdf
            
            obj.girf = H(1:end-1,:); % GIRF is contained in all but the last row
            obj.fieldOffsets = H(end,:); % last row contains the field offsets
            
            obj.gstf = fft_1D(obj.girf,1); % [lengthH, channels]    
        end % calcH_matrix_Tikhonov
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function calcH_matrix_Tikhonov_freqWeight(obj, inputMatrix, outputMatrix, lambda, alpha_array)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Calculate the GIRF by matrix inversion with the physics-informed
            % Tikhonov regularization
            
            % matrix dimensions:
            % inputMatrix: [timepoints, lengthH + 1]
            % output: [timepoints, channels]
            obj.dwelltime_compensated = 0;
            
            % Prepare regularization matrix
            R = eye(size(inputMatrix,2)); % identity matrix, size=lengthH+1
            % First, calculate Fourier-matrix
            FourierMatrix = zeros(obj.lengthH);
            for column=1:1:obj.lengthH
                omega = 2*pi*obj.f_axis(column);
                FourierMatrix(:,column) = squeeze(exp(1i*omega*obj.t_axis));
            end
            disp('        FourierMatrix created.')
            % Next, introduce different exponential decays for different frequency components
            FourierWithWeighting = FourierMatrix;
            for row=1:1:obj.lengthH
                weights = exp(alpha_array*obj.t_axis(row));
                FourierWithWeighting(row,:) = FourierWithWeighting(row,:).*weights;
            end
            R(1:obj.lengthH,1:obj.lengthH) = (FourierWithWeighting/FourierMatrix);
            clearvars FourierMatrix FourierWithWeighting;
            disp('        Regularization matrix with weighted Fourier operator created.')
            b = zeros(size(inputMatrix,2), size(outputMatrix,2));

            H = real(pinv([inputMatrix; lambda*R])*[outputMatrix; b]); % [lengthH + 1, channels]
            % https://www.math.uni-frankfurt.de/~harrach/talks/2014Bangalore_Lecture2.pdf
            
            obj.girf = H(1:end-1,:); % GIRF is contained in all but the last row
            obj.fieldOffsets = H(end,:); % last row contains the field offsets
            
            obj.gstf = fft_1D(obj.girf,1); % [lengthH, channels]
        end % calcH_matrix_Tikhonov_freqWeight
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function dwelltime_compensation(obj, dwelltime)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Perform the dwell time compensation as described by Stich et
            % al. (DOI: 10.1016/J.MRI.2020.06.005)
            if ~obj.dwelltime_compensated
                obj.gstf = obj.gstf./sinc(dwelltime*obj.f_axis).';
                gstf2 = fft_1D(obj.girf,1);
                gstf2 = gstf2./sinc(dwelltime*obj.f_axis).';
                obj.girf = real(ifft_1D(gstf2,1));
                obj.dwelltime_compensated = 1;
            end
        end % dwelltime_compensation
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function correct_GSTF_phase(obj, phaseAtZero, corr_delay)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Perform a phase correction to compensate for "looking ahead"
            % in the input matrix.
            % Also compensate that the impulse in the GIRF is at t=0, which
            % makes the phase of the FFT of the GIRF oscillate.
            linear_array = (1:1:size(obj.gstf,1))*pi;
            corrected_phase = angle(obj.gstf.*exp(1i*(obj.f_axis*2*pi*corr_delay+linear_array)).');
            corrected_phase_atZero = corrected_phase(floor(size(corrected_phase,1)/2)+3,2);
            if abs(corrected_phase_atZero - phaseAtZero) > 1
                linear_array = (0:1:size(obj.gstf,1)-1)*pi;
            end
            obj.gstf = obj.gstf.*exp(1i*(obj.f_axis*2*pi*corr_delay+linear_array)).';
        end % correct_GSTF_phase
        
    end % methods
    
end % classdef








