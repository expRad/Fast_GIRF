classdef GIRF_combined < handle
% This class defines a GIRF-object, in which the GSTF is the result of combining H_LF and H_HF.    
    properties
        gstf
        girf
        t_axis
        dt
        f_axis
        df
        fermi_weighting
        cutoffFreq
        lengthH
        fieldOffsets
    end % properties
    
    methods
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = GIRF_combined()
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % constructor
        end
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function combineGSTFs_cutoffFreq(obj, gstf_highRes, f_axis_highRes, gstf_lowRes, f_axis_lowRes, cutoff_freq, fieldOffsets, corr_delay)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Combine low-res and high-res GSTF at a cutoff frequency.
            % First, interpolate low-res GSTF (H_HF) to the high frequency resolution.
            gstf_interp = interp1(f_axis_lowRes, gstf_lowRes, f_axis_highRes);
            gstf_interp(isnan(gstf_interp)) = 0;
            % "Glue" both GSTFs together at the specified frequency.
            obj.cutoffFreq = cutoff_freq;
            weighting = zeros(size(f_axis_highRes)).';
            weighting(abs(f_axis_highRes)<cutoff_freq) = 1;
            obj.gstf = weighting.*gstf_highRes + (1-weighting).*gstf_interp;
            
            % Calculate time and frequency axes
            obj.df = f_axis_highRes(2) - f_axis_highRes(1);
            obj.f_axis = f_axis_highRes;
            obj.dt = 1/(f_axis_highRes(end) - f_axis_highRes(1));
            T = 1/obj.df;
            obj.t_axis = (0:obj.dt:T);
            if size(obj.t_axis,2)<size(obj.gstf,1)
                obj.t_axis = (0:obj.dt:T+obj.dt);
            end
            
            % Calculate combined GIRF
            % This needs some phase shifting to get the impulse at the start of the time axis.
            % (c.f. GIRF_matrix.correct_GSTF_phase())
            linear_array = (1:1:size(obj.gstf,1))*pi;
            corrected_phase = angle(obj.gstf.*exp(1i*(-obj.f_axis*2*pi*corr_delay+linear_array)).');
            gstf_to_girf = abs(obj.gstf).*exp(1i*corrected_phase);
            obj.girf = real(ifft_1D(gstf_to_girf,1));
            obj.fieldOffsets = fieldOffsets;
            
            obj.lengthH = size(obj.gstf,1);
        end % combineGSTFs_cutoffFreq
        
    end % methods
    
end % classdef







