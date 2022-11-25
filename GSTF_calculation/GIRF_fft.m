classdef GIRF_fft < handle
% This class defines a GIRF-object, in which the GSTF is calculated by division in the frequency domain.
    properties
        gstf
        girf
        t_axis
        dt
        f_axis
        df
        dwelltime_compensated
        fieldOffsets
    end % properties
    
    methods
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%
        function obj = GIRF_fft(dt)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
            % constructor
            obj.dt = dt;
            obj.dwelltime_compensated = 0;
        end
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function calcH_fft(obj, input, output, skipTriangles)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Calculate the GSTF by dividing the output spectrum by the
            % input spectrum.
            
            % Dimensions:
            % input: [readout-points, triangles]
            % output: [channels, readout-points, triangles]
            % skipTriangles is just a 1D array
            
            obj.gstf = zeros(size(output,2),size(output,1));
            obj.fieldOffsets = zeros(size(output,1),1);
            
            for channel=1:1:size(output,1)
                input_signal = input;
                output_signal = reshape(output(channel,:,:), [size(output,2),size(output,3)]);
                
                % Delete unwanted triangles
                for triang=size(input_signal,2):-1:1
                    if ismember(triang, skipTriangles)
                        input_signal(:,triang) = [];
                        output_signal(:,triang) = [];
                    end
                end % triang
                
                out_fft = fft_1D(output_signal,1);
                in_fft = fft_1D(input_signal,1);
                Norm = sum(abs(in_fft).^2,2);
                % Calculate least-squares-solution from all measured triangles
                obj.gstf(:,channel) = sum((conj(in_fft)).*(out_fft),2)./Norm;
            end % channel
            obj.girf = real(ifft_1D(obj.gstf,1));
            
            % Calculate the time and frequency axes
            F = 1/obj.dt;
            obj.df = F/(size(obj.gstf,1)-1);
            obj.f_axis = (-F/2:obj.df:F/2);
            T = 1/obj.df;
            obj.t_axis = (0:obj.dt:T);
        end % calcH_fft
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function dwelltime_compensation(obj, dwelltime)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Perform the dwell time compensation as described by Stich et
            % al. (DOI: 10.1016/J.MRI.2020.06.005)
            if ~obj.dwelltime_compensated
                obj.gstf = obj.gstf./sinc(dwelltime*obj.f_axis).';
                obj.girf = real(fft_1D(obj.gstf,1));
                obj.dwelltime_compensated = 1;
            end
        end % dwelltime_compensation
             
    end % methods
    
end % classdef







