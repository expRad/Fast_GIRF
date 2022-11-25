function [ out_signals_interp ] = interpolate_output_meas(out_signal, t_axis_ADC, t_axis_1us)
% This function interpolates the measured output gradient to a time grid with stepsize 1 us.
% Arguments:
%   out_signal: measured output gradient
%               dimensions: [channels, ADC readouts, read-out-points, test gradients]
%   t_axis_ADC: time grid of the ADC with step size 1 dwell time, corresponds to out_signal
%   t_axis_1us: time grid with step size of 1 microsecond
%
% Returns:
%   out_signals_interp: output gradient interpolated to t_axis_1us

%%
out_signals_interp = zeros(size(out_signal,1),size(out_signal,2),size(t_axis_1us,2),size(out_signal,4));

for channel=1:1:size(out_signal,1)
    for adc=1:1:size(out_signal,2)
        for triang=1:1:size(out_signal,4)
            if sum(squeeze(out_signal(channel,adc,:,triang)))~=0
                out_signal_new = interp1(t_axis_ADC, squeeze(out_signal(channel,adc,:,triang)), t_axis_1us, 'makima',0);
                out_signal_new(isnan(out_signal_new)) = 0;
                out_signals_interp(channel,adc,:,triang) = out_signal_new;
            end
        end
    end
end

end

