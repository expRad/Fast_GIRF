function  [diff_phase, magnitude_combined] = combineCoils(kspace, dwelltime, weigh_equal, singleCoil, numRepPerGIRF)
% This function combines the magnitude and phase data from different coil elements.
% Arguments:
%   kspace:        complex FID data, 
%                  size of kspace: [numROP, coils, numPE, numPE, slices, triangles, 1, 1, measurements(if>1)]
%   dwelltime:     dwell time in seconds
%   weigh_equal:   1 or 0, whether the different coil elements should be combined by simple averaging
%   singleCoil:    number of the coil element to be evaluated if a combination is not wanted
%   numRepPerGIRF: number of measurement repetitions to be averaged
%
% Returns:
%   diff_phase:         combined derivative of the measured phase data of the different coil elements
%   magnitude_combined: combined magnitude data

%% Extract phase and magnitude data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
diff_phase = unwrap(angle(kspace),[],1); % dimensions: [numROP, coils, numPE, numPE, slices, triangles, 1, 1, measurements(if>1)]
diff_phase(2:end,:,:,:,:,:,:,:,:) = diff(diff_phase,1,1)/dwelltime;
% Add a zero as the first phase derivative value, so the array size stays
% the same. Since the first few data points only record noise, this does
% not affect the gradient triangle.
diff_phase(1,:,:,:,:,:,:,:,:) = 0;

nMeas = size(diff_phase,9);

if singleCoil==0
    % If the used coil has multiple receive elements and we want to combine the data...
    if weigh_equal
        % ...combine the different coil elements by simple averaging if weigh_equal==1.
        diff_phase = mean(diff_phase, 2); % [numROP, 1, numPE, numPE, slices, triangles, 1, 1, measurements(if>1)]
    else
        % ...else, use squared magnitude as weights for the phase derivatives from the different coil elements.
        mag_mean = mean(abs(kspace(:,:,:,:,:,:,:,:,1:numRepPerGIRF)),9); % [numROP, coils, numPE, numPE, slices, triangles, 1, 1, 1]
        weights_squared = mag_mean .* mag_mean;
        clearvars mag_mean;
        sum_mag_squared = sum(weights_squared,2);
        for coil=1:1:size(weights_squared,2)
            weights_squared(:,coil,:,:,:,:,:,:,:) = weights_squared(:,coil,:,:,:,:,:,:,:)./sum_mag_squared;
        end
        clearvars sum_mag_squared;

        weights_squared = repmat(weights_squared, [1,1,1,1,1,1,1,1,nMeas]);
        diff_phase = weights_squared.*diff_phase;
        clearvars weights_squared;
        diff_phase = sum(diff_phase,2);
    end
    % Do a sum-of-squares combination for the magnitude data.
    magnitude_combined = sqrt( sum( abs(kspace).*abs(kspace), 2 ) ); % [numROP, 1, numPE, numPE, slices, triangles, 1, 1, measurements(if>1)]
    clearvars kspace;
else
    % If the used coil has only one receive channel or if we only want to
    % evaluate the data from a specific coil element...
    magnitude_combined = abs(kspace(:,singleCoil,:,:,:,:,:,:,:));
    diff_phase = diff_phase(:,singleCoil,:,:,:,:,:,:,:);
    clearvars kspace;
end

% Discard extra dimensions
magnitude_combined = reshape(magnitude_combined, [size(diff_phase,1),size(diff_phase,3),size(diff_phase,4),size(diff_phase,5),size(diff_phase,6),size(diff_phase,8),size(diff_phase,9)]);
diff_phase = reshape(diff_phase, [size(diff_phase,1),size(diff_phase,3),size(diff_phase,4),size(diff_phase,5),size(diff_phase,6),size(diff_phase,8),size(diff_phase,9)]);
% dimensions: [numROP, numPE, numPE, slices, triangles, 1, measurements(if>1)]

end
