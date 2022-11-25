function [ probingMatrix ] = createProbingMatrix( positions, calcChannels )
% This function creates the "probing matrix" that translates voxel positions 
% into basis functions, i.e. spherical harmonics, for the field expansion.
% Arguments:
%   positions:    voxel positions in mm
%   calcChannels: number of spherical harmonics to evaluate
%
% Returns:
%   probingMatrix: matrix that transforms between voxel positions and coefficients of sperical harmonics

%%
numVoxels = size(positions,1); % positions has size [numVoxels, 3]

% The number of voxels determines up to which order we can expand the measured field
probingMatrix = zeros(numVoxels, calcChannels);
% if numVoxels > 15: expansion up to third order possible -> calcChannels=16
% elseif numVoxels > 8: expansion up to second order possible -> calcChannels=9
% elseif numVoxels > 3: expansion only to first order possible -> calcChannels = 4

if numVoxels > 3
    for p=1:1:numVoxels
        x = positions(p,1);
        y = positions(p,2);
        z = positions(p,3);
        
        % The following formulas for the real-valued sperical harmonics were taken
        % from https://en.wikipedia.org/wiki/Table_of_spherical_harmonics#Real_spherical_harmonics
        probingMatrix(p,1) = 1;
        
        probingMatrix(p,2) = x;
        probingMatrix(p,3) = y;
        probingMatrix(p,4) = z;
        
        if ((numVoxels > 8) && (calcChannels > 4))
            probingMatrix(p,5) = x*y;
            probingMatrix(p,6) = y*z;
            probingMatrix(p,7) = (2*z*z-x*x-y*y);
            probingMatrix(p,8) = z*x;
            probingMatrix(p,9) = (x*x-y*y);
            
            if numVoxels > 15
                probingMatrix(p,10) = (3*x*x-y*y)*y;
                probingMatrix(p,11) = x*y*z;
                probingMatrix(p,12) = (4*z*z-x*x-y*y)*y;
                probingMatrix(p,13) = (2*z*z-3*x*x-3*y*y)*z;
                probingMatrix(p,14) = (4*z*z-x*x-y*y)*x;
                probingMatrix(p,15) = (x*x-y*y)*z;
                probingMatrix(p,16) = (x*x-3*y*y)*x;
            end
        end
    end
end

disp(['probingMatrix successfully created. size: ',num2str(size(probingMatrix))])

end










