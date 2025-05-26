function [bSample, bRef, bImagingAndRef] = CAIPIRINHAMasks(N1, N2, ...
                                                           Acceleration1, ...
                                                           Acceleration2, ...
                                                           RefLines1, ...
                                                           RefLines2, ...
                                                           CAIPIShift, ...
                                                           dimfast)
% Inputs:
%   N1, N2           - Size of k-space (rows x cols)
%   Acceleration1    - Acceleration factor along PE1 (typically slower dir)
%   Acceleration2    - Acceleration factor along PE2 (typically faster dir)
%   RefLines1,2      - Centered reference lines in each dimension
%   CAIPIShift       - Integer shift for CAIPIRINHA
%   dimfast          - 'PE' or '3D', determines CAIPI dimension ordering
% created by Jinyuan Zhang, 2025.05.25

    if nargin < 8
        dimfast = 'PE'; % default
    end
    
    % === Handle dimension flip if dimfast is '3D'
    bFlip = strcmpi(dimfast, '3D');
    if bFlip
        % Swap all parameters related to dimensions
        [N1           , N2           ] = deal(N2           , N1           );
        [Acceleration1, Acceleration2] = deal(Acceleration2, Acceleration1);
        [RefLines1    , RefLines2    ] = deal(RefLines2    , RefLines1    );
    end
    
    Dir1CenterIdx = Ordering.Utils.FFTCenterIndex(N1);
    Dir2CenterIdx = Ordering.Utils.FFTCenterIndex(N2);
    
    Dir1RefCenterIdx = Ordering.Utils.FFTCenterIndex(RefLines1);
    Dir2RefCenterIdx = Ordering.Utils.FFTCenterIndex(RefLines2);
    
    Dir1RefIdxMin = max(1 , Dir1CenterIdx + (1         - Dir1RefCenterIdx));
    Dir1RefIdxMax = min(N1, Dir1CenterIdx + (RefLines1 - Dir1RefCenterIdx));
    
    Dir2RefIdxMin = max(1 , Dir2CenterIdx + (1         - Dir2RefCenterIdx));
    Dir2RefIdxMax = min(N2, Dir2CenterIdx + (RefLines2 - Dir2RefCenterIdx));
    
    % Initialize masks
    bRef     = false(N1, N2);
    bImaging = false(N1, N2);
    
    % Build reference region mask
    bRef(Dir1RefIdxMin:Dir1RefIdxMax, Dir2RefIdxMin:Dir2RefIdxMax) = true;
    
    % Create CAIPIRINHA sampling mask
    for i = 1:N1
        if mod(i - Dir1CenterIdx, Acceleration1) == 0
            shift = mod((i - Dir1CenterIdx)/Acceleration1 * CAIPIShift, Acceleration2);
            for j = 1:N2
                if mod(j - Dir2CenterIdx - shift, Acceleration2) == 0
                  bImaging(i,j) = true;
                end
            end
        end
    end
    
    % Combine masks
    bSample        = bImaging | bRef;
    bImagingAndRef = bImaging & bRef;
    
    % === Flip back if necessary
    if bFlip
        bSample        = bSample.'       ;
        bRef           = bRef.'          ;
        bImagingAndRef = bImagingAndRef.';
    end
end
