function [seq, Label] = prep_NoiseScan(seq, Actual, PE3D, ADC, Label, sys)
    adc        = ADC.adc;
    %% Define sequence blocks
    % Start with a noise scan - only on the first repetition.
    % The line and partition of the noise scan must be on one of the "imaging"
    % scans. For simplicity we set it as the center of k-space because that
    % is always(?) sampled.
    
    % Round up adc total duration (including delays and dead time) to be on
    % the block duration raster time. Hopefully in the future we can remove the
    % extra dummy delay below of duration NoiseBlockDur. (This delay is in
    % parallel to the rest of the events in the block, not before them.)
    if strcmpi(Actual.NoiseScan, 'on')
        NoiseBlockDur = sys.blockDurationRaster * ceil(mr.calcDuration(adc)/sys.blockDurationRaster - 1e-6) ;
        seq.addBlock(adc, ...
                     mr.makeDelay(NoiseBlockDur), ... % Temp bug fix!!!
                     mr.makeLabel('SET', 'ONCE', 1), ... % only on 1st repetition.
                     mr.makeLabel('SET', 'LIN', PE3D.IdxCenter_PE - 1), ...
                     mr.makeLabel('SET', 'PAR', PE3D.IdxCenter_3D - 1), ...
                     mr.makeLabel('SET', 'NOISE', true), ...
                     Label.lblResetRefScan, Label.lblResetRefAndImaScan);
        seq.addBlock(mr.makeLabel('SET', 'NOISE', false), mr.makeLabel('SET', 'ONCE', 0)) ; % reset 'ONCE'
    end
end
