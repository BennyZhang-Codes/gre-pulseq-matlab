function [BWPerPixel_opt, ADCDwell_opt, ADCDuration_opt] = optimizeADCDuration(nRO, BWPerPixel_target, sys)
    % Set the search parameters
    minDuration = nRO * sys.adcRasterTime;  % Minimum possible duration
    maxDuration = 100e-3;                   % Max duration to search (adjust as needed)
    stepSize = 10e-6;                       % Fixed step in seconds (10 µs)

    durationCandidates = minDuration : stepSize : maxDuration;

    % Initialize best match
    bestError = 100;
    BWPerPixel_opt = NaN;
    ADCDwell_opt = NaN;
    ADCDuration_opt = NaN;

    for duration = durationCandidates
        dwell = duration / nRO;

        % Check if dwell is a multiple of adcRasterTime
        if mod(dwell, sys.adcRasterTime) == 0
            BW_candidate = 1 / (nRO * dwell);
            error = abs(BW_candidate - BWPerPixel_target);

            if error < bestError
                bestError = error;
                BWPerPixel_opt = BW_candidate;
                ADCDwell_opt = dwell;
                ADCDuration_opt = duration;
            end
        end
    end

    if isnan(ADCDwell_opt)
        error('No valid ADCDuration found. Try increasing maxDuration or adjusting constraints.');
    end
end
