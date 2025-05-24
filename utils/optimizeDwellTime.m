function [BWPerPixel_opt, ADCDwell_opt, ADCDuration_opt] = optimizeDwellTime(nRO, BWPerPixel_target, sys)
    % Define the search range for ADCDwell
    minDwell = sys.adcRasterTime;
    maxDwell = 1e-3;  % Maximum dwell time to consider (adjustable)
    dwellCandidates = minDwell : sys.adcRasterTime : maxDwell;

    % Initialize best match
    bestError = inf;
    BWPerPixel_opt = NaN;
    ADCDwell_opt = NaN;
    ADCDuration_opt = NaN;

    for dwell = dwellCandidates
        % Compute the total ADC duration for current dwell time
        ADCDuration = nRO * dwell;

        % Check if ADCDuration is an integer multiple of gradRasterTime
        if mod(ADCDuration, sys.gradRasterTime) == 0
            % Compute the corresponding bandwidth per pixel
            BW_candidate = 1 / (nRO * dwell);

            % Compute the absolute error from the target bandwidth
            error = abs(BW_candidate - BWPerPixel_target);

            % Update optimal values if error is smaller
            if error < bestError
                bestError = error;
                BWPerPixel_opt = BW_candidate;
                ADCDwell_opt = dwell;
                ADCDuration_opt = ADCDuration;
            end
        end
    end

    % If no valid dwell time is found, throw an error
    if isnan(ADCDwell_opt)
        error('No valid ADCDwell found. Try increasing maxDwell or relaxing timing constraints.');
    end
end