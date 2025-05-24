function [Grad, ADC, Delay] = prep_ADC(Actual, Grad, ADC, Delay, sys)

    ADCDwell = round(1/(Actual.nRO*Actual.BWPerPixel)/sys.adcRasterTime)*sys.adcRasterTime ; % [s] rounded to adcRasterTime
    ADCDuration = Actual.nRO * ADCDwell ; % Get ADC duration and the duration rounded up to the gradient raster
    Actual.BWPerPixel = 1/ADCDuration ; % [Hz/pixel] Update actual BW per pixel used
    
    % Find the ADC duration rounded to the twice(!) gradient raster. (We want
    % to ensure the center of the flat area is on the gradient raster, so TEs
    % can fall on the gradient raster
    % NOTE: We implicitly assume below (with the ADCDelay) that the RF raster
    %       is finer than the gradient raster. (That we do not have to worry
    ADCDurGradRaster = ceil(ADCDuration / (2*sys.gradRasterTime)) * (2*sys.gradRasterTime) ; % [s]
    
    % define RO gradient event
    ROAmp = 1 / Actual.FOV(1) / ADCDwell; % [Hz/m]
    Grad.GRO_acq = mr.makeTrapezoid('x', sys, 'flatTime', ADCDurGradRaster, 'amplitude', ROAmp) ;
    % define ADC event inc. delay relative to start of RO gradient.
    % Note: delay is in RF raster time, not ADC raster time.
    ADCDelay = Grad.GRO_acq.riseTime + (ADCDurGradRaster - ADCDuration)/2 ; % [s]
    ADCDelay = ceil(ADCDelay/sys.rfRasterTime) * sys.rfRasterTime ; % Ensure on RF(!) raster.
    ADC.adc  = mr.makeAdc(Actual.nRO, 'dwell', ADCDwell , 'delay', ADCDelay, 'system', sys) ; 
    
    % Round up total time of adc to the block raster time. In rare cases
    % the total time may be longer than the RO gradient and not be on the block
    % raster time. This will be used to slightly extend the block duration to be on the block raster.
    ROADCGrossDurBlockRaster = ceil(mr.calcDuration(ADC.adc) / sys.blockDurationRaster) * sys.blockDurationRaster ; % [s]
    % Create a matching "delay" (played in parallel to the other events in the block, so will just ensure its duration.
    Delay.ROADCDelay = mr.makeDelay(ROADCGrossDurBlockRaster) ;


    Delay.ADCDelay = ADCDelay;
    Delay.ADCDuration = ADCDuration;
end
