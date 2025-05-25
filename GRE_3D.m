% % clc; clear; close all;
%% Create a TSE sequence
addpath(genpath('pulseq'));
addpath(genpath('prep'  ));
addpath(genpath('check' ));
addpath(genpath('plot'  ));
addpath(genpath('utils' ));
% Instantiation and gradient limits
sys = mr.opts('MaxGrad', 40, 'GradUnit', 'mT/m', ...
    'MaxSlew', 200, 'SlewUnit', 'T/m/s', ...
    'rfRingdownTime', 100e-6, 'rfDeadTime', 100e-6, ...
    'adcDeadTime', 10e-6, 'adcRasterTime',100e-9,...
    'rfRasterTime',1e-6,'gradRasterTime',10e-6, ...
    'B0', 6.98, 'adcSamplesDivisor', 4);
seq = mr.Sequence(sys);
RF    = struct();
Grad  = struct();
ADC   = struct();
Delay = struct();
Label = struct();
%% Sequence parameters
params.NoiseScan        = 'on'; 
params.nDummy           = 5; 
params.nRep             = 1;

% 1
params.fovRO            = 200e-3;   %
params.fovPE            = 200e-3;
params.fov3D            = 160e-3;
params.FOV              = [params.fovRO, params.fovPE, params.fov3D]; % [m] RO x PE x 3D
params.SlabThickness    = params.FOV(3); % [m] empty mean no slab selection, i.e. excite all

params.nRO              = 100; 
params.nPE              = 100; 
params.n3D              = 80;
params.MatrixSize       = [params.nRO, params.nPE, params.n3D];       % [a.u.] RO x PE x 3D

params.AccelerationPE   = 3  ; % acceleration factor for phase direction.
params.Acceleration3D   = 3  ; % acceleration factor for partition direction.
params.nRefLinePE       = 24 ; % Number of fully sampled lines at center, along PE
params.nRefLine3D       = 24 ; % Number of fully sampled lines at center, along 3D
params.CAIPIShift       = 2  ;
params.DimFast          = 'PE' ; % 'PE' or '3D'; which dimension is stepped through first

% Array of TEs. Non-positive values will be interperted as using minimum TE possible.
params.TE               = [2.1, 3.6, 5.1]*1e-3; % non-positive values mean use shortest TE.
% 2.1, 3.6, 5.1   % 2.10, 3.12, 4.14 % -1, -1, -1
params.TR               = 25e-3; % [s]
params.BWPerPixel       = 1250 ;       % [Hz/pixel] updated later to obey ADC raster

params.nTE              = numel(params.TE);


% When acquiring with multiple TEs, should the RO gradient alternate in
% sign (faster) or always have the same sign (slower). Switching sign is
% refered to as bi-polar
params.bBipolarROGrads  = false ;


% RF pulse
% Define slab-selective excitation (in case it is requested)
% Original Definitions are based on Siemens a_gre (VE12U)
params.flipEx           = 10   ; % [deg]
params.tEx_Slab         = 2.0e-3 ; % [s]
params.tbpEx_Slab       = 8   ; % [a.u.]

% Define non-selective excitation 
% Original Definitions are based on Siemens a_gre (VE12U)
params.tEx_All          = 0.1e-3 ; % [s]


% RF & Gradient Spoiling
params.RFSpoilIncDeg    = 117 ; % [deg] RF spoiling increment

% Spoiler area can be defined either explicitly in mT*us/m or relative to params
% the net area of gradients required to achieve desired resolution. Two
% sets of variables are given at least one must be empty.
params.SpoilerArea_RO       = [] ; % [mT*us/m]
params.SpoilerArea_PE       = [] ; % [mT*us/m]
params.SpoilerArea_3D       = [] ; % [mT*us/m]
params.SpoilerAreaFactor_RO = 2 ;
params.SpoilerAreaFactor_PE = 2 ;
params.SpoilerAreaFactor_3D = 2 ;


% ---------
% Shuffling
% ---------
% Cutoff frequency (if relevant): 
% - [] (empty) - no shuffling.
% - 0 - full randomization.
% - > 0 - Cutoff to use for scrambling (if relevant)
params.CutoffFreq = [] ; 1/20 ; % [Hz] 

% Set trajectory
% set function handle to function used to order the trajectory.
% All functions accept the same inputs, although might not use them all.
params.fOrdering = @Ordering.Ordered_LocalShuffle ;
% params.fOrdering = @Ordering.SpiralRS_LocalShuffle ;
% params.fOrdering = @Ordering.Ordered_SegmentedShuffle ;
% params.fOrdering = @Ordering.SpiralSquare_LocalShuffle ;
% params.fOrdering = @Ordering.Gilbert ;
% params.fOrdering = @Ordering.Spiral_LocalShuffle ;
% Optional structure of extra parameters which are specific to the ordering
% function used.
params.OrderingExtraParamsStruct = [] ;
% params.OrderingExtraParamsStruct.bElliptic = true ;

% Random seed for shuffling
params.RandomSeed = 0 ;

%% Set orienation (non-oblique)
% Set axes (X/Y/Z vs. RO/PE/3D - oblique not supported) 
% In Siemens interpreter the defintions here must agree with the
% 'Orientation mapping' setting.

% mapping of RO/PE/3D to X/Y/Z
params.AxisRO = 'x' ;
params.AxisPE = 'y' ;
params.Axis3D = 'z' ; 

% Flip or not X/Y/Z to match patient positive/negative directions. Usefull
% if reconstruction is done by system to get correct orientation of images.
params.SignCorr.x = -1 ;
params.SignCorr.y = -1 ;
params.SignCorr.z = -1 ;

%%
Actual = params;
[Actual] = prep_TRTE(Actual, params, sys);

%% RF
[RF, Grad] = prep_RF(Actual, RF, Grad, sys);

%% RO: ADC and gradient
[Grad, ADC, Delay] = prep_ADC(Actual, Grad, ADC, Delay, sys);

%% Set ADC labels (Siemens MDH flags)
[Label] = prep_Label(Actual, Label, seq);

%% Set PE & 3D ordering
Actual.DimFast          = 'PE';
Actual.CutoffFreq       = [];
[PE3D] = prep_PE3DOrder(Actual);
[fig] = plot_PE3D(Actual, PE3D);
[fig] = plot_PE3DOrder(Actual, PE3D);

%% Prephaser & Rephaser
[Grad] = prep_PreRephaser(Actual, Grad, PE3D, sys);

%% Spoiler(s)
[Grad] = prep_Spoiler(Actual, Grad, sys);

% TR filler delay

% Initialize a zero delay that will be modified later on
% to ensure we get the desired TR (if possible).
Delay.delayTRFill = mr.makeDelay(1) ; % dummy time. (Delay of zero is not allowed)
Delay.delayTRFill.delay = 0 ; % [s] force delay of zero.

%% Noise Scan
[seq, Label] = prep_NoiseScan(seq, Actual, PE3D, ADC, Label, sys);

%%
% Initialize variables
RFSpoilInc = Actual.RFSpoilIncDeg * pi/180 ; % [rad]
% What is the actual minimum TR possible for current sequence (found when
% generating the sequence)
MinTRActual = 0 ;

% translate FOV from RO/PE/3D to X/Y/Z:
[~, FOVXYZOrder] = sort([Actual.AxisRO, Actual.AxisPE, Actual.Axis3D], 'ascend');
Actual.FOV       = Actual.FOV(FOVXYZOrder)                 ;

% Initialize sequence
RFSpoilPhase      = RF.rf_ex.phaseOffset; % [rad]
RFSpoilPhaseShift = 0                   ; % [rad] Next step of RF phase (updates)

% Store RO amplitude, so we can reset to this value at the start of each
% slice. (In case of a multi-echo bi-polar acquisition, the sign will
% alternate, so we wan to ensure we always start at the same sign.)
Grad.ROAmp0 = Grad.GRO_acq.amplitude ;
%%
seq.addBlock(mr.makeLabel('SET', 'REP', 0));
for irep = 1:Actual.nRep
% Loop over TR: non-positive counter values mean dummy scan (ky=kz=0). 
% Otherwise, advance in the PE3DOrder table of phase encodes (PE & 3D) 
for TRCounter = (-Actual.nDummy+1):size(PE3D.PE3DOrder, 1)

    % Reset duration of current TR
    TimeInTR = 0 ; % [s]

    % Set y (PE) and z (3D) indices
    if (TRCounter <= 0) % dummy scans)
        Idx_PE = 0                           ; % k = 0 for PE
        Idx_3D = 0                           ; % k = 0 for 3D
        bADCOn = false                       ; 
    else
        Idx_PE = PE3D.PE3DOrder(TRCounter, 1);
        Idx_3D = PE3D.PE3DOrder(TRCounter, 2);
        bADCOn = true                        ; 
    end
    
    % Update RF spoiling phase to use in this TR
    % update RF phase to use this round
    RFSpoilPhase = RFSpoilPhase + RFSpoilPhaseShift ; % [rad]
    
    % upadte RFSpoilPhaseShift for next round ;
    RFSpoilPhaseShift = RFSpoilPhaseShift + RFSpoilInc ; % [rad]
    
    % Excitation block
    % update RF spoiling phase of RF pulse
    RF.rf_ex.phaseOffset = RFSpoilPhase ; % [rad]
    % Add excitation block
    if (isempty(Actual.SlabThickness)) % non-selective excitation
        seq.addBlock(RF.rf_ex) ;
    else % slab selective excitation
        seq.addBlock(RF.rf_ex, Grad.G3D_ExAndRef) ;
    end

    % Update duration within TR
    TimeInTR = TimeInTR + seq.blockDurations(end) ;
    
    % Initialize time from TE: time between RF center (from start of block)
    % to end of the block.
    TimeFromExcite = seq.blockDurations(end) - (RF.rf_ex.delay + RF.rf_ex.shape_dur/2) ; % [s]


    % -----------------------------------------------------------------------
    % Slab refocusing(?) + prephasing block
    % -----------------------------------------------------------------------
    
    % Update y prephaser amplitude
    Grad.GPE_Prephase.amplitude =  Idx_PE * Grad.PrePhaseAmpStepPE ; % [Hz/m]
    % Update z prephaser amplitude
    Grad.G3D_Prephase.amplitude =  Idx_3D * Grad.PrePhaseAmpStep3D ; % [Hz/m]
    
    % Add block of prephasers
    seq.addBlock(Grad.GRO_Prephase, Grad.GPE_Prephase, Grad.G3D_Prephase) ;
    
    
    % Update duration within TR
    TimeInTR = TimeInTR + seq.blockDurations(end) ;
    % Update time from "excitation"
    TimeFromExcite = TimeFromExcite + seq.blockDurations(end) ;
    
    % -----------------------------------------------------------------------
    % multi-echo RO acquisition
    % -----------------------------------------------------------------------
    
    % Set phase of ADC to match RF spoiling phase. This way the RF spoiling
    % phase will not affect the acquired signal
    ADC.adc.phaseOffset = RF.rf_ex.phaseOffset ;

    
    % Reset RO amplitude.
    % (In case we have multi-echos with bi-polar gradients.)
    Grad.GRO_acq.amplitude  = Grad.ROAmp0 ;
    
    % Reset label used to mark if we should reverse the ADC.
    % (In bi-polar multi-echo, each time we flip the sign of the RO gradient). 
    LabelReverse = Label.lblResetRev ; % Don't reverse the first time

    for iTE = 1:Actual.nTE
        % Prepare for bi-polar or mono-polar RO grads in multi TE
        
        % For TEs beyond the first we have to either switch the RO direction
        % alternatingly, or insert a rephaser before we can re-use our RO gradient
        if (iTE > 1)
            if (Actual.bBipolarROGrads)
                % We are in bi-polar multi TE mode, so switch sign of RO gradient
                Grad.GRO_acq.amplitude = (-1)^(iTE-1) * Grad.ROAmp0 ;
                % Switch between reversing ADC or not.
                if (mod(iTE, 2) == 1) % odd (count from 1)
                    LabelReverse = Label.lblResetRev ;
                else % Even echo
                    LabelReverse = Label.lblSetRev   ;
                end
            
            else % monopolar case
                % We have to fully undo the last RO gradient before we can run the next.
                seq.addBlock(Grad.GRO_acqUndo) ;
                
                % Update duration within TR
                TimeInTR = TimeInTR + seq.blockDurations(end) ;
                % Update time from "excitation"
                TimeFromExcite = TimeFromExcite + seq.blockDurations(end) ;
            end
        end
        
        % Insert TE filler delay and the acqusition 
        % ------------------------------------------
        
        % Set filler delay to achieve requested TE (rounded up later)
        if (Actual.TE(iTE) < 0)
            % Negative TE is interperted as using the minimal TE possible.
            TEFill = 0 ;
        else % try and achieve requested TE
            TEFill = Actual.TE(iTE) - (TimeFromExcite + Delay.ADCDelay + Delay.ADCDuration/2) ;
        end
    
        % Sanity check
        if (TEFill < 0)
            error(['Cannot achieve desired TE[%d] = %f ms. ' 'Minimum possible is %f ms.'], ...
                iTE, 1e3*Actual.TE(iTE), 1e3*(Actual.TE(iTE) - TEFill))
        else
            TEFill = round(TEFill/sys.gradRasterTime)*sys.gradRasterTime ; % Round to gradient raster
        end

        % Set ADC labels (PE and partition). Note that the first index of each
        % label is zero (instead of marking k=0 position as zero index).
        % NOTE: We set the labels explicitly (using 'SET'), because the order
        %       may be arbitrary (depending on the order within PE3DOrder).
        Label_PE = mr.makeLabel('SET', 'LIN', Idx_PE - PE3D.IdxMin_PE) ; % PE
        Label_3D = mr.makeLabel('SET', 'PAR', Idx_3D - PE3D.IdxMin_3D) ; % 3D
    
        LabelRefUse       = Label.lblResetRefScan       ;
        LabelImaAndRefUse = Label.lblResetRefAndImaScan ;
        % TRCounter > 0 only when bADCOn, so we add that to our test
        if bADCOn  && PE3D.bRef(TRCounter) % parallel imaging reference line
            % Reference line, so mark it as such.
            LabelRefUse = Label.lblSetRefScan ; 
            if PE3D.bImaAndRef(TRCounter) % Also a regular imaging line
                % reference and(!) imaging line, so mark it as such.
                LabelImaAndRefUse = Label.lblSetRefAndImaScan ;
            end
        end
    
    
        % Add delay to events (and remove it after adding to block)
        Grad.GRO_acq.delay     = Grad.GRO_acq.delay     + TEFill ;
        ADC.adc.delay          = ADC.adc.delay          + TEFill ;
        Delay.ROADCDelay.delay = Delay.ROADCDelay.delay + TEFill ;
    
        if (bADCOn) % ADC used
            seq.addBlock(Grad.GRO_acq, ADC.adc, ...
                       Delay.ROADCDelay, ... ensure we are on the block raster
                       Label_PE, Label_3D, Label.lblSetEchos(iTE), ...
                       LabelRefUse, LabelImaAndRefUse, ...
                       LabelReverse) ; 
        else % no ADC
            seq.addBlock(Grad.GRO_acq, Delay.ROADCDelay) ; % ensure consistancy + the block raster
        end
    
        % remove extra delay (for next round)
        Grad.GRO_acq.delay     = Grad.GRO_acq.delay     - TEFill ;
        ADC.adc.delay          = ADC.adc.delay          - TEFill ;
        Delay.ROADCDelay.delay = Delay.ROADCDelay.delay - TEFill ;

        % update actual TE
        Actual.TE(iTE) = TimeFromExcite + TEFill + Delay.ADCDelay + Delay.ADCDuration/2 ;
        
        % Update duration within TR
        TimeInTR = TimeInTR + seq.blockDurations(end) ;
        % Update time from "excitation"
        TimeFromExcite = TimeFromExcite + seq.blockDurations(end) ;
    end

  
    % -----------------------------------------------------------------------
    % Rephasers
    % -----------------------------------------------------------------------

    % Update PE rephaser amplitude according to prephaser
    Grad.GPE_Rephase.amplitude =  -Grad.GPE_Prephase.amplitude ; % [Hz/m]
    % Update 3D rephaser amplitude according to prephaser
    Grad.G3D_Rephase.amplitude =  -Grad.G3D_Prephase.amplitude ; % [Hz/m]
    % Update RO rephaser amplitude according to prephaser
    if (Actual.bBipolarROGrads) % bi-polar RO gradients (alternating sign)
        Grad.GRO_Rephase.amplitude =  (-1)^(Actual.nTE+1) * Grad.GRO_Prephase.amplitude ; % [Hz/m]
    else % mono-polar
        % Same as prephaser, because their sum should cancel the RO gradient.
        Grad.GRO_Rephase.amplitude =  Grad.GRO_Prephase.amplitude ; % [Hz/m]
    end

    % Add block of rephasers. 
    seq.addBlock(Grad.GRO_Rephase, Grad.GPE_Rephase, Grad.G3D_Rephase) ;

    % Update duration within TR
    TimeInTR = TimeInTR + seq.blockDurations(end) ;

    % -----------------------------------------------------------------------
    % Spoilers
    % -----------------------------------------------------------------------
    if (Grad.bSpoilers)  
        seq.addBlock(Grad.GRO_Spoil, Grad.GPE_Spoil, Grad.G3D_Spoil) ; % Add block of spoilers. 
        TimeInTR = TimeInTR + seq.blockDurations(end) ;                % Update duration within TR
    end

    % -----------------------------------------------------------------------
    % Update minimum possible TR (for information only)
    % -----------------------------------------------------------------------
    % update minimum possible TR (before adding the TR fill time)
    MinTRActual = max(MinTRActual, TimeInTR) ;
    
    % -----------------------------------------------------------------------
    % TR Fill block
    % -----------------------------------------------------------------------
    % Set filler delay to achieve requested TR (rounded up later)
    TRFill = Actual.TR - TimeInTR ;
    
    % Sanity check
    if (TRFill < -eps(0))
        error(['Total time (%f ms) of blocks within current TR (#%d) is ' ...
        'longer than desired TR (%f ms)!'], 1e3*TimeInTR, TRCounter, 1e3*Actual.TR) ;
    end
    
    Delay.delayTRFill.delay = TRFill ; % update delay of eTRFill
    seq.addBlock(Delay.delayTRFill)  ;  % Add delay to the sequence
    
    % Update duration within TR
    TimeInTR = TimeInTR + seq.blockDurations(end) ;
end
seq.addBlock(mr.makeLabel('INC', 'REP', 1));
end

%% timing & PNS & definition
[seq] = check_Timing(seq);
% [seq] = check_PNS(seq);
[seq, str_res, str_mat, str_r] = prep_Definition(seq, Actual, PE3D);

outpath = 'E:/pulseq/idea/pulseq_150/GRE/';

seqname = sprintf('GRE_%s_%s_%s_tr%s_%ste_Mono', str_res, str_mat, str_r, ...
    num2str(Actual.TR*1e3), num2str(Actual.nTE));
seq.write(strcat(outpath, seqname,'.seq'));
save(strcat(outpath, seqname),'params','Actual', 'PE3D');
% fig = seq.plot('Label', 'LIN,SLC,ECO,REP');
% print(fig.f, '-dpng', '-loose', '-r300', '-image', sprintf('%s_allTR.png', seqname));

fig = seq.plot('Label', 'LIN,PAR,ECO,REP', 'timeRange', [0, 7*Actual.TR]);
% print(fig.f, '-dpng', '-loose', '-r300', '-image', sprintf('%s_1TR.png', seqname));
%% k-space trajectory calculation
% [ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();
% 
% % fig = plot_kspace(ktraj, ktraj_adc);
% figure; 
% plot(t_ktraj.', ktraj.') ;
% hold on ;
% xlabel('t [s]') ;
% ylabel('k [1/m]') ;
% title('k-space components as functions of time') ;
% % calculateKspacePP() should return physical x, y, z and not logical RO
% % PE and 3D.
% legend('k_x', 'k_y', 'k_z') ;
% 
% % Plot k-space trajectory (3D)
% % ---------------------------
% figure; 
% plot3(ktraj(1,:), ktraj(2,:), ktraj(3,:),'b') ;
% hold on ; 
% plot3(ktraj_adc(1,:), ktraj_adc(2,:), ktraj_adc(3,:), 'r.') ; 
% 
% % calculateKspacePP() should return physical x, y, z and not logical RO
% % PE and 3D.
% xlabel('k_x [1/m]') ;
% ylabel('k_y [1/m]') ;
% zlabel('k_z [1/m]') ;
% grid on ;
% title('3D k-space') ;
% legend('trajectory', 'ADC')
% 
%% evaluate label settings more specifically

lbls=seq.evalLabels('evolution','adc');
lbl_names=fieldnames(lbls);
figure; hold on;
for n=1:length(lbl_names)
    plot(lbls.(lbl_names{n}));
end
legend(lbl_names(:));
title('evolution of labels/counters/flags');
xlabel('adc number');



