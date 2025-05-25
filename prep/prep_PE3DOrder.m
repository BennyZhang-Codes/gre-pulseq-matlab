function [PE3D] = prep_PE3DOrder(Actual)
    %% Helper (anonymous) functions
    
    % Checks if a number is practically an integer:
    IsIntValue = @(n) abs(rem(n,1)) <= eps(1) ;

    IdxCenter_PE      = Ordering.Utils.FFTCenterIndex(Actual.nPE) ; 
    IdxMax_PE         = Actual.nPE - IdxCenter_PE        ; % Find start and end indices of y phase encoding.
    IdxMin_PE         = IdxMax_PE  - Actual.nPE + 1      ;
    IdxAbsMax_PE      = max(abs([IdxMax_PE, IdxMin_PE])) ; % Find largest index of the two, in absolute value.

    IdxCenter_3D      = Ordering.Utils.FFTCenterIndex(Actual.n3D) ; 
    IdxMax_3D         = Actual.n3D - IdxCenter_3D        ; % Find start and end indices of z phase encoding.
    IdxMin_3D         = IdxMax_3D  - Actual.n3D + 1      ;
    IdxAbsMax_3D      = max(abs([IdxMax_3D, IdxMin_3D])) ; % Find largest index of the two, in absolute value.

    % Is there (valid) acceleration along PE or 3D
    bAccelerationPE = false ; % intialize no acceleration along PE
    bAcceleration3D = false ; % intialize no acceleration along 3D
    if (isfield(Actual, 'AccelerationPE') && ~isempty(Actual.AccelerationPE))
        if (~IsIntValue(Actual.AccelerationPE) || Actual.AccelerationPE < 1)
            error('AccelerationPE is expected to be a positive integer.')
        end
        if (Actual.AccelerationPE > 1)
            bAccelerationPE = true ;
        end
    end
    if (isfield(Actual, 'Acceleration3D') && ~isempty(Actual.Acceleration3D))
        if (~IsIntValue(Actual.Acceleration3D) || Actual.Acceleration3D < 1)
            error('Acceleration3D is expected to be a positive integer.')
        end
        if (Actual.Acceleration3D > 1)
            bAcceleration3D = true ;
        end
    end
    
    % Switch from PE and 3D to fast (#1) and slow (#2) dimensions. (Which is fast/slow depends on 'DimFast'.)
    Acceleration1 = 1 ;  % initialize
    nRefLine1     = [];
    Acceleration2 = 1 ;
    nRefLine2     = [];
    
    switch lower(Actual.DimFast)
        case 'pe' % Step through PE first (and then 3D)
            % stepped through first
            n1                = Actual.nPE            ; 
            IdxMin_1          = IdxMin_PE             ;
            IdxMax_1          = IdxMax_PE             ;
            IdxCenter_1       = IdxCenter_PE          ;
            bAcceleration1    = bAccelerationPE       ;
            if (bAcceleration1)
                Acceleration1 = Actual.AccelerationPE ;
                nRefLine1     = Actual.nRefLinePE     ;
            end
        
            % stepped through second
            n2                = Actual.n3D            ; 
            IdxMin_2          = IdxMin_3D             ;
            IdxMax_2          = IdxMax_3D             ;
            IdxCenter_2       = IdxCenter_3D          ;
            bAcceleration2    = bAcceleration3D       ;
            if (bAcceleration2)
                Acceleration2 = Actual.Acceleration3D ;
                nRefLine2     = Actual.nRefLine3D     ;
            end
        case '3d'  % Step through 3D first (and then PE)
            % stepped through first
            n1                = Actual.n3D            ; 
            IdxMin_1          = IdxMin_3D             ;
            IdxMax_1          = IdxMax_3D             ;
            IdxCenter_1       = IdxCenter_3D          ;
            bAcceleration1    = bAcceleration3D       ;
            if (bAcceleration1)
                Acceleration1 = Actual.Acceleration3D ;
                nRefLine1     = Actual.nRefLine3D     ;
            end
            
            % stepped through second
            n2                = Actual.nPE            ; 
            IdxMin_2          = IdxMin_PE             ;
            IdxMax_2          = IdxMax_PE             ;
            IdxCenter_2       = IdxCenter_PE          ;
            bAcceleration2    = bAccelerationPE       ;
            if (bAcceleration2)
                Acceleration2 = Actual.AccelerationPE ;
                nRefLine2     = Actual.nRefLinePE     ;
            end 
        otherwise
            error('DimFast is ''%s'' but should be either ''PE'' or ''3D'' ', Actual.DimFast);
    end
    
    % Set Sampling masks 
    if (bAcceleration1 || bAcceleration2)
        % [bSample, bRefFull, bImaAndRefFull] = Ordering.Utils.GrappaMasks(n1, n2, Acceleration1, Acceleration2, nRefLine1, nRefLine2) ;
        [bSample, bRefFull, bImaAndRefFull] = CAIPIRINHAMasks(n1, n2, Acceleration1, Acceleration2, nRefLine1, nRefLine2, Actual.CAIPIShift) ;
    else
        [bSample, bRefFull, bImaAndRefFull] = Ordering.Utils.GrappaMasks(n1, n2, 1, 1) ;
    end

    %%
    % Set order of sampling
    % Set handle fOrdering to get the desired order. NOTE: IOut and JOut are zero for k=0.
    
    % product of cutoff frequency and TR
    TRtimesCutoffFreq = Actual.CutoffFreq * Actual.TR ;
    
    [orderMat, IOut, JOut, SampledReorder] = Actual.fOrdering(n1, n2, ...
        TRtimesCutoffFreq, bSample, Actual.OrderingExtraParamsStruct, Actual.RandomSeed) ;
    
    % Extract subset of Ref/ImaAndRef from actually sampled samples to
    % match the length of IOut and Jout. Recall that IOut and JOut are possibly
    % reorderd relative to original bRefFull and bImaAndRefFull, so we have
    % to reorder them and cut them short. This is what we have orderOut for.
    bRef       = bRefFull(bSample(:))          ;
    bRef       = bRef(SampledReorder(:))       ;
    bImaAndRef = bImaAndRefFull(bSample(:))    ;
    bImaAndRef = bImaAndRef(SampledReorder(:)) ;
    
    % conversion to k-space representation format (indices are zero for k=0).
    IOut = IOut - Ordering.Utils.FFTCenterIndex(n1);
    JOut = JOut - Ordering.Utils.FFTCenterIndex(n2);
    
    % Get list of PE and 3D indices to sample: return from using fast dimension (1) and slow dimension (2).
    if strcmpi(Actual.DimFast, 'PE') % Step through Y (#1) first and then Z (#2)
        PE3DOrder = [IOut(:), JOut(:)] ;
    else
        PE3DOrder = [JOut(:), IOut(:)] ;
    end
    
    % % DEBUG: Show bRef and bImaAndRef are marked correctly.



    PE3D.bSample        = bSample       ;
    PE3D.bRefFull       = bRefFull      ;
    PE3D.bImaAndRefFull = bImaAndRefFull;

    PE3D.IdxCenter_PE   = IdxCenter_PE  ; 
    PE3D.IdxMax_PE      = IdxMax_PE     ;
    PE3D.IdxMin_PE      = IdxMin_PE     ;
    PE3D.IdxAbsMax_PE   = IdxAbsMax_PE  ; 

    PE3D.IdxCenter_3D   = IdxCenter_3D  ; 
    PE3D.IdxMax_3D      = IdxMax_3D     ;
    PE3D.IdxMin_3D      = IdxMin_3D     ;
    PE3D.IdxAbsMax_3D   = IdxAbsMax_3D  ;

    PE3D.bRef           = bRef          ;
    PE3D.bImaAndRef     = bImaAndRef    ;
    PE3D.PE3DOrder      = PE3DOrder     ;
end
