function [Grad] = prep_Spoiler(Actual, Grad, sys)
    % mapping of RO/PE/3D to X/Y/Z
    AxisRO = Actual.AxisRO ;
    AxisPE = Actual.AxisPE ;
    Axis3D = Actual.Axis3D ; % Slab/Paritions
    
    % Flip or not X/Y/Z to match patient positive/negative directions. Usefull
    % if reconstruction is done by system to get correct orientation of images.
    SignCorr = Actual.SignCorr ; % a structure


    % Similar to prephaser/rephasers we will first find the area in each
    % direction, then combine as root of sum of squares, then define a single
    % gradient to cover that and finally make copies of per direction with
    % matching amplitudes.
    
    % 1. RO spoiler area
    if (~isfield(Actual, 'SpoilerArea_RO'))      , Actual.SpoilerArea_RO       = [] ; end
    if (~isfield(Actual, 'SpoilerAreaFactor_RO')), Actual.SpoilerAreaFactor_RO = [] ; end
    % How is spoiler defined? Absolute area or relative area
    bSpoilerAreaDefinedRO   = ~isempty(Actual.SpoilerArea_RO)       ;
    bSpoilerFactorDefinedRO = ~isempty(Actual.SpoilerAreaFactor_RO) ;
    % Set area of spoiler
    if (bSpoilerAreaDefinedRO && ~bSpoilerFactorDefinedRO) 
        % Area explcitly defined in mT*us/m. Translate to Hz*s/m
        AreaSpoiler_RO = Actual.SpoilerArea_RO * 1e-3 * SystemDefault.gamma ; % [Hz*s/m]
    elseif (bSpoilerFactorDefinedRO && ~bSpoilerAreaDefinedRO )
        AreaSpoiler_RO = Actual.SpoilerAreaFactor_RO * Actual.nRO/Actual.FOV(1) ; % [Hz*s/m = 1/m]
    elseif (~bSpoilerAreaDefinedRO && ~bSpoilerFactorDefinedRO)
        warning('No RO spoiler is defined. Neither ''SpoilerArea_RO'' nor ''SpoilerAreaFactor_RO''. Assuming zero.');
        AreaSpoiler_RO = 0 ;
    else % both are defined
        error('RO spoiler is defined twice, both in ''SpoilerArea_RO'' and in ''SpoilerAreaFactor_RO''. Set one to be empty.');
    end
    
    % 2. PE spoiler area
    if (~isfield(Actual, 'SpoilerArea_PE'))      , Actual.SpoilerArea_PE       = [] ; end
    if (~isfield(Actual, 'SpoilerAreaFactor_PE')), Actual.SpoilerAreaFactor_PE = [] ; end
    % How is spoiler defined? Absolute area or relative area
    bSpoilerAreaDefinedPE   = ~isempty(Actual.SpoilerArea_PE)       ;
    bSpoilerFactorDefinedPE = ~isempty(Actual.SpoilerAreaFactor_PE) ;
    % Set area of spoiler
    if (bSpoilerAreaDefinedPE && ~bSpoilerFactorDefinedPE) 
        % Area explcitly defined in mT*us/m. Translate to Hz*s/m
        AreaSpoiler_PE = Actual.SpoilerArea_PE * 1e-3 * SystemDefault.gamma ; % [Hz*s/m]
    elseif (bSpoilerFactorDefinedPE && ~bSpoilerAreaDefinedPE )
        AreaSpoiler_PE = Actual.SpoilerAreaFactor_PE * Actual.nPE/Actual.FOV(2) ; % [Hz*s/m = 1/m]
    elseif (~bSpoilerAreaDefinedPE && ~bSpoilerFactorDefinedPE)
        warning('No PE spoiler is defined. Neither ''SpoilerArea_PE'' nor ''SpoilerAreaFactor_PE''. Assuming zero.');
        AreaSpoiler_PE = 0 ;
    else % both are defined
          error('PE spoiler is defined twice, both in ''SpoilerArea_PE'' and in ''SpoilerAreaFactor_PE''. Set one to be empty.');
    end
    
    % 3. 3D spoiler area
    if (~isfield(Actual, 'SpoilerArea_3D'))      , Actual.SpoilerArea_3D       = [] ; end
    if (~isfield(Actual, 'SpoilerAreaFactor_3D')), Actual.SpoilerAreaFactor_3D = [] ; end
    % How is spoiler defined? Absolute area or relative area
    bSpoilerAreaDefined3D   = ~isempty(Actual.SpoilerArea_3D)       ;
    bSpoilerFactorDefined3D = ~isempty(Actual.SpoilerAreaFactor_3D) ;
    % Set area of spoiler
    if (bSpoilerAreaDefined3D && ~bSpoilerFactorDefined3D) 
        % Area explcitly defined in mT*us/m. Translate to Hz*s/m
        AreaSpoiler_3D = Actual.SpoilerArea_3D * 1e-3 * SystemDefault.gamma ; % [Hz*s/m]
    elseif (bSpoilerFactorDefined3D && ~bSpoilerAreaDefined3D )
        AreaSpoiler_3D = Actual.SpoilerAreaFactor_3D * Actual.n3D/Actual.FOV(3) ; % [Hz*s/m = 1/m]
    elseif (~bSpoilerAreaDefined3D && ~bSpoilerFactorDefined3D)
        warning('No 3D spoiler is defined. Neither ''SpoilerArea_3D'' nor ''SpoilerAreaFactor_3D''. Assuming zero.');
        AreaSpoiler_3D = 0 ;
    else % both are defined
        error('3D spoiler is defined twice, both in ''SpoilerArea_3D'' and in ''SpoilerAreaFactor_3D''. Set one to be empty.');
    end
    
    
    % Single total spoiler event (supports any oblique)
    % -------------------------------------------------
    
    % We start by defining a gradient that supports the combined area
    % (root-sum-of-squares) of all directions. The AxisRO direction used is
    % arbitrary.
    SpoilMaxAbsArea = sqrt(AreaSpoiler_RO^2 + AreaSpoiler_PE^2 + AreaSpoiler_3D^2) ;
    if (SpoilMaxAbsArea == 0 )
        Grad.bSpoilers = false ; % Mark that no spoiler is used
    else
        Grad.bSpoilers = true  ; % Mark that spoilers are used
        
        % Define dummy "RO" spoiler event
        Grad.Spoil = mr.makeTrapezoid(AxisRO, sys, 'area', SpoilMaxAbsArea);
        SpoilDuration = Grad.Spoil.riseTime + Grad.Spoil.flatTime + Grad.Spoil.fallTime ;
        
        % RO spoiler event
        Grad.GRO_Spoil =  mr.makeTrapezoid(AxisRO, sys, 'duration', SpoilDuration, ...
        'riseTime', Grad.Spoil.riseTime, 'fallTime', Grad.Spoil.fallTime, 'area', SignCorr.(AxisRO) * AreaSpoiler_RO) ;
        
        % PE spoiler event
        Grad.GPE_Spoil =  mr.makeTrapezoid(AxisPE, sys, 'duration', SpoilDuration, ...
        'riseTime', Grad.Spoil.riseTime, 'fallTime', Grad.Spoil.fallTime, 'area', SignCorr.(AxisPE) * AreaSpoiler_PE) ;
        
        % 3D spoiler event
        Grad.G3D_Spoil =  mr.makeTrapezoid(Axis3D, sys, 'duration', SpoilDuration, ...
        'riseTime', Grad.Spoil.riseTime, 'fallTime', Grad.Spoil.fallTime, 'area', SignCorr.(Axis3D) * AreaSpoiler_3D) ;
    end

end
