function [RF, Grad] = prep_RF(Actual, RF, Grad, sys)

    % mapping of RO/PE/3D to X/Y/Z
    % AxisRO   = Actual.AxisRO   ;
    % AxisPE   = Actual.AxisPE   ;
    Axis3D   = Actual.Axis3D   ; 
    SignCorr = Actual.SignCorr ; 

    fa = Actual.flipEx * pi/180 ;
    if (isempty(Actual.SlabThickness)) % non-selective excitation
        [RF.rf_ex] = mr.makeBlockPulse(fa, sys, 'duration', Actual.tEx_All, 'use', 'excitation') ;
    else                               % slab selective excitation
        [RF.rf_ex, Grad.G3D_Ex, Grad.G3D_Ex_Ref] = mr.makeSincPulse(fa, sys, 'duration', Actual.tEx_Slab, ...
            'timeBwProduct', Actual.tbpEx_Slab, 'apodization', 0.5, 'sliceThickness', Actual.SlabThickness, 'use', 'excitation') ;
        % Fix direction & sign (of two gradients)
        Grad.G3D_Ex.channel           = Axis3D                                        ;
        Grad.G3D_Ex.amplitude         = SignCorr.(Axis3D) * Grad.G3D_Ex.amplitude     ;
        Grad.G3D_Ex.area              = SignCorr.(Axis3D) * Grad.G3D_Ex.area          ;
        Grad.G3D_Ex_Ref.channel       = Axis3D                                        ;
        Grad.G3D_Ex_Ref.amplitude     = SignCorr.(Axis3D) * Grad.G3D_Ex_Ref.amplitude ;
        Grad.G3D_Ex_Ref.area          = SignCorr.(Axis3D) * Grad.G3D_Ex_Ref.area      ;
        % Update Actual parameters
        Actual.RFExciteSlabDuration   = RF.rf_ex.shape_dur                            ;
        Actual.RFExciteSlabNumSamples = numel(RF.rf_ex.t)                             ;
    end
    Grad.G3D_ExAndRef = concatGrads({Grad.G3D_Ex, Grad.G3D_Ex_Ref}, sys);
end
