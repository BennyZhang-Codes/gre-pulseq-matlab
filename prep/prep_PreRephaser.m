function [Grad] = prep_PreRephaser(Actual, Grad, PE3D, sys)
    % mapping of RO/PE/3D to X/Y/Z
    AxisRO = Actual.AxisRO ;
    AxisPE = Actual.AxisPE ;
    Axis3D = Actual.Axis3D ; % Slab/Paritions
    
    % Flip or not X/Y/Z to match patient positive/negative directions. Usefull
    % if reconstruction is done by system to get correct orientation of images.
    SignCorr = Actual.SignCorr ; % a structure


    % Max. RO prephase
    GradPrephaseMaxAbsArea_RO = abs(Grad.GRO_acq.area/2) ; % [1/m] ; % Max absolute area of the RO prephase(!) gradient
    
    % Max. PE prehase
    GradMaxAbsArea_PE = PE3D.IdxAbsMax_PE/Actual.FOV(2)       ; % [1/m] % Max absolute area of the prephase gradient
    
    % Max. 3D prehase
    GradMaxAbsArea_3D = PE3D.IdxAbsMax_3D/Actual.FOV(3)       ; % [1/m] % Max absolute area of the prephase gradient
    
    % -------------------------------------------------------------------------
    % Single total prephaser event (supports any oblique)
    % -------------------------------------------------------------------------
    
    % Our prephasers (RO, PE, and 3D) are going to have the same timings (ramp
    % up, flat duration, and ramp down) and will also supprt any oblique
    % direction. For this we'll first define a trapazoid with an area that
    % combines all the prephasers (combines areas via a root of the sum of
    % squares). Once we have that we'll defined three separate gradients with
    % the same timing with the respective area (or max area).
    
    % We start by defining a gradient that supports the combined area
    % (root-sum-of-squares) of all directions. The AxisRO direction used is
    % arbitrary.
    PrephaseMaxAbsArea = sqrt(GradPrephaseMaxAbsArea_RO^2 + GradMaxAbsArea_PE^2 + GradMaxAbsArea_3D^2) ;
    Grad.Prephase = mr.makeTrapezoid(AxisRO, sys, 'area', PrephaseMaxAbsArea) ;
    
    % -------------------------------------------------------------------------
    % PE prephaser/rephaser
    % -------------------------------------------------------------------------
    % Define PE prephase (Amplitude Will be updated later)
    Grad.GPE_Prephase         = Grad.Prephase ;
    Grad.GPE_Prephase.channel = AxisPE        ;
    
    % Define PE rephaser. (Amplitude Will be updated later according to actual GPE_Prephase used)
    Grad.GPE_Rephase = Grad.GPE_Prephase ;
    
    % The amplitude of the prephase gradient for 1 k-space step, for the
    % gradient just defined. (We will use an amplitude which is an integer multiple of this.)
    Grad.PrePhaseAmpStepPE = SignCorr.(AxisPE) * (GradMaxAbsArea_PE/PrephaseMaxAbsArea * ...
        Grad.GPE_Prephase.amplitude) / PE3D.IdxAbsMax_PE ; % [Hz/m]
    
    % -------------------------------------------------------------------------
    % 3D prephaser/rephaser
    % -------------------------------------------------------------------------
    
    % Define 3D prephase (Amplitude Will be updated later)
    Grad.G3D_Prephase         = Grad.Prephase ;
    Grad.G3D_Prephase.channel = Axis3D        ;
    
    % Define 3D rephaser. (Amplitude Will be updated later according to actual
    % Grad.G3D_Prephase used)
    Grad.G3D_Rephase = Grad.G3D_Prephase ;
    
    % The amplitude of the prephase gradient for 1 k-space step, for the
    % gradient just defined. (We will use an amplitude which is an integer
    % multiple of this.)
    Grad.PrePhaseAmpStep3D = SignCorr.(Axis3D) * (GradMaxAbsArea_3D/PrephaseMaxAbsArea * ...
        Grad.G3D_Prephase.amplitude) / PE3D.IdxAbsMax_3D ; % [Hz/m]
    
    % -------------------------------------------------------------------------
    % RO prephaser/rephaser (and TE filler delays)
    % -------------------------------------------------------------------------
    
    % Define RO prephaser event
    % Its area should be half the area of the RO gradient (we shifted the ADC to be centered at the center of the RO gradient).
    Grad.GRO_Prephase         = Grad.Prephase ;
    Grad.GRO_Prephase.channel = AxisRO        ;
    % Set ampilitude. Because it depends on GROacq, whose sign has already been corrected, we do not need to correct sign here.
    Grad.GRO_Prephase.amplitude = (-Grad.GRO_acq.area/2) * Grad.GRO_Prephase.amplitude/PrephaseMaxAbsArea ;
    
    % Define RO rephaser. 
    % (Same as prepahser! the sum of prephaser and rephaser should cancel the RO so they are both half the RO with reversed sign.)
    Grad.GRO_Rephase = Grad.GRO_Prephase ;
    
    % Define RO prephaser for multiple mono-polar TEs (bBipolarROGrads is
    % false). In this case after the acquisition RO gradient we have to negate
    % it before we can run it again for the next TE.
    Grad.GRO_acqUndo = mr.makeTrapezoid(AxisRO, sys, 'area', -Grad.GRO_acq.area) ;

end
