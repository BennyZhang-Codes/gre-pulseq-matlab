function [seq, str_res, str_mat, str_r] = prep_Definition(seq, Actual, PE3D)

    % prepare sequence export
    res = round(1e3*Actual.fovRO/Actual.nRO, 2);
    a = fix(res);
    b = (res - a)*100;
    if mod(b, 10) == 0
        b = b/10;
    end
    str_res_RO = [num2str(a),'p',num2str(b)];

    res = round(1e3*Actual.fovPE/Actual.nPE, 2);
    a = fix(res);
    b = (res - a)*100;
    if mod(b, 10) == 0
        b = b/10;
    end
    str_res_PE = [num2str(a),'p',num2str(b)];

    res = round(1e3*Actual.fov3D/Actual.n3D, 2);
    a = fix(res);
    b = (res - a)*100;
    if mod(b, 10) == 0
        b = b/10;
    end
    str_res_3D = [num2str(a),'p',num2str(b)];

    str_res    = [str_res_RO, 'x', str_res_PE, 'x', str_res_3D];
    str_mat    = [num2str(Actual.nRO), 'x', num2str(Actual.nPE), 'x', num2str(Actual.n3D)];

    str_r      = [num2str(Actual.AccelerationPE), 'x', num2str(Actual.Acceleration3D)];


    % readout oversampling 
    % seq.setDefinition('ReadoutOversamplingFactor', readoutOS             );
    
    % % sequence definitions: additional information required by GRAPPA
    % seq.setDefinition('kSpaceCenterLine'     , kSpaceCenterLine          ); % PE center line index
    % seq.setDefinition('PhaseResolution'      , (fovRead/nX)/(fovPhase/nY)); % phase resolution
    seq.setDefinition('AccelerationFactorPE' , Actual.AccelerationPE     );   
    seq.setDefinition('AccelerationFactor3D' , Actual.Acceleration3D     );    
    seq.setDefinition('kSpaceCenterLine'     , PE3D.IdxCenter_PE - 1     ); % PE center line index
    seq.setDefinition('kSpaceCenterPartition', PE3D.IdxCenter_3D - 1     ); % PE center line index

    % seq.setDefinition('FirstFourierLine'     , FirstFourierLine          );  
    % seq.setDefinition('FirstRefLine'         , FirstRefLine              ); 
    % seq.setDefinition('FirstFourierLine'     , FirstFourier3D            );  
    % seq.setDefinition('FirstRefLine'         , FirstRef3D                ); 
    % 
    % seq.setDefinition('nRefLine'             , nRef                      ); % number of ACS line
    
    seq.setDefinition('FOV'                  , Actual.FOV                  );
    % seq.setDefinition('FOV'                  , [Actual.fovRO Actual.fovPE Actual.fov3D] );
    seq.setDefinition('MatrixSize'           , [Actual.nRO Actual.nPE Actual.n3D]       );
    seq.setDefinition('TR'                   , Actual.TR                        );
    seq.setDefinition('TE'                   , Actual.TE                 );
    seq.setDefinition('Excit_FlipAngle'      , Actual.flipEx             );
    seq.setDefinition('Excit_Duration'       , Actual.tEx_Slab           );
    seq.setDefinition('Excit_TBP'            , Actual.tbpEx_Slab         );

    seq.setDefinition('nDummy'               , Actual.nDummy             );
    seq.setDefinition('BW'                   , Actual.BWPerPixel         );
    seq.setDefinition('nRep'                 , Actual.nRep               );
    seq.setDefinition('nEcho'                , Actual.nTE                );

    
    if Actual.bBipolarROGrads
        ReadoutMode = 'Bipolar';
    else
        ReadoutMode = 'Monopolar';
    end
    seq.setDefinition('ReadoutMode'          , ReadoutMode               );

    seq.setDefinition('DimFast'              , Actual.DimFast            );

    seq.setDefinition('Developer'            , 'Jinyuan Zhang'           );
    seq.setDefinition('Name'                 , 'gre'                     );
end
