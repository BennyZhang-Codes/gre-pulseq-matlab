function [Actual] = prep_TRTE(Actual, params, sys)

    % For Siemens Recon - ensure nTE is consistent with UI
    Actual.nPE = round((round(params.nPE/params.nRO*100)/100) * params.nRO) ;
    if (Actual.nPE ~= params.nPE)
        warning('nPE updated to be consistent with Siemens UI (integer percent of nRO): %d --> %d', params.nPE, Actual.nPE) ;
    end
    
    % Ensure TE & TR are on gradient raster time
    % Round TEs (supports multi TEs for multiple echos)
    bTEPos = (Actual.TE > 0) ; % non-positive TEs mean use shortest TE
    Actual.TE(bTEPos) = round(params.TE(bTEPos) / sys.gradRasterTime) * sys.gradRasterTime ;
    for iTE = 1:Actual.nTE % Warn user if this has an effect.
        if (bTEPos(iTE) && abs(Actual.TE(iTE) - params.TE(iTE)) > eps(sys.gradRasterTime))
            warning('TE(%d) updated to be on the gradient raster: %f ms --> %f ms', iTE, params.TE(iTE)*1e3, Actual.TE(iTE)*1e3) ;
        end
    end
    
    % Round user defined TR
    Actual.TR = round(params.TR / sys.gradRasterTime) * sys.gradRasterTime ;
    if (abs(Actual.TR - params.TR) > eps(sys.gradRasterTime)) % Warn user if this has an effect.
        warning('TR updated to be on the gradient raster: %f ms --> %f ms', params.TR*1e3, Actual.TR*1e3) ;
    end

end
