function [seq] = check_PNS(seq, Actual)

    switch lower(Actual.ScannerType)
        case 'terra-xr'
            asc_file = 'MP_GPA_K2298_2250V_793A_SC72CD_EGA.asc';
        case 'terra-xj'
            asc_file = 'MP_GPA_K2259_2000V_650A_SC72CD_EGA.asc';
        otherwise
            error('Unsupported ScannerType: %s', Actual.ScannerType);
    end
    
    [pns_ok, pns_n, pns_c, tpns] = seq.calcPNS(asc_file);

    if (pns_ok)
        fprintf('PNS check passed successfully\n');
    else
        fprintf('PNS check failed! The sequence will probably be stopped by the Gradient Watchdog\n');
    end
end
