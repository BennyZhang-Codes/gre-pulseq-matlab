function [seq] = check_PNS(seq)
    [pns_ok, pns_n, pns_c, tpns] = seq.calcPNS('MP_GPA_K2259_2000V_650A_SC72CD_EGA.asc'); % TERRA-XJ

    if (pns_ok)
        fprintf('PNS check passed successfully\n');
    else
        fprintf('PNS check failed! The sequence will probably be stopped by the Gradient Watchdog\n');
    end
end
