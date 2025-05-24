function [g, t, T_total] = design_gradient_min_time(A, tMax, G0, G1, Gmax, Smax, GradRasterTime)
    tMin = 0;
    iter_count = 0;
    max_iter = 1000;
    T_total = (tMin + tMax)/2;
    
    Continue = true;
    while Continue && iter_count <= max_iter
        errorOccurred = false;
        try
            [g, t] = design_gradient_waveform(A, T_total, G0, G1, Gmax, Smax, GradRasterTime);
        catch exception
            % fprintf('T_total %f is not feasible: ', T_total);
            % disp(exception.message);
            errorOccurred = true; 
        end
        if errorOccurred
            tMin = T_total;
            T_total = round((tMax + tMin)/2 / GradRasterTime) * GradRasterTime;
        else
            tMax = T_total;
            if (tMax - tMin) < GradRasterTime
                Continue = false;
            else
                T_total = round((tMax + tMin)/2 / GradRasterTime) * GradRasterTime;
            end
        end
        iter_count = iter_count + 1;
    end

    % fprintf('# design_gradient_min_time >>> : %.0f us | %.2f 1/m\n', t(end)*1e6, A);
end