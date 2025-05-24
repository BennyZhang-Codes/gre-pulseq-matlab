function [g, t] = design_gradient_waveform(A, T_total, G0, G1, Gmax, Smax, GradRasterTime)
% DESIGN_GRADIENT_WAVEFORM Generates trapezoidal gradient waveform with sign-sensitive handling
% Inputs:
%   A              : Target gradient area (mT/m·s), sign determines polarity
%   T_total        : Total waveform duration (s)
%   G0             : Initial gradient (mT/m), can be positive/negative
%   G1             : Final gradient (mT/m), can be positive/negative
%   Gmax           : Maximum gradient amplitude (absolute value, mT/m)
%   Smax           : Maximum slew rate (absolute value, T/m/s)
%   GradRasterTime : System gradient raster time (s)

    t_start = 0; % Waveform start time (s)

    % Validate gradient polarity constraints
    if abs(G0) > Gmax || abs(G1) > Gmax
        error('Initial/Final gradient exceeds maximum limit');
    end
    
    % Align total duration to raster grid
    n_total = round(T_total / GradRasterTime);
    T_total_adj = n_total * GradRasterTime;
    
    [A_min, A_max] = calculate_min_max_area(T_total, G0, G1, Gmax, Smax, GradRasterTime);
    if A < A_min || A > A_max
        error('Requested area %.3f is outside the feasible range [%.3f, %.3f]', A, A_min, A_max);
    end
    
    % Solve trapezoid parameters with sign handling
    [Gp, t_ramp1, t_plat, t_ramp2, exitflag] = find_platform_gradient(...
        A, T_total_adj, G0, G1, Gmax, Smax, GradRasterTime);
    
    if exitflag > 0
        % Construct time points sequence
        t0 = t_start;
        t1 = t0 + t_ramp1;
        t2 = t1 + t_plat;
        t3 = t2 + t_ramp2;
        
        t = [t0, t1, t2, t3];
        g = [G0, Gp, Gp, G1];

        % fprintf('# design_gradient_waveform >>> : %.0f us | %.2f 1/m\n', t(end)*1e6, A);
        return;
    end
    
    error('No valid waveform found within constraints');
end


function [A_min, A_max] = calculate_min_max_area(T_total, G0, G1, Gmax, Smax, dt)
% Calculate the maximum possible positive and negative areas
    [~, A_max] = compute_extreme_area(Gmax, T_total, G0, G1, Smax, dt);
    [~, A_min] = compute_extreme_area(-Gmax, T_total, G0, G1, Smax, dt);
end

function [feasible, A] = compute_extreme_area(Gp, T_total, G0, G1, Smax, dt)
% Calculate the actual area for the given Gp
    t_ramp1 = round(abs(Gp - G0)/(Smax*dt))*dt;
    t_ramp2 = round(abs(G1 - Gp)/(Smax*dt))*dt;
    t_plat = T_total - t_ramp1 - t_ramp2;
    
    if t_plat < 0
        feasible = false;
        A = 0;
    else
        A = 0.5*(G0 + Gp)*t_ramp1 + Gp*t_plat + 0.5*(Gp + G1)*t_ramp2;
        feasible = true;
    end
end


function [Gp, t_ramp1, t_plat, t_ramp2, exitflag] = find_platform_gradient(...
    A, T_total, G0, G1, Gmax, Smax, dt)

    % Intelligently determine the search interval
    [initial_lb, initial_ub] = determine_search_bounds(A, G0, G1, Gmax);
    
    % Sign-sensitive numerical solver
    options = optimset('Display','off', 'TolX', 1e-6);
    try
        [Gp, ~, exitflag] = fminbnd(@(Gp) abs(area_eq(Gp, A, T_total, G0, G1, Smax, dt)-A),...
            initial_lb, initial_ub, options);
    catch ME
        if contains(ME.message, 'function values at interval endpoints')
            error('No solution found. Please check whether the input parameters are within the physically feasible range.');
        else
            rethrow(ME);
        end
    end
    
    if exitflag <= 0
        t_ramp1 = 0; t_plat = 0; t_ramp2 = 0; 
        return;
    end
    
    % Direction-aware ramp time calculation
    delta1 = Gp - G0;
    t_ramp1 = max(ceil(abs(delta1)/(Smax*dt)), 1) * dt;
    if sign(delta1) ~= sign(Gp - G0)
        error('Ramp1 direction inconsistency detected');
    end
    
    delta2 = G1 - Gp;
    t_ramp2 = max(ceil(abs(delta2)/(Smax*dt)), 1) * dt;
    if sign(delta2) ~= sign(G1 - Gp)
        error('Ramp2 direction inconsistency detected');
    end
    
    t_plat = T_total - t_ramp1 - t_ramp2;
    t_plat = max(round(t_plat/dt), 0) * dt;
    
    % Verify temporal constraints
    T_total_new = t_ramp1 + t_plat + t_ramp2;
    if abs(T_total_new - T_total) > dt/2
        exitflag = -1;
        return;
    end
    
    % Recalculate with sign preservation
    a = 0.5*(t_ramp1 + t_ramp2) + t_plat;
    b = 0.5*(G0*t_ramp1 + G1*t_ramp2);
    Gp_new = (A - b)/a;

    Gp = Gp_new;
    exitflag = 1;
end

function F = area_eq(Gp, A, T_total, G0, G1, Smax, dt)
% Sign-sensitive area calculation

    % Direction-consistent ramp durations
    t_ramp1 = max(ceil(abs(Gp - G0)/(Smax*dt)), 1)*dt;
    t_ramp2 = max(ceil(abs(G1 - Gp)/(Smax*dt)), 1)*dt;
    t_plat = T_total - t_ramp1 - t_ramp2;
    
    if t_plat < 0
        F = inf;
    else
        % Signed area components
        A_ramp1 = 0.5*(G0 + Gp)*t_ramp1;
        A_plat = Gp*t_plat;
        A_ramp2 = 0.5*(Gp + G1)*t_ramp2;
        
        F = (A_ramp1 + A_plat + A_ramp2) - A;
    end
end

function [lb, ub] = determine_search_bounds(A, G0, G1, Gmax)

    if A > 0
        % lb = max(-Gmax, min([G0, G1]) - 0.2*Gmax);
        lb = -Gmax;
        ub = Gmax;
    else
        lb = -Gmax;
        ub = Gmax;
        % ub = min(Gmax, max([G0, G1]) + 0.2*Gmax);
    end
    
    % make sure the bounds
    if lb >= ub
        lb = -Gmax;
        ub = Gmax;
    end
end