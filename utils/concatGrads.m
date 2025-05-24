function [g] = concatGrads(g_list, sys)
    g_t = [];
    g_a = [];
    t_offset = 0;
    g_1 = g_list{1};
    
    for i = 1:length(g_list)
        g_i = g_list{i};

        switch lower(g_i.type)
            case 'grad'
                tt  = g_i.tt';
                a   = g_i.waveform';
                dur = g_i.shape_dur; 
                delay = g_i.delay;
            case 'trap'
                if g_i.flatTime > 0
                    tt  = [0 g_i.riseTime  g_i.riseTime+g_i.flatTime g_i.riseTime+g_i.flatTime+g_i.fallTime];
                    a   = [0 g_i.amplitude g_i.amplitude             0                                     ];
                else
                    tt  = [0 g_i.riseTime  g_i.riseTime+g_i.fallTime];
                    a   = [0 g_i.amplitude 0                        ];
                end
                dur = g_i.riseTime+g_i.flatTime+g_i.fallTime;
                delay = g_i.delay;
            otherwise
                error('invalid grad type');
        end

        if i == 1
            g_t = [g_t tt];
            g_a = [g_a a];
        else
            if delay == 0
                g_t = [g_t t_offset + tt(2:end)];
                g_a = [g_a a(2:end)];
            else
                t_offset = t_offset + delay;
                g_t = [g_t t_offset  t_offset + tt(2:end)];
                g_a = [g_a 0         a(2:end)];
            end

        end
        t_offset = t_offset + dur;
    end

    

    g = mr.makeExtendedTrapezoid(g_1.channel, 'times', g_t, 'amplitudes', g_a, 'system', sys);
    g.delay = g_1.delay;
end