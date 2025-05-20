
% for Nx = [200] % 200, 280, 332, 400
%     for roDuration = [1.0e-3, 1.2e-3, 1.5e-3, 1.6e-3, 1.8e-3]
Nx = 400;
roDuration = 1.2e-3;
% set system limits
sys = mr.opts('MaxGrad', 40, 'GradUnit', 'mT/m', ...
    'MaxSlew', 180, 'SlewUnit', 'T/m/s', ... 
    'rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6);

% sys = mr.opts('MaxGrad', 50, 'GradUnit', 'mT/m', ...
%     'MaxSlew', 140, 'SlewUnit', 'T/m/s', ... 
%     'rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6);

% basic parameters
seq=mr.Sequence(sys);           % Create a new sequence object
fov=200e-3; 
% Nx=400; c
Ny=Nx;      % Define FOV and resolution
alpha=10;                       % flip angle
sliceThickness=2e-3;            % slice
TR=25e-3;                       % TR, a single value
prepscans = 40; % number of dummy preparation scans
nTE = 6;
TE1 = 3.06*1e-3;  % first echo time 
esp = 1.02*1e-3;
TE=(0:nTE-1) * esp + TE1;
%TODO: play with TEs to make them really minimal

% more in-depth parameters
rfSpoilingInc=117;              % RF spoiling increment
rfDuration=2.0e-3;
% roDuration=1.0e-3;              % not all values are possible, watch out for the checkTiming output
% roDuration=1.2e-3;
adcDwell = round(roDuration/Nx/seq.adcRasterTime) * seq.adcRasterTime;
roDuration = adcDwell * Nx;
BWPerPixel = 1/roDuration;
% disp(roDuration)
% disp(BWPerPixel)
%%

% Create alpha-degree slice selection pulse and corresponding gradients 
[rf, gz, gzReph] = mr.makeSincPulse(alpha*pi/180,'Duration',rfDuration,...
    'SliceThickness',sliceThickness,'apodization',0.42,'timeBwProduct',4,'system',sys);
% disp(gz.fallTime)
% disp(gz.flatTime)
%%
% Define other gradients and ADC events
deltak=1/fov; % Pulseq default units for k-space are inverse meters
gxp = mr.makeTrapezoid('x','FlatArea',Nx*deltak,'FlatTime',roDuration,'system',sys); % Pulseq default units for gradient amplitudes are 1/Hz
adc = mr.makeAdc(Nx,'Duration',gxp.flatTime,'Delay',gxp.riseTime,'system',sys);
gxPre = mr.makeTrapezoid('x','Area',-gxp.area/2,'system',sys); % if no 'Duration' is provided shortest possible duration will be used
phaseAreas = ((0:Ny-1)-Ny/2)*deltak;

% gradient spoiling
spSign=1;
% gxSpoil=mr.makeTrapezoid('x','Area',2*Nx*deltak*spSign,'system',sys);      % 2 cycles over the voxel size in X
gxSpoil=mr.makeExtendedTrapezoidArea('x',gxp.amplitude*spSign,0,2*Nx*deltak*spSign,sys); 
gzSpoil=mr.makeTrapezoid('z','Area',4/sliceThickness,'Delay',gxp.delay+gxp.riseTime+gxp.flatTime,'system',sys); % 4 cycles over the slice thickness

gx=mr.makeExtendedTrapezoid(gxp.channel,...
    'times', [gxp.delay gxp.riseTime+gxp.delay gxSpoil.tt+gxp.delay+gxp.riseTime+gxp.flatTime],...
    'amplitudes',[0 gxp.amplitude  gxSpoil.waveform],'system',sys);
% Calculate timing (need to decide on the block structure already)
delayTE = zeros(size(TE)) ;
for c=1:length(TE)
    delayTE(c)=TE(c) - ceil((gz.fallTime + gz.flatTime/2 + gxp.riseTime + gxp.flatTime/2)/seq.gradRasterTime)*seq.gradRasterTime;
end
assert(all(delayTE>=mr.calcDuration(gxPre,gzReph)));
%%
delayTR=round((TR - mr.calcDuration(gz) - delayTE)/seq.gradRasterTime)*seq.gradRasterTime;

assert(all(delayTR>=mr.calcDuration(gxSpoil,gzSpoil)));

% initialize the RF spoling counters 
rf_phase=0;
rf_inc=0;
for c=1:length(TE) % loop over TEs
    if (c==1)
        seq.addBlock(mr.makeLabel('SET','ECO',0));
        % define dummy preparation blocks
        for i=1:prepscans
            rf.phaseOffset=rf_phase/180*pi;
            adc.phaseOffset=rf_phase/180*pi;
            rf_inc=mod(rf_inc+rfSpoilingInc, 360.0);
            rf_phase=mod(rf_phase+rf_inc, 360.0);
            %
            seq.addBlock(rf,gz);
            gyPre = mr.makeTrapezoid('y','Area',0,'Duration',mr.calcDuration(gxPre),'system',sys);
            seq.addBlock(mr.align('left', mr.makeDelay(delayTE(c)),gyPre,gzReph,'right',gxPre)); 
            gyRep = mr.scaleGrad(gyPre,-1);
            gyRep.delay = gxp.delay+gxp.riseTime+gxp.flatTime;
            seq.addBlock(mr.makeDelay(delayTR(c)), gx, gzSpoil, gyRep);
        end
    end
    seq.addBlock(mr.makeLabel('SET','LIN',0));

    % define sequence blocks
    for i=1:Ny % loop over phase encodes
        rf.phaseOffset=rf_phase/180*pi;
        adc.phaseOffset=rf_phase/180*pi;
        rf_inc=mod(rf_inc+rfSpoilingInc, 360.0);
        rf_phase=mod(rf_phase+rf_inc, 360.0);
        %
        seq.addBlock(rf,gz);
        gyPre = mr.makeTrapezoid('y','Area',phaseAreas(i),'Duration',mr.calcDuration(gxPre),'system',sys);
        seq.addBlock(mr.align('left', mr.makeDelay(delayTE(c)),gyPre,gzReph,'right',gxPre)); 
        
        gyRep = mr.scaleGrad(gyPre,-1);
        gyRep.delay = gxp.delay+gxp.riseTime+gxp.flatTime;
        seq.addBlock(mr.makeDelay(delayTR(c)), gx,adc, gzSpoil, gyRep);
        seq.addBlock(mr.makeLabel('INC','LIN',1));
    end
    seq.addBlock(mr.makeLabel('INC','ECO',1));
end
%% check whether the timing of the sequence is correct
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

%% prepare sequence export
res = round(1e3*fov/Nx, 2);
a = fix(res);
b = (res - a)*100;
if mod(b, 10) == 0
    b = b/10;
end
prefix = [num2str(a),'p',num2str(b),'_',num2str(Nx),'_r1'];

seq.setDefinition('FOV', [fov fov sliceThickness]);
seq.setDefinition('MatrixSize', [Nx Ny 1]);
seq.setDefinition('FlipAngle', alpha);
seq.setDefinition('SliceThickness', sliceThickness);
seq.setDefinition('TR', TR);
seq.setDefinition('TE', TE);
seq.setDefinition('ReadoutMode', 'Separate');
seq.setDefinition('Dummy', prepscans);
seq.setDefinition('ESP', esp);
seq.setDefinition('BW', BWPerPixel);
seq.setDefinition('ADC_DwellTime', adcDwell);
seq.setDefinition('Name', sprintf('s%s_%sp%s', num2str(nTE), num2str(a), num2str(b)));
seq.setDefinition('Developer', 'Jinyuan Zhang');

seq.write(sprintf('gres%se_%s_tr%s_fa%s_bw%s.seq', num2str(nTE), prefix, num2str(TR*1e3), num2str(alpha), num2str(round(BWPerPixel))));

%% plot sequence and k-space diagrams
seq.plot('TimeDisp', 'ms', 'Label', 'LIN,ECO');
% seq.plot('timeRange', [0 10]*TR, 'TimeDisp', 'ms', 'Label', 'LIN,ECO');
% seq.plot('timeRange', [0 30]*TR, 'TimeDisp', 'ms', 'Label', 'LIN,ECO');
% ('timeRange', [0 nTRs]*this.seq_params.TR, 'TimeDisp', 'ms', 'label', 'lin')
% k-space trajectory calculation
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();

% plot k-spaces
figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D k-space plot
axis('equal'); % enforce aspect ratio for the correct trajectory display
hold;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % plot the sampling points
title('full k-space trajectory (k_x x k_y)');

%% PNS calc
warning('OFF', 'mr:restoreShape');
[pns_ok, pns_n, pns_c, tpns] = seq.calcPNS('MP_GPA_K2259_2000V_650A_SC72CD_EGA.asc'); % TERRA-XJ

if (pns_ok)
    fprintf('PNS check passed successfully\n');
else
    fprintf('PNS check failed! The sequence will probably be stopped by the Gradient Watchdog\n');
end

%% very optional slow step, but useful for testing during development e.g. for the real TE, TR or for staying within slewrate limits  
rep = seq.testReport;
fprintf([rep{:}]);
%     end
% end

