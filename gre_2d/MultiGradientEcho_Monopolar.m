
Nx = 400;
roDuration = 1.0e-3;
% set system limits
sys = mr.opts('MaxGrad', 50, 'GradUnit', 'mT/m', ...
    'MaxSlew', 140, 'SlewUnit', 'T/m/s', ... 
    'rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6);

% basic parameters
seq=mr.Sequence(sys);           % Create a new sequence object
fov=200e-3; 
% Nx=200; 
Ny=Nx;      % Define FOV and resolution
alpha=10;                       % flip angle
sliceThickness=1e-3;            % slice
TR=25e-3;                       % TR, a single value
prepscans = 10; % number of dummy preparation scans
nTE = 6;
TE1 = 4.08*1e-3;  % first echo time 
esp = 2.04*1e-3;
TE=(0:nTE-1) * esp + TE1;               % give a vector here to have multiple TEs % TODO: reduce TEs
%TODO: play with TEs to make them really minimal

% more in-depth parameters
rfSpoilingInc=117;              % RF spoiling increment
rfDuration=2e-3;
% roDuration=1.8e-3;              % not all values are possible, watch out for the checkTiming output
adcDwell = round(roDuration/Nx/seq.adcRasterTime)* seq.adcRasterTime;
roDuration = adcDwell * Nx;
BWPerPixel = 1/roDuration;
disp(BWPerPixel)
%%

% Create alpha-degree slice selection pulse and corresponding gradients 
[rf, gz, gzReph] = mr.makeSincPulse(alpha*pi/180,'Duration',rfDuration,...
    'SliceThickness',sliceThickness,'apodization',0.42,'timeBwProduct',4,'system',sys);

% Define other gradients and ADC events
deltak=1/fov; % Pulseq default units for k-space are inverse meters
% gx = mr.makeTrapezoid('x','FlatArea',Nx*deltak,'FlatTime',roDuration,'system',sys, 'riseTime', (esp-roDuration)/2, 'fallTime', (esp-roDuration)/2); % Pulseq default units for gradient amplitudes are 1/Hz
gx = mr.makeTrapezoid('x','FlatArea',Nx*deltak,'FlatTime',roDuration,'system',sys); % Pulseq default units for gradient amplitudes are 1/Hz

gxFlyBack = mr.makeTrapezoid('x','Area',-gx.area,'system',sys);
adc = mr.makeAdc(Nx,'Duration',gx.flatTime,'Delay',gx.riseTime,'system',sys);
gxPre = mr.makeTrapezoid('x','Area',-gx.area/2,'system',sys); % if no 'Duration' is provided shortest possible duration will be used
phaseAreas = ((0:Ny-1)-Ny/2)*deltak;

% gradient spoiling
% gxSpoil=mr.makeTrapezoid('x','Area',2*Nx*deltak*spSign,'system',sys);      % 2 cycles over the voxel size in X
gxSpoil=mr.makeExtendedTrapezoidArea('x',gx.amplitude,0,2*Nx*deltak*spSign,sys); 
gzSpoil=mr.makeTrapezoid('z','Area',4/sliceThickness,'Delay',gx.delay+gx.riseTime+gx.flatTime,'system',sys); % 4 cycles over the slice thickness

% Calculate timing (need to decide on the block structure already)
helperT=ceil((gz.fallTime + gz.flatTime/2 + gx.riseTime + gx.flatTime/2)/seq.gradRasterTime)*seq.gradRasterTime;

delayTE = zeros(size(TE)) ;

for c=1:length(TE)
    delayTE(c)=TE(c) - helperT;
    helperT=helperT+delayTE(c)+mr.calcDuration(gx);
end

%%
assert(all(delayTE(1)>=mr.calcDuration(gxPre,gzReph)));
assert(all(delayTE(2:end)>=0));
delayTR=round((TR - mr.calcDuration(gz) - sum(delayTE) ...
    - mr.calcDuration(gx)*(length(TE)-1))/seq.gradRasterTime)*seq.gradRasterTime;
assert(all(delayTR>=mr.calcDuration(gxSpoil,gzSpoil)));

% initialize the RF spoling counters 
rf_phase=0;
rf_inc=0;


% define sequence blocks
for i=-(prepscans-1):Ny % loop over phase encodes
    rf.phaseOffset=rf_phase/180*pi;
    adc.phaseOffset=rf_phase/180*pi;
    rf_inc=mod(rf_inc+rfSpoilingInc, 360.0);
    rf_phase=mod(rf_phase+rf_inc, 360.0);
    %
    seq.addBlock(rf,gz);
    if (i>0)
        gyPre = mr.makeTrapezoid('y','Area',phaseAreas(i),'Duration',mr.calcDuration(gxPre),'system',sys);
    else
        gyPre = mr.makeTrapezoid('y','Area',0,'Duration',mr.calcDuration(gxPre),'system',sys);
    end
    for c=1:nTE % loop over TEs
        if (c==1)
            seq.addBlock(mr.align('left', mr.makeDelay(delayTE(c)),gyPre,gzReph,'right',gxPre));
            seq.addBlock(mr.makeLabel('SET','ECO',0));
        else
            seq.addBlock(mr.align('left', mr.makeDelay(delayTE(c)),'right',gxFlyBack));
        end
        
        if (c==nTE)
            gx1=mr.makeExtendedTrapezoid(gx.channel,...
                'times', [gx.delay gx.riseTime+gx.delay gxSpoil.tt+gx.delay+gx.riseTime+gx.flatTime],...
                'amplitudes',[0 gx.amplitude  gxSpoil.waveform],'system',sys);
            gyRep = mr.scaleGrad(gyPre,-1);
            gyRep.delay = gx.delay+gx.riseTime+gx.flatTime;
            if (i>0), seq.addBlock(mr.makeDelay(delayTR),gx1,adc, gzSpoil, gyRep); else, seq.addBlock(mr.makeDelay(delayTR),gx,gzSpoil, gyRep); end
        else
            if (i>0), seq.addBlock(gx,adc); else, seq.addBlock(gx); end
        end
        seq.addBlock(mr.makeLabel('INC','ECO',1));
        % to check/debug TE calculation with seq.testReport() comment out
        % the above line and uncommend the line below this comment block; 
        % change 'c==3' statement to select the echo to test
        %if c==2, seq.addBlock(gx,adc); else, seq.addBlock(gx); end
    end
    if (i>0), seq.addBlock(mr.makeLabel('INC','LIN',1)); end
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
seq.setDefinition('ReadoutMode', 'Monopolar');
seq.setDefinition('Dummy', prepscans);
seq.setDefinition('ESP', esp);
seq.setDefinition('BW', BWPerPixel);
seq.setDefinition('ADC_DwellTime', adcDwell);
seq.setDefinition('Name', sprintf('m%s_%sp%s', num2str(nTE), num2str(a), num2str(b)));
seq.setDefinition('Developer', 'Jinyuan Zhang');

seq.write(sprintf('grem%se_%s_tr%s_fa%s_bw%s.seq', num2str(nTE), prefix, num2str(TR*1e3), num2str(alpha), num2str(round(BWPerPixel))));

%% plot sequence and k-space diagrams
seq.plot('timeRange', [6 15]*TR, 'TimeDisp', 'ms', 'Label', 'LIN,ECO');
% seq.plot('TimeDisp', 'ms', 'Label', 'LIN,ECO');
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
% rep = seq.testReport;
% fprintf([rep{:}]);
