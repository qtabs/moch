function pars = loadParameters()

    % use the estimulus function to define the stimulus parameters
    pars.est          = defineStimulus();

    % RUNNING PARAMETERS
    pars.onlySubcort  = false;      % used to run only the subcortical stage
    pars.verb         = false;      % controls the verbosity of the libs
    pars.freqInterval = [33, 2000]; % [minFreq, maxFreq] (Hz)
    pars.N            = 250;        % number of exc/inh populations 
    pars.cortFs       = 1000;       % sample rate in cortex (samples per sec)

    % NORMALISATION
    pars.subCortAff = 4;     % number of afferents in the subcortical pathway
    pars.regularise = 1;     % whether to regularise (1) or not (0) the SACF
    pars.mu0        = 75;    % expected peak activity after normalisation (Hz)
    pars.tauSubthal = 20;    % subthalamic lowpass tau (ms)
    pars.solvOnset  = 1;     % Solvability onset for the delay channels (ms)
    pars.SACFGround = 0.35;  % SACF baseline ([0,1]; -1 for average signal)
    pars.subDelay   = 50;    % Delay introduced by subcortical processing (ms)
    pars.subDelayDy = 75;    % subDelay for dyads (ms)

    % CONNECTIVITIES (nA)
    % Two cortical regions: decoder (P) and sustainer (Q)

    % conductivities from subcortical aras to P/E
    pars.Jtn  = 0.00;      % NMDA to P-excitatory
    pars.Jta  = 2.7;       % AMPA to P-excitatory 
    
    % conductivities P to P
    pars.Jee  = 0.14;      % NMDA to excitatory
    pars.Jae  = 9.9024E-4; % AMPA to excitatory [Wong2006]
    pars.Jie  = 0.53;      % GABA to excitatory
    pars.Jei  = 0.17;      % NMDA to inhibitory
    pars.Jai  = 6.517E-5;  % AMPA to inhibitory [Wong2006]
    pars.Jii  = 0.11;      % GABA to inhibitory

    % conductivities Q to Q
    pars.Jqee = 0.25;      % NMDA to inhibitory
    pars.Jqae = 9.9024E-4; % AMPA to inhibitory [Wong2006]
    pars.Jqie = 0.80;      % GABA to excitatory 
    pars.Jqei = 0.00;      % NMDA to inhibitory
    pars.Jqai = 6.517E-5;  % GABA to inhibitory [Wong2006]
    pars.Jqii = 0.00;      % NMDA to excitatory
    pars.QI0i = 0.26;      % Constant input at Q populations
    pars.QI0e = 0.18;      % Constant input at Q populations

    % connectivities from P to Q (bottom-up)
    pars.Jpqen = 0.00;     % NMDA to excitatory
    pars.Jpqea = 0.35;     % AMPA to excitatory
    pars.Jpqeg = 0.00;     % GABA to excitatory
    pars.Jpqin = 0.00;     % NMDA to inhibitory
    pars.Jpqia = 0.00;     % AMPA to inhibitory 
    pars.Jpqig = 0.45;     % GABA to inhibitory

    % connectivities from Q to P (top-down)
    pars.Jqpen = 0;        % NMDA to excitatory
    pars.Jqpea = 0;        % AMPA to excitatory
    pars.Jqpeg = 0;        % GABA to excitatory
    pars.Jqpin = 0.45;     % NMDA to inhibitory
    pars.Jqpia = 0;        % AMPA to inhibitory
    pars.Jqpig = 0;        % GABA to inhibitory

    % structural connectivity parameters
    pars.cc  = 0.10;   % Cie par: baseline uniform inhibition
    pars.cS  = 0.10;   % Cie par: self uniform inhibition
    pars.kie = 66;     % Cie par: number of target ensembles
    pars.kei = 3;      % Cei par: number of target ensembles
    pars     = connectivities(pars); % Matrix generator (some hardcoded pars)


    % ENSEMBLE AND SYNAPTIC DYNAMICS

    % adaptation and delays
    pars.ensembleDelays = true; % delay in the connect. between ensembles?
    pars.selfDelay      = 0;    % delay between closest populations (ms)
    pars.maxDelay       = 0;    % delay between farest populations (ms)
    pars.tauAdapt       = 100;  % adaptation time constant (ms)
    pars.adaptStr       = 1e-5;    % adaptation strength (0 for no adaptation)

    % time constants
    pars.tauNMDA    = 30;    % NMDA-gating time constants (ms)
    pars.tauGABA    = 5;     % GABA-gating time constants (ms) [Brunel 2001]
    pars.tauAMPA    = 2;     % AMPA-gating time constants (ms) [Brunel 2001]
    pars.tauEff     = true;  % activates dynamic taus          [Ostojic 2011]
    pars.tauPop     = 10;    % population time constants  (ms) [Ostojic 2011]
    pars.tauSACF    = -1;    % sacf tau (-1 for Wiegribe's taus; ms)

    % population parameters
    pars.gamma = 0.641/1000; % NMDA coupling   [Brunel 2001]
    pars.sigma = 0.0007;   % noise amplitude [Wong 2006]
    pars.I0e   = 0.3146;     % constant population input, exc (nA)
    pars.I0i   = 0.1501;     % constant population input, inh (nA)
    pars.ae    = 310;        % excit. non-lin. pars (1/VnC) [Wong 2006]
    pars.be    = 125;        % excit. non-lin. pars (Hz) [Wong 2006]
    pars.de    = 0.16;       % excit. non-lin. pars (s) [Wong 2006]
    pars.ai    = 615;        % inhib. non-linearity pars  (1/VnC)[Wong 2006]
    pars.bi    = 177;        % inhib. non-linearity pars (Hz) [Wong 2006]
    pars.di    = 0.087;      % inhib. non-linearity pars [Wong 2006]

end


function est = defineStimulus()
  
    % estimulus parameters 
    est.filename = -1;    % path to wav file (-1 to generate the stimulus)
    est.loud     = 70;    % loudness rescaling (-1 to use original; dB SPL)
    est.interval = -1;    % wavefile interval (False to use whole file; ms)
    est.onset    = 0;     % adds a silence before stimulus onset (ms)
    est.tail     = 0;     % adds a stimulus silent tail (ms)
    est.bandpass = false; % bandpass filtering (false for none; [,] Hz)

    % if filename == -1, the lib creates a estimulus
    % Available stimuli: PT (pure tone)
    %                    HCT (harmonic complex tone in sin phase)
    %                    altHCT (harmonic complex tone in alternate phase)
    %                    IRN (iterated rippled noise)
    %                    noiseIRN (IRN preceeded by gaussian noise)
    %                    noise (gaussian noise)
    %                    CT (click train; not tested)
    %                    chord/HCTchord/IRNchord (chords, use "notes")
    %                    FMsweep (est.f = avg. freq; est.shift = freq gap)
    est.type      = 'IRN';
    est.f         = 200;     % fundamental/equivalent frequency (Hz)
    est.dur       = 400;     % duration (ms)
    est.harms     = 15:20;   % list of the harmonics (only for HTC*) 
    est.harmFact  = 1;       % attenuation multip. factor (only for HTC*)  
    est.shift     = 0;       % frequency shift in HCTs / FM sweeps (Hz)
    est.maskNoise = 0;       % ratio of masking noise
    est.nOfIts    = 32;      % number of interations (only for IRN*)
    est.notes     = [0 1];   % list of notes (only for *chords)
    est.tuning    = 'just'; % chord intonation ('equal', 'just')
    est.noiseOff  = 0;     % duration of noise offset (only for IRNnose)
    est.save      = 1;       % 1 for saving the stimulus .wav, 0 otherwise

end


function pars = connectivities(pars)

    % 0. Generating lagSpace
    lagMin = 1000 / pars.freqInterval(2);
    lagMax = 1000 / pars.freqInterval(1);
    lagSpace = linspace(lagMax, lagMin, pars.N)';


    % 1. Cei
    Cei = zeros(pars.N);
    for i = 1:pars.N
        for r = 1:pars.kei
            [dif, ind] = min((lagSpace/lagSpace(i) - r).^2);
            if (lagSpace(max(ind - 1, 1)) > r * lagSpace(i)) %&& dif < 1
                Cei(i, ind) = 1;
            end
        end
        Cei(i, i) = 1;
    end


    % 2. Cii
    for i = 1:pars.N
        for j = 1:pars.N
            if abs(i - j) == 0
                Cii(i, j) = 0;
            elseif i > j
                Cii(i, j) = 1;
            else
                Cii(i, j) = 1;
            end
        end
    end

    % 3. Cie
    Cie = zeros(pars.N);
    delta = 2 * (lagSpace(1) - lagSpace(2));
    for j = 1:pars.N
        for i = 1:pars.N
            for r = 2:pars.kie    
                diss = abs(lagSpace(i) - r*lagSpace(j)) / (delta * min(r, 5));
                if diss < 1
                    Cie(i, j) = sqrt(1 - diss);
                end
            end
        end
    end
    Cie(Cie < pars.cc) = pars.cc;
    Cie = Cie + (pars.cS - pars.cc) * eye(size(Cie));

    % 4. Cee
    Cee = eye(pars.N);

    % 5. Connectivities for Q
    Cqei = eye(pars.N);
    Cqii = eye(pars.N);
    Cqie = eye(pars.N);
    Cqee = eye(pars.N);
    Cpq  = eye(pars.N);


    % Visualisation
    if false
        subplot(221); 
        imagesc(lagSpace,lagSpace,Cei);title('Cei');colorbar();caxis([0 1]);
        subplot(222); 
        imagesc(lagSpace,lagSpace,Cii);title('Cii');colorbar();caxis([0 1]);
        subplot(223); 
        imagesc(lagSpace,lagSpace,Cie);title('Cie');colorbar();caxis([0 1]);
        subplot(224); 
        imagesc(lagSpace,lagSpace,Cee);title('Cee');colorbar();caxis([0 1]);
    end

    pars.Cee = Cee;
    pars.Cei = Cei;
    pars.Cie = Cie;
    pars.Cii = Cii;
    pars.Cqee = Cqee;
    pars.Cqei = Cqei;
    pars.Cqie = Cqie;
    pars.Cqii = Cqii;

end