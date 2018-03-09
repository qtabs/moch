N        = 60;
notes    = [0:1:12];
dyads    = {'unison', 'minor second', 'second', 'minor third', ...
            'third', 'fourth', 'tritone', 'perfect fifth', ...
            'minor sixth', 'sixth', 'minor seventh', 'seventh', 'octave'};

parbase = loadParameters();
parbase.est.dur      = 200; 
parbase.est.type     = 'IRNchordSS';
parbase.est.f        = 160;
parbase.est.nOfIts   = 8;
parbase.est.bandpass = [125, 2000];
parbase.est.tuning   = 'just';

onset(2:length(notes)) = parbase.subDelayDy;
onset(1)               = parbase.subDelay;
parbase.subDelay   = 0;
parbase.subDelayDy = 0;

for i = 1:length(notes)
    pars{i} = parbase;
    pars{i}.est.notes = [0, notes(i)]; 
end

label = {'100Hz', '200Hz', '32its', 'equal'};
for i = 1:length(label)
    parMat{i} = pars;
end

for j = 1:length(notes)
    parMat{1}{j}.est.f = 100;
    parMat{2}{j}.est.f = 200;
    parMat{3}{j}.est.tuning = 'equal'; 
    parMat{4}{j}.est.nOfIts = 32;
end

[~, r] = tdoch(pars{1});
lagSpace  = r.lagSpace;
timeSpace = r.timeSpace;

for ix = 1:length(label)

    pars = parMat{ix};

    for i = 1:length(notes)
        tt = tic;
        fprintf(' - %d of %d ...\n', i, length(notes));
        parfor n = 1:N
            [s, r]  = tdoch(pars{i});
            ACPar{n}(i, :)    = mean(r.A(176:200, :), 1);
            DePar{n}(i, :)    = mean(s.p.He(176:200, :), 1);
            DiPar{n}(i, :)    = mean(s.p.Hi(176:200, :), 1);
            SePar{n}(i, :)    = mean(s.q.He(176:200, :), 1);
            SiPar{n}(i, :)    = mean(s.q.Hi(176:200, :), 1);
            [~, latPar{n}(i)] = max(mean(s.p.He, 2));
        end
        fprintf('done! time: %.0fm\n', toc(tt)/60);
    end 

    for n = 1:N
        ACMat(:, :, n) = ACPar{n};
        DeMat(:, :, n) = DePar{n};
        DiMat(:, :, n) = DiPar{n};
        SeMat(:, :, n) = SePar{n};
        SiMat(:, :, n) = SiPar{n};
        lat0(:, n)     = latPar{n};
    end

    AC{ix} = mean(ACMat, 3);
    De{ix} = mean(DeMat, 3);
    Di{ix} = mean(DiMat, 3);
    Se{ix} = mean(SeMat, 3);
    Si{ix} = mean(SiMat, 3);

    for i = 1:length(notes)
        lat{ix}(i, :) = lat0(i, :) + onset(i);
    end 

    save('figS9ix4.mat')

end

for ix = 1:length(label)

    % Psychoacoustics
    fig1 = figure;

    
    maximum = 25 * ceil(max([De{ix}(:); Di{ix}(:); Se{ix}(:)]) / 25);

    subplot(2,3,1)
    imagesc(lagSpace, notes, AC{ix})
    xlabel('regularised SACF characteristic delay (ms)');
    ylabel('stimulus period (ms)');
    caxis([0, maximum]); colormap(parula);

    subplot(2,3,2)
    imagesc(lagSpace, notes, De{ix})
    xlabel('decoder excitatory characteristic delay (ms)');
    ylabel('stimulus period (ms)');
    caxis([0, maximum]); colormap(parula);

    subplot(2,3,3)
    imagesc(lagSpace, notes, Di{ix})
    xlabel('decoder inhibitory characteristic delay (ms)');
    ylabel('stimulus period (ms)');
    caxis([0, maximum]); colormap(parula);

    subplot(2,3,4)
    imagesc(lagSpace, notes, Se{ix})
    xlabel('sustainer excitatory characteristic delay (ms)');
    ylabel('stimulus period (ms)');
    caxis([0, maximum]); colormap(parula);

    subplot(2,3,5)
    imagesc(lagSpace, notes, Si{ix})
    caxis([0, maximum]); colormap(parula);
    c = colorbar(); ylabel(c, 'average population activity (Hz)');

    fig1.PaperPosition = [0 0 10 6];
    print(fig1, sprintf('psyDyads%s.svg', label{ix}), '-dsvg');

end

% POR predictions
fig2 = figure;

% Psychophysics (Bidelman 2014)
notes = 0:12;
dAvg = [.409,1.00,.730,.573,.617,.540,.730,.461, ...
                                            .695,.557,.756,.895,.401];
dAvg = 1 - dAvg + (1 - max(1 - dAvg));
dErr = [.093,.042,.064,.051,.065,.043,.065,.081, ...
                                            .062,.065,.074,.069,.073];

for ix = 1:length(label)

    latAvg = (mean(lat{ix}, 2))';
    latErr = (std(lat{ix}, 0, 2)   / sqrt(N))';

    subplot(3, 2, ix)
    hold off; 
    ax = plotyy(notes, latAvg, notes, dAvg);
    hold(ax(1), 'on')
    ax(1).ColorOrderIndex = 1;
    errorbar(ax(1), notes, latAvg, latErr, 'lineWidth', 2);
    set(ax(1), 'YTick', 145:5:180)
    hold(ax(2), 'on')
    ax(2).ColorOrderIndex = 1;
    errorbar(ax(2), notes, dAvg, dErr);
    set(ax(2), 'YTick', 0.3:0.1:1)
    ax(2).YDir = 'reverse';
    ylabel(ax(1), 'POR latency (ms)')
    ylabel(ax(2), 'perceived consonance (au)'), 
    xlim(ax(1), [-0.5, 12.5])
    xlim(ax(2), [-0.5, 12.5])  
    ylim(ax(1), [140, 185])    %%  These values should be adjusted 
    ylim(ax(2), [0.2, 1.11])   %%  independently for each stim type!
    set(gca, 'XTick', 0:12);
    set(gca, 'XTickLabel', dyads);
    set(gca, 'XTickLabelRotation', 45);
    title(label{ix})

end

fig2.PaperPosition = [0 0 10 7];
print(fig2, 'latsOtherDyads.svg', '-dsvg');


% Stats
for ix = 1:length(label)
    [rCons, pCons] = corrcoef(dAvg, mean(lat{ix}'));
    fprintf('Corr with consonance type %s: ', label{ix});
    fprintf('R = %.3f, p = %.4f\n\n\n', rCons(1,2), pCons(1,2));
end]