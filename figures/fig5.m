% Fig 5 - Results with dyads

N  = 60;
notes    = [0:1:12];
dyads    = {'unison', 'minor second', 'second', 'minor third', ...
            'third', 'fourth', 'tritone', 'perfect fifth', ...
            'minor sixth', 'sixth', 'minor seventh', 'seventh', 'octave'};
its      = 8;
dur      = 200;
f0       = 160;
bandpass = [125, 2000];
tuning   = 'just';
stimType = 'IRNchordSS';

parbase  = loadParameters();
onset    = repmat(parbase.subDelayDy, size(notes));
onset(1) = parbase.subDelay;
parbase.subDelayDy   = 0;
parbase.subDelay     = 0;
parbase.est.dur      = dur; 
parbase.est.type     = stimType;
parbase.est.f        = f0;
parbase.est.nOfIts   = its;
parbase.est.bandpass = bandpass;
parbase.est.tuning   = tuning;

for i = 1:length(notes)
    pars{i} = parbase;
    pars{i}.est.notes = [0, notes(i)]; 
end

[~, r] = tdoch(pars{1}]);
lagSpace  = r.lagSpace;
timeSpace = r.timeSpace;

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
    fprintf('done! time left: %.0fm\n', (length(notes)-i) * toc(tt)/60);
end 

for n = 1:N
    ACMat(:, :, n) = ACPar{n};
    DeMat(:, :, n) = DePar{n};
    DiMat(:, :, n) = DiPar{n};
    SeMat(:, :, n) = SePar{n};
    SiMat(:, :, n) = SiPar{n};
    lat0(:, n)     = latPar{n};
end

% Psychoacoustics
fig1 = figure;

AC = mean(ACMat, 3);
De = mean(DeMat, 3);
Di = mean(DiMat, 3);
Se = mean(SeMat, 3);
Si = mean(SiMat, 3);

maximum = 25 * ceil(max([De(:); Di(:); Se(:)]) / 25);

subplot(2,3,1)
imagesc(lagSpace,notes, AC)
xlabel('regularised SACF characteristic delay (ms)');
ylabel('stimulus period (ms)');
caxis([0, maximum]); colormap(parula);

subplot(2,3,2)
imagesc(lagSpace, notes, De)
xlabel('decoder excitatory characteristic delay (ms)');
ylabel('stimulus period (ms)');
caxis([0, maximum]); colormap(parula);

subplot(2,3,3)
imagesc(lagSpace, notes, Di)
xlabel('decoder inhibitory characteristic delay (ms)');
ylabel('stimulus period (ms)');
caxis([0, maximum]); colormap(parula);

subplot(2,3,4)
imagesc(lagSpace, notes, Se)
xlabel('sustainer excitatory characteristic delay (ms)');
ylabel('stimulus period (ms)');
caxis([0, maximum]); colormap(parula);

subplot(2,3,5)
imagesc(lagSpace, notes, Si)
caxis([0, maximum]); colormap(parula);
c = colorbar(); ylabel(c, 'average population activity (Hz)');

fig1.PaperPosition = [0 0 10 6];
print(fig1, sprintf('fig5-0.svg', j), '-dsvg');



% POR predictions
fig2 = figure;

for i = 1:length(notes)
    lat(i, :) = lat0(i, :) + onset(i);
end 

latAvg = (mean(lat, 2))';
latErr = (std(lat, 0, 2)   / sqrt(N))';

% MEG fields
datapath = '~/Cloud/Projects/TDoCh/Doc/figs/data/';
aefs     = load([datapath, 'aefRes.mat']);
n1Lats   = load([datapath, 'n1Lat.mat']);
eLat     = n1Lats.n1Lat;
eNotes   = aefs.notes;
eLatAvg  = aefs.n1LatAvg;
eLatErr  = aefs.n1LatErr;

subplot(121)
hold off; 
errorbar(eNotes, latAvg(eNotes+1), latErr(eNotes+1));
hold on;
errorbar(eNotes, eLatAvg, eLatErr)
set(gca, 'XTick', eNotes);
set(gca, 'YTick', 145:5:180)
for i = 1:length(eNotes), eDyads{i} = dyads{eNotes(i) + 1}; end;
set(gca, 'XTickLabel', eDyads);
set(gca, 'XTickLabelRotation', 45);
xlim([-0.5, 10.5]);
ylim([140, 185])
ylabel('POR latency (ms)')
legend('model predictions', 'observer latency')

subplot(122)
hold off; 
errorbar(notes, latAvg, latErr);
set('YTick', 145:5:180)
ylabel('POR latency (ms)')
xlim([-0.5, 12.5])
ylim([140, 185])
set(gca, 'XTick', 0:12);
set(gca, 'XTickLabel', dyads);
set(gca, 'XTickLabelRotation', 45);

fig2.PaperPosition = [0 0 10 3];
print(fig2, 'fig5-1.svg', '-dsvg');


% Stats
% 1. Correlation between latency predictions and observations
[r, p] = corrcoef(eLat, lat(eNotes + 1, 1:size(eLat, 2)));
fprintf('Corr with MEG data: R = %.3f, p = %.4f\n\n', r(1,2), p(1,2));

% 2. p-values for the differencial response
diss = 2:2:6; cons = 1:2:5;
for i = 1:length(diss)
    for j = 1:length(cons)
        ix = eNotes(diss(i)) + 1; jx = eNotes(cons(j)) + 1;
        [pExp(i, j), ~, uStruct] = ranksum(eLat(diss(i), :), ...
                                                       eLat(cons(j), :));
        uExp(i, j) = uStruct.ranksum;
        [pMod(i, j), ~, uStruct] = ranksum(lat(ix, :), ...
                                           lat(jx, :),  'tail', 'right');
        uMod(i, j) = uStruct.ranksum;
    end
end

fprintf('Experimental ranksum p-values for dissonance <> consonance:\n');
fprintf('    |    P1    |    M3    |    P5    |\n');
fprintf(' m2 | %.6f | %.6f | %.6f |\n'  , pExp(1, :));
fprintf(' TT | %.6f | %.6f | %.6f |\n'  , pExp(2, :));
fprintf(' m7 | %.6f | %.6f | %.6f |\n\n', pExp(3, :));

fprintf('Experimental ranksum U-values for dissonance > consonance:\n');
fprintf('    |   P1   |   M3   |   P5   |\n');
fprintf(' m2 | %6.0f | %6.0f | %6.0f |\n'  , uExp(1, :));
fprintf(' TT | %6.0f | %6.0f | %6.0f |\n'  , uExp(2, :));
fprintf(' m7 | %6.0f | %6.0f | %6.0f |\n\n', uExp(3, :));

fprintf('Model (1-tail) ranksum p-values for dissonance > consonance:\n');
fprintf('    |    P1    |    M3    |    P5    |\n');
fprintf(' m2 | %.6f | %.6f | %.6f |\n'  , pMod(1, :));
fprintf(' TT | %.6f | %.6f | %.6f |\n'  , pMod(2, :));
fprintf(' m7 | %.6f | %.6f | %.6f |\n\n', pMod(3, :));

fprintf('Model (1-tail) ranksum U-values for dissonance > consonance:\n');
fprintf('    |   P1   |   M3   |   P5   |\n');
fprintf(' m2 | %6.0f | %6.0f | %6.0f |\n'  , uMod(1, :));
fprintf(' TT | %6.0f | %6.0f | %6.0f |\n'  , uMod(2, :));
fprintf(' m7 | %6.0f | %6.0f | %6.0f |\n\n', uMod(3, :));

% 2. p-values for the extended dyads
diss = [1 2 6 10 11] + 1; cons = [0 4 5 7 12] + 1;
for i = 1:length(diss)
    for j = 1:length(cons)
        [pExt(i, j), ~, uStruct] = ranksum(lat(diss(i), :), ...
                                       lat(cons(j), :), 'tail', 'right');
        uExt(i, j) = uStruct.ranksum;
    end
end

fprintf('Ext (1-tail) ranksum p-values for dissonance > consonance:\n');
fprintf('    |    P1    |    M3    |    P4    |    P5    |    P8    |\n');
fprintf(' m2 | %.6f | %.6f | %.6f | %.6f | %.6f |\n'  , pExt(1, :));
fprintf(' M2 | %.6f | %.6f | %.6f | %.6f | %.6f |\n'  , pExt(2, :));
fprintf(' TT | %.6f | %.6f | %.6f | %.6f | %.6f |\n'  , pExt(3, :));
fprintf(' m7 | %.6f | %.6f | %.6f | %.6f | %.6f |\n'  , pExt(4, :));
fprintf(' M7 | %.6f | %.6f | %.6f | %.6f | %.6f |\n\n', pExt(5, :));

fprintf('Ext (1-tail) ranksum U-values for dissonance > consonance:\n');
fprintf('    |   P1   |   M3   |   P4   |   P5   |   P8   |\n');
fprintf(' m2 | %6.0f | %6.0f | %6.0f | %6.0f | %6.0f |\n'  , uExt(1, :));
fprintf(' M2 | %6.0f | %6.0f | %6.0f | %6.0f | %6.0f |\n'  , uExt(2, :));
fprintf(' TT | %6.0f | %6.0f | %6.0f | %6.0f | %6.0f |\n'  , uExt(3, :));
fprintf(' m7 | %6.0f | %6.0f | %6.0f | %6.0f | %6.0f |\n'  , uExt(4, :));
fprintf(' M7 | %6.0f | %6.0f | %6.0f | %6.0f | %6.0f |\n\n', uExt(5, :));



% Comparison of the simulated dipole moments with MEG data
clear

datapath   = '../data/';
aefs       = load([datapath, 'aefRes.mat']);
notes      = aefs.notes;
dyadDefs   = {'unison', 'minor second', 'second', 'minor third', ...
              'third', 'fourth', 'tritone', 'perfect fifth', ...
           'minor sixth', 'sixth', 'minor seventh', 'seventh', 'octave'};

for i = 1:length(notes), dyads{i} = dyadDefs(notes(i) + 1); end;

N                    = 60;
parbase              = loadParameters();
parbase.est.dur      = 275; 
parbase.est.type     = 'IRNchordSS';
parbase.est.f        = 160;
parbase.est.nOfIts   = 8;
parbase.est.bandpass = [125, 2000];
parbase.est.tuning   = 'just';

for i = 1:length(notes)
    pars{i} = parbase;
    pars{i}.est.notes = [0, notes(i)]; 
end

pars{1}.est.tail = parbase.subDelayDy - parbase.subDelay;

for i = 1:length(notes)
    tt = tic;
    fprintf(' - %d of %d ...\n', i, length(notes));
    parfor n = 1:N
        disp(n)
        [s, r] = tdoch(pars{i});
        spHe(:, n) = mean(s.p.He, 2);
    end
    modPor{i} = mean(spHe, 2)';
    modErr{i} = std(spHe, 0, 2)';
    fprintf('done! time left: %.0fm\n', (length(notes)-i) * toc(tt)/60);
end 

[~, r] = tdoch(pars{1});

timeSpace = r.timeSpace;
megOnset = 1256;

for i = 1:length(notes)
    megPor{i} = aefs.n1fAvg(i, megOnset:(megOnset + length(timeSpace) - 1));
    megErr{i} = aefs.n1fErr(i, megOnset:(megOnset + length(timeSpace) - 1));
end

addpath(genpath('~/Apps/matlab/boundedline'));
close all
fig = figure;
for i = 1:6
    subplot(2, 3, i)
    hold off; 
    ax = plotyy(timeSpace, modPor{i}, timeSpace, megPor{i});
    hold(ax(2), 'on');
    boundedline(timeSpace, modPor{i}, modErr{i}, ax(1), 'b');
    boundedline(timeSpace, megPor{i}, megErr{i}, ax(2), 'r');
    ax(1).YDir = 'reverse';
    ylabel(ax(1), 'decoder avg activity (Hz)'), 
    ylabel(ax(2), 'dipole moment (nAm)')
    ax(1).XTick = 0:50:300;
    ax(2).XTick = 0:50:300;
    ax(1).YTick = 0:1:6;
    ax(2).YTick = -30:5:5;
    xlim(ax(1), [0, 350]);
    ylim(ax(1), [-1, 6.5]);
    xlim(ax(2), [0, 350]);
    ylim(ax(2), [-27, 5]);
    if i == 1, ylim(ax(1), [-0.2, 3]); ylim(ax(2), [-32, 5]); end
    title(dyads{i})
end

legend('experimental data', 'simulated fields');

fig.PaperPosition = [0 0 10 4];
print(fig, sprintf('fig5-2.svg', j), '-dsvg');

