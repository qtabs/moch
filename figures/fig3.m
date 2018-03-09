% Fig 3 - Results with simple IRN

N      = 60;  
delays = 4:16;
its    = [2 4 8 16 32];

% Parameters
parbase = loadParameters();
onset = parbase.subDelay; 

% Fixed nOfIts, varying pitch
parbase.est.bandpass = [800, 3200];
parbase.est.nOfIts   = 16;
parbase.est.dur      = 200;
parbase.est.noiseOff = 0;
parbase.est.type     = 'IRN';
parbase.subDelay     = 0;

[~, ~, lagSpace, timeSpace] = tdoch(parbase);

for i = 1:length(delays)
    pars{i} = parbase;
    pars{i}.est.f = 1000 / delays(i);
end

% Coputing IRN stuff
for i = 1:length(delays)
    tt = tic;
    fprintf(' - %d of %d ...', i, length(delays));
    parfor n = 1:N
        [s{n}, r{n}] = tdoch(pars{i});
    end
    for n = 1:N
        [~, latMat(n,i)] = max(mean(s{n}.p.He, 2));
    end
    expS{i} = s{1};
    expR{i} = r{1};
    fprintf(' time: %.0fs\n', toc(tt));
end

latAvg = mean(latMat);
latErr = std(latMat);


% Fixed pitch, varying its 
parbase2 = loadParameters();
parbase2.est.bandpass = [800, 3200];
parbase2.est.f = 1000 / 16;
parbase2.est.noiseOff = 0;
parbase2.est.type     = 'IRN';
parbase2.est.dur      = 200;
parbase2.subDelay     = 0;

for i = 1:length(its)
    pars2{i} = parbase2;
    pars2{i}.est.nOfIts = its(i);
end

% Coputing IRN stuff

for i = 1:length(its)
    tt = tic;
    fprintf(' - %d of %d ...', i, length(its));
    parfor n = 1:N
        [s{n}, r{n}] = tdoch(pars2{i});
    end
    for n = 1:N
        [~, latMat2(n, i)] = max(mean(s{n}.p.He, 2)); 
    end

    fprintf(' time: %.0fs\n', toc(tt));
end

lat2Avg = mean(latMat2);
lat2Err = std(latMat2);


% Latency predictions
fig1 = figure;

subplot(121)
hold off
errorbar(delays, latAvg + onset, latErr);
hold on
krumbDelays = [64 32 16 8, 4];
krumbLats = [362.5, 253.5, 164.0, 140.5, 128.5];
krumbErrs = [11.5, 10.5, 10.0, 10.5, 10.5];
errorbar(krumbDelays(3:end), krumbLats(3:end), krumbErrs(3:end));
xlim([min(delays) - 0.5, max(delays) + 0.5]);
ylim([115, 185])
ylim([115, 215])
legend('predicted POR latency', 'experimental POR latency')
xlabel('IRN delay (ms)')
ylabel('latency (ms)')

subplot(122)
hold off
errorbar(1:5, lat2Avg + onset, lat2Err);
hold on
eLats = [179, 197, 185, 160, 158];
eLErr = [10, 10, 10, 10, 10];
ax = errorbar(1:5, eLats, eLErr);
ylim([115, 215])
xlim([0.5, 5.5])
set(gca, 'XTick', 1:5)
set(gca, 'XTickLabel', {'2', '4', '8', '16', '32'})

legend('predicted POR latency', 'experimental POR latency')
xlabel('number of iterations')
ylabel('latency (ms)')

fig1.PaperPosition = [0 0 10 3];
print(fig1, sprintf('fig3-0.svg', j), '-dsvg');


% Example with real data
pars = loadParameters();
pars.est.bandpass = [800, 3200];
pars.est.f        = 1000 / 8;
pars.subCortAff   = 4;
pars.nOfIts       = 16;
pars.est.noiseOff = 0;
pars.est.type     = 'IRN';
pars.sigma        = 0.0001;
pars.est.dur      = 250 - pars.subDelay;

[s, r, lagSpace, timeSpace] = tdoch(pars);

dd = load('~/Cloud/Projects/TDoCh/Data/simpleIRN/irn_exp1.mat');
expFields = mean(mean(dd.in_d8_i16_g200, 3), 2);
expFields = expFields(101:350);
expErrors = std(mean(dd.in_d8_i16_g200, 2), 0, 3);
expErrors = expErrors(101:350) / sqrt(size(dd.in_d8_i16_g200, 3));
modFields = mean(s.p.He, 2);

fig2 = figure;
ax = plotyy(timeSpace, modFields, timeSpace, expFields);
boundedline(timeSpace, expFields, expErrors, ax(2), 'b');

set(ax(1), 'Ydir','reverse')
xlim(ax(1), [0, 250]);
xlim(ax(2), [0, 250]);
ylim(ax(1), [0, 4]);
ylim(ax(2), [-30, 2]);
yticks(ax(2), [-30:10:0])
yticks(ax(1), [0:1:4])
xlabel(ax(1), 'time after tone onset (ms)');
ylabel(ax(1), 'collective exc. activation in the decoder (Hz)')
ylabel(ax(2), 'equivalen dipole moment in the POR generator (nAm)')

fig2.PaperPosition = [0 0 9 2.2];
print(fig2, sprintf('fig3-1.svg', j), '-dsvg');


% Psychoacoustics
clearvars -except lagSpace

parbase3 = loadParameters();
parbase3.est.dur = 250;
parbase3.est.its   = 16;
parbase3.est.type  = 'IRN';
parbase3.est.bandpass = [800, 3200];
parbase3.subDelay = 0;

lags     = lagSpace(1:2:end);

for i = 1:length(lags)
    pars{i} = parbase3;
    pars{i}.est.f = 1000 / lags(i);
end

parfor i = 1:length(lags)
    fprintf('-- lag %d/%d', i, length(lags))
    tt = tic;
    [s{i}, r{i}] = tdoch(pars{i});
    time = toc(tt) / 60;
    fprintf(' done! t = %.1f\n', time);
end

interval = 201:250;
for i = 1:length(lags)
    sacf(i, :) = mean(r{i}.A(interval, :), 1);
    pexc(i, :) = mean(s{i}.p.He(interval, :), 1);
    pinh(i, :) = mean(s{i}.p.Hi(interval, :), 1);
    qexc(i, :) = mean(s{i}.q.He(interval, :), 1);
    qinh(i, :) = mean(s{i}.q.Hi(interval, :), 1);
end

maxInt = 25 * ceil(max([pexc(:); pinh(:); qexc(:)]) / 25);

fig3 = figure;

subplot(231); 
imagesc(lagSpace, lags, sacf);
title('regularised SACF')
caxis([0, maxInt]); colormap(fake_parula);
xlabel('characteristic delay (ms)')
ylabel('stimulus period (ms)')

subplot(232); 
imagesc(lagSpace, lags, pexc);
title('decoder excitatory')
caxis([0, maxInt]); colormap(fake_parula);
xlabel('characteristic delay (ms)')
ylabel('stimulus period (ms)')

subplot(233); 
imagesc(lagSpace, lags, pinh);
title('decoder inhibitory')
caxis([0, maxInt]); colormap(fake_parula);
xlabel('characteristic delay (ms)')
ylabel('stimulus period (ms)')

subplot(234); 
imagesc(lagSpace, lags, qexc);
title('sustainer excitatory')
caxis([0, maxInt]); colormap(fake_parula);
xlabel('characteristic delay (ms)')
ylabel('stimulus period (ms)')

subplot(235); 
imagesc(lagSpace, lags, qinh);
title('sustainer inhibitory')
caxis([0, maxInt]); colormap(fake_parula);
xlabel('characteristic delay (ms)')

subplot(236); 
imagesc(lagSpace, lags, qinh);
caxis([0, maxInt]); colormap(fake_parula);
c = colorbar();
ylabel(c, 'average population activity (Hz)')

fig3.PaperPosition = [0 0 10 6];
print(fig3, sprintf('fig3-2.svg', j), '-dsvg');
close(fig3);