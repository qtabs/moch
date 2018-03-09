freq = [100, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 1000];

N = 10;  
t0 = tic;
parbase = loadParameters();
parbase.est.type     = 'PT';
parbase.est.bandpass = false;
parbase.est.dur      = 100;

for i = 1:length(freq)
    pars{i}       = parbase;
    pars{i}.est.f = freq(i);
end

for n = 1:N
    tt = tic;
    fprintf(' - %d of %d\n', n, N);
    parfor i = 1:length(freq)
        fprintf(' --- %d of %d\n', i, length(freq));
        [s{i}{n}, r{i}{n}] = tdoch(pars{i});
    end
    fprintf(' ...done! time: %.0fm\n', toc(tt)/60);
end


clear latAvgK latErrK; lat = 0; 
for i = 1:length(freq)
    for n = 1:N
        [~, lat(n)] = max(mean(s{i}{n}.p.He, 2)); 
    end
    latAvg(i) = mean(lat);
    latErr(i) = std(lat) / sqrt(N);
end

fig = figure;
%subplot(121)
%hold off; errorbar(freq, latAvg, latErr);
latL = [122,114,109.5,108,106,108,107.5,106,104,102.5,104,100.5,101.5];
errL = [3.0, 4.0, 4.5, 4.5, 3.5, 4.0, 5.5, 1.5, 4.5, 3.0, 3.0, 2.0, 3.5];
latR = [116.5,109,108,109.5,104,103.5,103.5,104,100.5,101.5,101,101,99.5];
errR = [4.5, 2.5, 4.0, 4.5, 2.0, 2.5, 2.5, 3.0, 2.5, 1.5, 2.0, 1.5, 2.0];
%hold on, 
%errorbar(freq, latL, errL);
%set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), 1))
%errorbar(freq, latR, errR, '--');
%legend('model POR prediction', ...
%       'experimental N100 (left)', 'experimental N100 (right)')
%xlabel('stimulus frequency (Hz)')
%ylabel('latency (ms)')
%xlim([75, 1025])
%ylim([95, 165])

%subplot(122)
latN100 = 0.5 * (latAvg + 100);
hold off; errorbar(freq, latL, errL); hold on
set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), 1))
errorbar(freq, latR, errR, '--');
errorbar(freq, latN100, 0.5 * latErr, 'k');
legend('experimental N100 (left)', 'experimental N100 (right)', ...
       'model N100 prediction')
xlabel('stimulus frequency (Hz)')
ylabel('latency (ms)')
xlim([75, 1025])
ylim([96, 127])

fig.PaperPosition = [0 0 6 4]; 
print(fig, sprintf('PTlats.svg', j), '-dsvg');
