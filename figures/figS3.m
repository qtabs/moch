N      = 10;  
delays = 4:16;

% Parameters
for j = 1:4 
    parbase{j} = loadParameters();
    parbase{j}.est.dur  = 200;
    parbase{j}.est.bandpass = [125, 2000];
    parbase{j}.est.type = 'IRN';
    onset = parbase{j}.subDelay;
    parbase{j}.subDelay = 0;
end

parbase{1}.est.its = 8;
label{1} = '4 iterations';

parbase{2}.est.its = 16;
label{2} = '8 iterations';

parbase{3}.est.its = 32;
label{3} = '32 iterations';

parbase{4}.est.its = 64;
label{4} = '64 iterations';

for j = 1:length(label)

    TT = tic;
    fprintf('Family %d/%d\n', j, length(label));

    for i = 1:length(delays)
        pars{i} = parbase{j};
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
            [~, latMat{j}(n, i)] = max(mean(s{n}.p.He, 2)); 
        end

        fprintf(' time: %.0fs\n', toc(tt));
    end

    latAvg{j} = mean(latMat{j});
    latErr{j} = std(latMat{j});

    TT = toc(TT) * (length(label) - j) / 3600;
    fprintf('... done! Time left: %.1fh\n\n', TT);

end

fig = figure;

for j = 1:4
    subplot(2,2,j)
    hold off;
    errorbar(delays, latAvg{j} + onset, latErr{j});
    hold on
    krumbDelays = [64 32 16 8, 4];
    krumbLats = [362.5, 253.5, 164.0, 140.5, 128.5];
    krumbErrs = [11.5, 10.5, 10.0, 10.5, 10.5];
    xlim([min(delays) - 0.5, max(delays) + 0.5]);
    ylim([120, 182])
    xlabel('IRN delay (ms)')
    ylabel('latency (ms)')
    title(label{j})
end

fig.PaperPosition = [0 0 10 6];
print(fig, sprintf('IRNextraLats1252000.svg', j), '-dsvg');

% Stats
for j = 1:4
    fprintf('Parametrisation: %s\n', label{j});
    [P, S] = polyfit(delays, latAvg{j}, 1);
    ste = sqrt(diag(inv(S.R) * inv(S.R')) .* S.normr .^ 2 ./ S.df);
    fitRes = [P(1), ste(1); P(2), ste(2)];
    fprintf('Lat fit with delay: m = %.1f%s%.1f, ', P(1), char(177), ste(1));
    fprintf('n = %.1f%s%.1f\n\n', P(2), char(177), ste(2));
end