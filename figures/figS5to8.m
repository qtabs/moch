
if exist('lagSpace') == 0
    try 
        lagSpace = r{1}{1}.lagSpace;
    catch
        [~, ~, lagSpace, ~] = tdoch();
    end
end

lags = lagSpace(1:2:end);

for i = 1:length(lags)

    for j = 1:8
        pars{j}{i} = loadParameters();
        pars{j}{i}.est.f = 1000 / lags(i);
        pars{j}{i}.est.dur = 250;
        pars{j}{i}.est.bandpass = false;
        % this parameter does not affect psycoacoustic results
        % but severely increases the computational cost
        pars{j}{i}.subCortAff = 1; 
        pars{j}{i}.subDelay   = 0; 
    end
    
    pars{1}{i}.est.type  = 'PT';
    label{1} = 'pure tones';

    pars{2}{i}.est.type  = 'HCT';
    pars{2}{i}.est.harms = 0:5;
    pars{2}{i}.est.harmFact = 1;
    label{2} = 'HCTs, harmonics 0-5';

    pars{3}{i}.est.type  = 'CT';
    label{3} = 'click trains';

    pars{4}{i}.est.type  = 'HCT';
    pars{4}{i}.est.harms = 1:4;
    pars{4}{i}.est.harmFact = 1;
    label{4} = 'HCTs, harmonics 1-4';

    pars{5}{i}.est.type  = 'HCT';
    pars{5}{i}.est.harms = 0:50;
    pars{5}{i}.est.harmFact = 0.5;
    pars{5}{i}.est.bandpass = [3200, 5000];
    label{5} = 'HCTs, unresolved harmonics';

    pars{6}{i}.est.type  = 'IRN';
    pars{6}{i}.est.its   = 32;
    label{6} = 'IRNs, 32 iterations';
        
    pars{7}{i}.est.type  = 'IRN';
    pars{7}{i}.est.its   = 4;
    label{7} = 'IRNs, 4 iterations';

    pars{8}{i}.est.type  = 'IRN';
    pars{8}{i}.est.its   = 8;
    pars{8}{i}.est.bandpass = [125, 2000];
    label{8} = 'IRNs, 8 iterations, bandpass 125Hz-2kHz';

end


for j = 1:length(pars)

    fprintf('Label %d/%d\n', j, length(pars))
    rt = tic;
    pTmp = pars{j};

    parfor i = 1:length(lags)
        fprintf('-- lag %d/%d', i, length(lags))
        tt = tic;
        [sTmp{i}, rTmp{i}] = tdoch(pTmp{i});
        time = toc(tt) / 60;
        fprintf(' done! t = %.1f\n', time);
    end
    for i = 1:length(lags)
        r{j}{i} = rTmp{i};
        interval = 201:250;
        sacf{j}(i, :) = mean(r{j}{i}.A(interval, :), 1);
        pexc{j}(i, :) = mean(sTmp{i}.p.He(interval, :), 1);
        pinh{j}(i, :) = mean(sTmp{i}.p.Hi(interval, :), 1);
        qexc{j}(i, :) = mean(sTmp{i}.q.He(interval, :), 1);
        qinh{j}(i, :) = mean(sTmp{i}.q.Hi(interval, :), 1);
    end

    rt = toc(rt) / 60;
    rtl = rt * (length(label) - j);
    fprintf('Label done! Real time: %.1fm, left: %.1fm\n', rt, rtl)
    save('panicPsycho.mat', '-v7.3')

end


for j = 1:length(pars)

    fig{j} = figure;
    maximum = 25 * ceil(max([pexc{j}(:); pinh{j}(:); qexc{j}(:)]) / 25);
        
    subplot(231); 
    imagesc(lagSpace, lags, sacf{j});
    title('regularised SACF')
    caxis([0, maximum]);  colormap(parula);
    xlabel('characteristic delay (ms)')
    ylabel('stimulus period (ms)')
    title(label{j}) 

    subplot(232); 
    imagesc(lagSpace, lags, pexc{j});
    title('decoder excitatory')
    caxis([0, maximum]); colormap(parula);
    xlabel('characteristic delay (ms)')
    ylabel('stimulus period (ms)')
    
    subplot(233); 
    imagesc(lagSpace, lags, pinh{j});
    title('decoder inhibitory')
    caxis([0, maximum]); colormap(parula);
    xlabel('characteristic delay (ms)')
    ylabel('stimulus period (ms)')
    
    subplot(234); 
    imagesc(lagSpace, lags, qexc{j});
    title('sustainer excitatory')
    caxis([0, maximum]); colormap(parula);
    xlabel('characteristic delay (ms)')
    ylabel('stimulus period (ms)')
    
    subplot(235); 
    imagesc(lagSpace, lags, qinh{j});
    title('sustainer inhibitory')
    caxis([0, maximum]); colormap(parula);
    xlabel('characteristic delay (ms)')
    
    subplot(236); 
    imagesc(lagSpace, lags, qinh{j});
    title(label{j})
    caxis([0, maximum]); colormap(parula);
    c = colorbar();
    ylabel(c, 'average population activity (Hz)')

    fig{j}.PaperPosition = [0 0 10 6];
    print(fig{j}, sprintf('psy%d.svg', j), '-dsvg');
    close(fig{j});

end