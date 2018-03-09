% Fig 1 -- Model's diagram

delays = 4:4:12;

pars = loadParameters();
pars.est.dur      = 300;
pars.est.type     = 'IRN';
pars.est.nOfIts   = 16;
pars.est.noiseOff = 0;
pars.est.bandpass = [800, 3200];  
pars.sigma = 0;  % No cortical noise for the examples -> clearer plots


for i = 1:length(delays)
    disp(i)
    pars.est.f = 1000 ./ delays(i);
    [s, r, lagSpace, timeSpace] = tdoch(pars);
    He{i} = mean(s.q.He(251:end, :));
    Ac{i} = mean(r.A(251:end, :));
    por(:, i) = mean(s.p.He, 2);
end

fig0 = figure;
l = [4.1, 3.1, 2.9];
for i = 1:length(delays)
    subplot(1, length(delays), i)
    plot((1:length(por(:, i))), por(:, i), 'k')
    ylabel('excitatory activity in the decoder network (Hz)')
    xlabel('time after tone onset (ms)')
    xlim([0, 300])
    ylim([0, l(i)])
    set(gca,'Ydir','reverse')
end

fig0.PaperPosition = [0 0 20 2]; 
print(fig0, 'fig1-0.svg', '-dsvg');
close(fig0);

fig = figure;   

for i = 1:length(delays)
    subplot(4, 3, i)
    plot(lagSpace, He{i}, 'k')
    title(sprintf('perceived pitch: %.0fms', 1000 / delays(i)))
    xlabel('characteristic period of the population (ms)')
    ylabel('average firing rate (Hz)')
    xlim([0, 30])
    ylim([0, 80])
end

subplot(4, 3, 6)
imagesc(lagSpace, lagSpace, pars.Cei)
title('decoder exc-to-inh')
xticks(5:5:30)
colormap(parula)
colorbar()
caxis([0, 1])

subplot(4, 3, 9)
imagesc(lagSpace, lagSpace, pars.Cie)
title('decoder inh-to-exc')
xticks(5:5:30)
colorbar()
colormap(parula)
caxis([0, 1])

for i = 1:length(delays)
    subplot(4, 3, 9 + i)
    plot(lagSpace, Ac{i}, 'k')
    xlabel('characteristic period of the population (ms)')
    ylabel('average firing rate (Hz)')
    xlim([0, 30])
    ylim([0, 80])
end

subplot(4, 3, [6 9])
imagesc(lagSpace, lagSpace, pars.Cie)
title('decoder inh-to-exc')
colorbar()
colormap(parula)
caxis([0, 1])
fig.PaperPosition = [0 0 8 11]; 
print(fig, 'fig1-1.svg', '-dsvg');
close(fig);