% Fig 2 -- Details of the decoding dynamics

pars = loadParameters();
pars.est.bandpass = [800, 3200];  
pars.est.type     = 'IRN';
pars.est.f        = 200;   
pars.est.dur      = 400;   
pars.est.nOfIts   = 16;    

[s, r, lagSpace, timeSpace] = tdoch(pars);
sSamp(1, :) = mean(s.p.He, 2);

parfor i = 2:5
    s0 = tdochCortex(r, pars);
    sSamp(i, :) = mean(s0.p.He, 2);
end

pars.sigma = 0;
sAvg = tdochCortex(r, pars);

fig = figure;
xlims = [1, 250];

subplot(6,4, 1:3 ), 
imagesc(timeSpace, lagSpace, r.A');    
ylabel({'subcortical'; 'characteristic period (ms)'}); 
colormap(parula); 
xlim(xlims)
ylim([1, 20]); 
caxis([0, 150])

subplot(6,4, 5:7 ), 
imagesc(timeSpace, lagSpace, s.p.He'); 
ylabel({'decoder';'excitatory'; 'characteristic period (ms)'}); 
colormap(parula); 
xlim(xlims)
ylim([1, 20]); 
caxis([0, 150])

subplot(6,4, 9:11 ), 
imagesc(timeSpace, lagSpace, s.p.Hi'); 
ylabel({'decoder';'inhibitory'; 'characteristic period (ms)'}); 
colormap(parula); 
xlim(xlims)
ylim([1, 20]); 
caxis([0, 150])

subplot(6,4, 13:15 ), 
imagesc(timeSpace, lagSpace, s.q.He'); 
ylabel({'sustainer';'excitatory'; 'characteristic period (ms)'}); 
colormap(parula); 
xlim(xlims)
ylim([1, 20]); 
caxis([0, 150])

subplot(6,4, 17:19 ), 
imagesc(timeSpace, lagSpace, s.q.Hi'); 
xlabel('time (ms)');
ylabel({'sustainer';'inhibitory'; 'characteristic period (ms)'}); 
colormap(parula); 
xlim(xlims)
ylim([1, 20]); 
caxis([0, 150])

subplot(6,4, 4*(1:4)); 
imagesc(timeSpace, lagSpace, r.A'); 
colormap(parula);
caxis([0, 150]);  
colb = colorbar(); 
ylabel(colb, 'population activity (Hz)')

subplot(6,4, 21:23)
hold off;
plot(timeSpace, -mean(sSamp), 'k'); ax = gca;
hold on; 
for i = 1:5,
    ax.ColorOrderIndex = 1; 
    plot(timeSpace, -sSamp(i, :)); 
    ylim([-5, 0])
    xlim(xlims)
end;
hold off;
ylabel({'average decoder'; 'excitatory activity (Hz)'})
xlabel('time (ms)');
legend('average response', 'response in a single trial')

fig.PaperPosition = [0 0 12 12]; 
print(fig, sprintf('fig2.svg', j), '-dsvg');