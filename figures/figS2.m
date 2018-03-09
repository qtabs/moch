pars = loadParameters();
pars.est.bandpass  = [800, 3200];
pars.est.type      = 'IRNseq';
pars.est.notes     = [0, 2];
pars.est.nOfIts    = 16; 
pars.est.f         = 200;
pars.est.dur       = 300;
pars.est.tail      = 250; 
pars.sigma         = 0;

[s, r, lagSpace, timeSpace] = tdoch(pars);

fig = figure;
L = 25*ceil(max([r.A(:); s.p.He(:); s.p.Hi(:); s.q.He(:); s.q.Hi(:)])/25);

subplot(5, 4, [ 1  2]), imagesc(timeSpace, lagSpace, r.A');    
ylim([1, 15]); ylabel({'periodicity'; 'detectors'}); 
colormap(parula); caxis([0, L])
subplot(5, 4, [ 5  6]), imagesc(timeSpace, lagSpace, s.p.He'); 
ylim([1, 15]); ylabel('excitatory'); 
colormap(parula); caxis([0, L])
subplot(5, 4, [9 10]), imagesc(timeSpace, lagSpace, s.p.Hi'); 
ylim([1, 15]); ylabel('inhbitory'); 
colormap(parula); caxis([0, L])
subplot(5, 4, [13 14]), imagesc(timeSpace, lagSpace, s.q.He'); 
ylim([1, 15]); ylabel('excitatory'); 
colormap(parula); caxis([0, L])
subplot(5, 4, [17 18]), imagesc(timeSpace, lagSpace, s.q.Hi'); 
ylim([1, 15]); ylabel('inhibitory'); 
xlabel('time (ms)'); xlabel('characteristic lag (ms)')
colormap(parula); caxis([0, L])


X1 = [s.p.He];

[coeff, score, lat, tsq, percentage] = pca(X1);
Y1 = X1 * coeff(:, 1:10);

Z1 = Y1(:, 2);
Z2 = Y1(1:t1, 1) + Y1(1:t1, 3);

[~, ind0] = min( (lagSpace - 1000/pars.est.f).^2 );
[~, ind1] = min( (lagSpace - 1000/(9 * pars.est.f / 8)).^2 );

Y2 = [s.q.He(:, ind0), s.q.Hi(:, ind0)];
Y3 = [s.q.He(:, ind1), s.q.Hi(:, ind1)];


Z1 = Y1(:, 2) + 0.3 * Y1(:, 3);
Z2 = Y1(:, 1);

t1 = 50; t2 = 85+50; t3 = 300+50; t4 = 435 + 50; t5 = 600 + 50;
subplot(5, 4, [3 4 7 8])
hold off
plot(Z1(1:t1), Z2(1:t1), 'o')
hold on
plot(Z1(t1:t2), Z2(t1:t2), '.')
plot(Z1(t2:t3), Z2(t2:t3), '.')
plot(Z1(t3:t4), Z2(t3:t4), '.')
plot(Z1(t4:t5), Z2(t4:t5), '.')
plot(Z1(t5:end), Z2(t5:end), '.')
xlabel('decoder variables - PCA 1');
ylabel('decoder variables - PCA 2');
xlim([-5, 115])
ylim([-105, 130])
legend('before onset', 'between onset and convergence', ...
       'between convergence and pitch change', ...
       'between pitch change and new convergence', ...
       'between convergence and offset', 'after offset');

subplot(5, 4, [11 12])
hold off
plot(Y2(1:t1, 2), Y2(1:t1, 1), 'o')
hold on
plot(Y2(t1:t2, 2), Y2(t1:t2, 1), '.')
plot(Y2(t2:t3, 2), Y2(t2:t3, 1), '.')
plot(Y2(t3:t4, 2), Y2(t3:t4, 1), '.')
plot(Y2(t4:t5, 2), Y2(t4:t5, 1), '.')
plot(Y2(t5:end, 2), Y2(t5:end, 1), '.')
ylabel('sustainer excitatory (Hz)');
xlim([-2, 63])
ylim([-5, 90])

subplot(5, 4, [15 16])
hold off
plot(Y3(1:t1, 2), Y3(1:t1, 1), 'o')
hold on
plot(Y3(t1:t2, 2), Y3(t1:t2, 1), '.')
plot(Y3(t2:t3, 2), Y3(t2:t3, 1), '.')
plot(Y3(t3:t4, 2), Y3(t3:t4, 1), '.')
plot(Y3(t4:t5, 2), Y3(t4:t5, 1), '.')
plot(Y3(t5:end, 2), Y3(t5:end, 1), '.')
xlabel('sustainer inhibitory (Hz)');
xlim([-2, 63])
ylim([-5, 90])

fig.PaperPosition = [0 0 10 10];
print(fig, sprintf('pitchChage.svg', j), '-dsvg');

subplot(5, 4, [1 2 5 6 9 10 13 14]); imagesc(timeSpace, lagSpace, r.A'); 
caxis([0, L]); 
colb = colorbar(); ylabel(colb, 'population activity (Hz)')
fig.PaperPosition = [0 0 5 10];
print(fig, sprintf('pitchChageBar.svg', j), '-dsvg');