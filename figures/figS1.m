
pars = loadParameters();
pars.est.bandpass  = [800, 3200];
pars.est.type      = 'IRN';
pars.est.nOfIts    = 16; 
pars.subCortAff    = 4;
pars.est.f         = 200;
pars.est.dur       = 350;
pars.est.tail      = 350; 
pars.sigma         = 0;

[s, r, lagSpace, timeSpace] = tdoch(pars);

X1 = [s.p.He];

[coeff, score, lat, tsq, percentage] = pca(X1);
Y1 = X1 * coeff(:, 1:10);

[~, ind] = min( (lagSpace - 1000/pars.est.f).^2 );
Y2 = [s.q.He(:, ind), s.q.Hi(:, ind)];

t1 = 51; t2 = 125+50; t3 = 350+50;

c = get(gca,'ColorOrder');

fig = figure;

subplot(1, 2, 1), hold off;
plot(Y1(1:t1, 2), Y1(1:t1, 1), 'o', 'Color', c(1, :)); hold on
plot(Y1(t1:t2, 2), Y1(t1:t2, 1), '.', 'Color', c(2, :))
plot(Y1(t2:t3, 2), Y1(t2:t3, 1), '.', 'Color', c(3, :))
plot(Y1(t3:end, 2), Y1(t3:end, 1), '.', 'Color', c(4, :))
xlabel('decoder variables - PCA 1');
ylabel('decoder variables - PCA 2');
xlim([-38, 90]) 
ylim([-10, 170])

subplot(1, 2, 2), hold off;
plot(Y2(1:t1, 2), Y2(1:t1, 1), 'o', 'Color', c(1, :)); hold on
plot(Y2(t1:t2, 2), Y2(t1:t2, 1), '.', 'Color', c(2, :));
plot(Y2(t2:t3, 2), Y2(t2:t3, 1), '.', 'Color', c(3, :));
plot(Y2(t3:end, 2), Y2(t3:end, 1), '.', 'Color', c(4, :)); hold on
xlabel('sustainer inhibitory (Hz)');
ylabel('sustainer excitatory (Hz)');
xlim([-2, 63])
ylim([-5, 85])

h = zeros(4, 1);
h(1) = plot(NaN, NaN, 'o', 'Color', c(1, :));
h(2) = plot(NaN, NaN, '.', 'Color', c(2, :)); 
h(3) = plot(NaN, NaN, '.', 'Color', c(3, :));
h(4) = plot(NaN, NaN, '.', 'Color', c(4, :));
legend(h, 'before cortical onset', ...
          'between onset and convergence', ...
          'between convergence and offset', 'after offset', ...
          'AutoUpdate','off');
legend('Location','northwest')

fig.PaperPosition = [0 0 10 4];
print(fig, sprintf('attractorSingleIRN.svg', j), '-dsvg');
