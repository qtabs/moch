
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

for t = 1:size(s.p.He, 1)

    subplot(6,3,[1,2,4,5,7,8]), hold off;
    plot(lagSpace, [s.p.He(t, :); s.p.Hi(t, :); r.A(t, :)]);
    ylabel('population firing rate (Hz)')
    % xlabel('characteristic lag (ms)')
    legend('decoder excitatory', 'decoder inhibitory', ...
           'peridocity detectors')
    xlim([0.5, 30])
    ylim([0, 123])
    text(1.1, 118, 'a)', 'FontSize', 14)
    
    text(2.5, 117.5, sprintf('time = %dms', t), 'FontSize', 12);

    subplot(6,3,[10, 11, 13, 14, 16, 17]), hold off;
    plot(lagSpace, [s.q.He(t, :); s.q.Hi(t, :)]);
    ylabel('population firing rate (Hz)')
    xlabel('characteristic lag (ms)')
    legend('sustainer excitatory', 'sustainer inhibitory')
    xlim([0.5, 30])
    ylim([0, 83])
    text(1.1, 79, 'b)', 'FontSize', 14)

    subplot(6,3,[3,6]), hold off;
    h = zeros(4, 1);
    h(1) = plot(NaN, NaN, 'o', 'Color', c(1, :)); hold on;
    h(2) = plot(NaN, NaN, '.', 'Color', c(2, :)); 
    h(3) = plot(NaN, NaN, '.', 'Color', c(3, :));
    h(4) = plot(NaN, NaN, '.', 'Color', c(4, :));

    legend(h, 'before cortical onset', ...
              'between onset and convergence', ...
              'between convergence and offset', 'after offset', ...
              'AutoUpdate','off');
    legend('Location','southeast')

    if t < t1       
        plot(Y1(1:t, 2), Y1(1:t, 1), 'o', 'Color', c(1, :))
    elseif t < t2
        plot(Y1(1:t1, 2), Y1(1:t1, 1), 'o', 'Color', c(1, :))
        plot(Y1(t1:t, 2), Y1(t1:t, 1), '.', 'Color', c(2, :))
    elseif t < t3
        plot(Y1(1:t1, 2), Y1(1:t1, 1), 'o', 'Color', c(1, :))
        plot(Y1(t1:t2, 2), Y1(t1:t2, 1), '.', 'Color', c(2, :))
        plot(Y1(t2:t, 2), Y1(t2:t, 1), '.', 'Color', c(3, :))
    else
        plot(Y1(1:t1, 2), Y1(1:t1, 1), 'o', 'Color', c(1, :))
        plot(Y1(t1:t2, 2), Y1(t1:t2, 1), '.', 'Color', c(2, :))
        plot(Y1(t2:t3, 2), Y1(t2:t3, 1), '.', 'Color', c(3, :))
        plot(Y1(t3:t, 2), Y1(t3:t, 1), '.', 'Color', c(4, :))
    end
    xlabel('decoder variables - PCA 1');
    ylabel('decoder variables - PCA 2');
    xlim([-35, 90]) 
    ylim([-10, 170])
    text(-30, 155, 'c)', 'FontSize', 14)


    subplot(6,3,[9, 12]), hold off;
    if t < t1               
        plot(Y2(1:t, 2), Y2(1:t, 1), 'o', 'Color', c(1, :));
    elseif t < t2
        plot(Y2(1:t1, 2), Y2(1:t1, 1), 'o', 'Color', c(1, :)); hold on
        plot(Y2(t1:t, 2), Y2(t1:t, 1), '.', 'Color', c(2, :));
    elseif t < t3
        plot(Y2(1:t1, 2), Y2(1:t1, 1), 'o', 'Color', c(1, :)); hold on
        plot(Y2(t1:t2, 2), Y2(t1:t2, 1), '.', 'Color', c(2, :));
        plot(Y2(t2:t, 2), Y2(t2:t, 1), '.', 'Color', c(3, :));
    else
        plot(Y2(1:t1, 2), Y2(1:t1, 1), 'o', 'Color', c(1, :)); hold on
        plot(Y2(t1:t2, 2), Y2(t1:t2, 1), '.', 'Color', c(2, :));
        plot(Y2(t2:t3, 2), Y2(t2:t3, 1), '.', 'Color', c(3, :));
        plot(Y2(t3:t, 2), Y2(t3:t, 1), '.', 'Color', c(4, :)); hold on
    end
    xlabel('sustainer inhibitory (Hz)');
    ylabel('sustainer excitatory (Hz)');
    xlim([-2, 63])
    ylim([-5, 85])
    text(-0, 77, 'd)', 'FontSize', 14)
    
    subplot(6,3,[15, 18]), 
    hold off; plot(timeSpace(1:t), mean(s.p.He(1:t, :), 2), 'lineWidth', 2);
    hold on;  plot(timeSpace(1:t), mean(s.q.He(1:t, :), 2), 'lineWidth', 2);
    ylabel('average excitatory firing rate (Hz)');
    xlabel('time (ms)');
    legend('decoder', 'sustainer');
    xlim([0, 750])
    ylim([0, 3.8])
    text(12, 3.5,'e)', 'FontSize', 14)

    F(t) = getframe(fig);

end

vid = VideoWriter('decodingDynamicsIRN.avi');
vid.Quality = 100;
vid.FrameRate = 15;
open(vid);
writeVideo(vid, F);
close(vid);
