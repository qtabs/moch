pars = loadParameters();
pars.sigma        = 0;
pars.est.dur       = 325;
pars.est.tail      = 0; 
pars.est.type     = 'IRNchordSS';
pars.est.f        = 160;
pars.est.nOfIts   = 8;
pars.est.bandpass = [125, 2000];
pars.est.tuning   = 'just';

notes = [1, 7];

pars.est.notes = [0, notes(1)];
for i = 1:2
    pars.est.notes = [0, notes(i)];
    [s{i}, r{i}, lagSpace, timeSpace] = tdoch(pars);
    X1 = [s{i}.p.He];
    [coeff, score, lat, tsq, percentage] = pca(X1);
    Y{i} = X1 * coeff(:, 1:10);
end

c = get(gca,'ColorOrder');

fig = figure;

for t = 1:length(timeSpace)

    subplot(6,3,[1,2,4,5,7,8]), hold off;
    plot(lagSpace, [s{1}.p.He(t, :); s{1}.p.Hi(t, :); ...
                    r{1}.A(t, :); s{1}.q.He(t, :)]);
    ylabel('population firing rate (Hz)')
    xlabel('characteristic lag (ms)')
    legend('decoder excitatory', 'decoder inhibitory', ...
           'peridocity detectors', 'sustainer excitatory')
    xlim([0.5, 30])
    ylim([0, 105]) 
    text(1, 99, 'a) minor second', 'FontSize', 14)
    text(6.5, 98.5, sprintf('time = %dms', t), 'FontSize', 12);

    subplot(6,3,[10, 11, 13, 14, 16, 17]), hold off;
    plot(lagSpace, [s{2}.p.He(t, :); s{2}.p.Hi(t, :); ...
                            r{2}.A(t, :); s{2}.q.He(t, :)]);
    ylabel('population firing rate (Hz)')
    xlabel('characteristic lag (ms)')
    legend('decoder excitatory', 'decoder inhibitory', ...
           'peridocity detectors', 'sustainer excitatory')
    xlim([0.5, 30])
    ylim([0, 115])
    text(1, 108, 'b) perfect fifth', 'FontSize', 14)

    subplot(6,3,[3,6]), hold off;
    h = zeros(3, 1);
    h(1) = plot(NaN, NaN, 'o', 'Color', c(1, :)); hold on;
    h(2) = plot(NaN, NaN, '.', 'Color', c(2, :)); 
    h(3) = plot(NaN, NaN, '.', 'Color', c(3, :));

    legend(h, 'before onset', ...
              'before convergence', ...
              'after convergence', 'AutoUpdate','off');
    legend('Location','southeast')

    t1 = 75; t2 = 75+155; t3 = 350+75;

    if t < t1       
        plot(Y{1}(1:t, 2), Y{1}(1:t, 1), 'o', 'Color', c(1, :))
    elseif t < t2
        plot(Y{1}(1:t1, 2), Y{1}(1:t1, 1), 'o', 'Color', c(1, :)), hold on
        plot(Y{1}(t1:t, 2), Y{1}(t1:t, 1), '.', 'Color', c(2, :))
    elseif t < t3
        plot(Y{1}(1:t1, 2), Y{1}(1:t1, 1), 'o', 'Color', c(1, :)), hold on
        plot(Y{1}(t1:t2, 2), Y{1}(t1:t2, 1), '.', 'Color', c(2, :))
        plot(Y{1}(t2:t, 2), Y{1}(t2:t, 1), '.', 'Color', c(3, :))
    else
        plot(Y{1}(1:t1, 2), Y{1}(1:t1, 1) + Y{1}(1:t1, 3), 'o', 'Color', c(1, :)), hold on
        plot(Y{1}(t1:t2, 2), Y{1}(t1:t2, 1) + Y{1}(t1:t2, 3), '.', 'Color', c(2, :))
        plot(Y{1}(t2:t3, 2), Y{1}(t2:t3, 1) + Y{1}(t2:t3, 3), '.', 'Color', c(3, :))
        plot(Y{1}(t3:t, 2), Y{1}(t3:t, 1) + Y{1}(t3:t, 3), '.', 'Color', c(4, :))
    end
    xlabel('decoder variables - PCA 1');
    ylabel('decoder variables - PCA 2');
    xlim([-25, 65]) 
    ylim([-5, 225])
    text(-23, 208, 'c) minor second', 'FontSize', 14)


    t1 = 75; t2 = 75+135; t3 = 350+75;
    subplot(6,3,[9, 12]), hold off;
    if t < t1       
        plot(Y{2}(1:t, 2), Y{2}(1:t, 1), 'o', 'Color', c(1, :))
    elseif t < t2
        plot(Y{2}(1:t1, 2), Y{2}(1:t1, 1), 'o', 'Color', c(1, :)), hold on
        plot(Y{2}(t1:t, 2), Y{2}(t1:t, 1), '.', 'Color', c(2, :))
    elseif t < t3
        plot(Y{2}(1:t1, 2), Y{2}(1:t1, 1), 'o', 'Color', c(1, :)), hold on
        plot(Y{2}(t1:t2, 2), Y{2}(t1:t2, 1), '.', 'Color', c(2, :))
        plot(Y{2}(t2:t, 2), Y{2}(t2:t, 1), '.', 'Color', c(3, :))
    else
        hold off; plot(Y{2}(1:t1, 2), Y{2}(1:t1, 1), 'o', 'Color', c(1, :)), hold on
        plot(Y{2}(t1:t2, 2), Y{2}(t1:t2, 1), '.', 'Color', c(2, :))
        plot(Y{2}(t2:t3, 2), Y{2}(t2:t3, 1), '.', 'Color', c(3, :))
        plot(Y{2}(t3:t, 2), Y{2}(t3:t, 1), '.', 'Color', c(4, :))
    end
    xlabel('decoder variables - PCA 1');
    ylabel('decoder variables - PCA 2');
    xlim([-55, 135]) 
    ylim([-5, 220])
    text(-51, 203, 'd) perfect fifth', 'FontSize', 14)
    
    subplot(6,3,[15, 18]), 
    hold off; plot(timeSpace(1:t), mean(s{1}.p.He(1:t, :), 2), 'LineWidth', 2);
    hold on;  plot(timeSpace(1:t), mean(s{2}.p.He(1:t, :), 2), 'LineWidth', 2);
    ylabel('avg firing rate in decoder excit (Hz)');
    xlabel('time (ms)');
    legend('minor second', 'perfect fifth');
    xlim([0, length(timeSpace)])
    ylim([0, 6.75])
    text(9, 6.3,'e)', 'FontSize', 14)

    F(t) = getframe(fig);

end

vid = VideoWriter('decodingDynamicsDyads.avi');
vid.Quality = 100;
vid.FrameRate = 15;
open(vid);
writeVideo(vid, F);
close(vid);