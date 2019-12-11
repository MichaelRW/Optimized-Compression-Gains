function [ fgs ] = compare_plot( psth1, psth2, filename )
% COMPARE_PLOT Plots together the Power Ratio, Phase Response, Box Plot,
% and Histogram Analysis of two neurograms for easy comparison.
% 
% fgs = compare_plot( Neurogram_Struct1, Neurogram_Struct2, filename )
% The 'normal' neurogram should be the first and the 'impaired' should be
% second in the input argument list. Figures are saved in the current
% directory as PWRR, PHSR, BOXP, and HIST with 'filename' appended.


if strcmp(psth1.type, 'FINE')
%% PWRR

    PWRR = figure;
    plot(psth1.psth_freq,psth1.pwrr,'LineWidth', 2);
    ylabel('Synchrony')
    xlabel('Centre Frequencies (Hz)')
    title(['Power Ratio'])
    ylim([0 1]);
    hold on;
    led = 'legend(gca';
    for i = 1:length(psth1.data_struct.approx_formants)
        num = psth1.data_struct.approx_formants(i);
        line([num num], [0 1], 'LineStyle', ':', 'Color', 'k');
        led = [led ',''' num2str(num) ' Hz'''];
    end
    led = [led ');'];
    eval(led)
    plot(psth2.psth_freq,psth2.pwrr, '-.', 'LineWidth', 2);
    hold off;
    wysiwyg;
    if exist('filename')
        hgsave(PWRR,['PWRR_' filename]);
    end
    
%% PHSR

    PHSR = figure;
    hold on;
    
    [ phsr_freq1, phsr_freq2, phsr1, phsr2 ] = same_phsr_range( psth1, psth2 );
    
    mx = max([phsr1{:} phsr2{:}]);
    mn = min([phsr1{:} phsr2{:}]);
    
    for i = 1:length(psth1.data_struct.approx_formants)
        num = psth1.data_struct.approx_formants(i);
        plot(phsr_freq1{i},phsr1{i},'b')
        plot(phsr_freq2{i},phsr2{i},'r')
        line([num num], [mn mx], 'LineStyle', ':', 'Color', 'k');
    end
    ylim([mn mx]);
    xlim([0 8000]);
    ylabel('Relative Phase Shift (deg)');
    xlabel('Centre Frequencies (Hz)');
    title(['Phase Response']);
    wysiwyg;
    if exist('filename')
        hgsave(PHSR, ['PHSR_' filename]);
    end
    
%% BOXP
    BOXP = figure;
    h2 = fill([0.1 10 10 0.1],[0.1*2^(1/2) 10*2^(1/2) 10/2^(1/2) 0.1/2^(1/2)],0.9*ones(1,3));
    set(h2,'edgecolor','none')
    set(gca,'xscale','log','yscale','log')
    title('Box Plot')
    xlabel('BF (kHz)')
    ylabel('Frequency (kHz)')
    axis([0.1 10 0.095 4.5])
    set(gca,'xticklabel','0.1|1.0|10')
    set(gca,'yticklabel','0.1|1.0')
    hold on
    plot([0.5 0.5],[0.095 4.5],'k--')
    plot([1.7 1.7],[0.095 4.5],'k--')
    plot([2.5 2.5],[0.095 4.5],'k--')
    plot([0.1 10],[0.5 0.5],'k--')
    plot([0.1 10],[1.7 1.7],'k--')
    plot([0.1 10],[2.5 2.5],'k--')

    small_box = ((psth1.boxp - psth2.boxp) <= 0) .* psth1.boxp;
    zro = ((psth1.boxp - psth2.boxp) == 0);
    small_box(zro) = small_box(zro)*(-1);
    figure(BOXP)
    for i = 1:length(psth1.psth_freq)
    h3 = loglog(psth1.psth_freq(i)/1e3,psth1.boxp_freq/1e3,'ks','markerfacecolor','k','markeredgecolor',[0 0 0]);
        for lp=1:length(h3)
            sz = psth1.boxp(i,lp);
            if sz < 1
                set(h3(lp),'marker','none')
            else
                set(h3(lp),'markersize',2*sz);
            end
        end

    h4 = loglog(psth2.psth_freq(i)/1e3,psth2.boxp_freq/1e3,'ks','markerfacecolor','r','markeredgecolor',[1 0 0]);
        for lp=1:length(h4)
            sz = psth2.boxp(i,lp);
            if sz < 1
                set(h4(lp),'marker','none')
            else
                set(h4(lp),'markersize',2*sz);
            end
        end

    h5 = loglog(psth1.psth_freq(i)/1e3,psth1.boxp_freq/1e3,'ks','markerfacecolor','k','markeredgecolor',[0 0 0]);
        for lp=1:length(h5)
            sz = small_box(i,lp);
            if sz < 0
                set(h5(lp),'markersize',-2*sz, 'markerfacecolor', 'g', 'markeredgecolor', 'g' );
            elseif abs(sz) < 1
                set(h5(lp),'marker','none')
            else
                set(h5(lp),'markersize',2*sz);
            end
        end
    end
    wysiwyg;
    if exist('filename')
        hgsave(BOXP,['BOXP_' filename]);
    end
    
    fgs.PWRR = PWRR;
    fgs.PHSR = PHSR;
    fgs.BOXP = BOXP;
    
elseif strcmp(psth1.type, 'AVG')
%% HIST
    HIST = figure;
    title('Histogram Response')
    mx = max([psth1.hist_mnmx psth2.hist_mnmx]);
    mn = min([psth1.hist_mnmx psth2.hist_mnmx]);
    
    subplot(211)
    imagesc(psth1.hist_bins, psth1.psth_freq, psth1.hist, [mn mx]); axis xy
    set(gca, 'YScale', 'log');
    set(gca, 'YTick', [250 500 750 1000 1500 2000 3000 4000 6000 8000]);
    set(gca, 'YTickLabel', [250 500 750 1000 1500 2000 3000 4000 6000 8000]);
    ylabel('Center Frequency (Hz)');
    xlabel('Rate Bins (spikes/sec/fiber)');
    axis square
    
    subplot(212)
    imagesc(psth2.hist_bins, psth2.psth_freq, psth2.hist, [mn mx]); axis xy
    set(gca, 'YScale', 'log');
    set(gca, 'YTick', [250 500 750 1000 1500 2000 3000 4000 6000 8000]);
    set(gca, 'YTickLabel', [250 500 750 1000 1500 2000 3000 4000 6000 8000]);
    ylabel('Center Frequency (Hz)');
    xlabel('Rate Bins (spikes/sec/fiber)');
    axis square
    wysiwyg;
    % Scale bins by 40 and relable to 0 -> 200
    if exist('filename')
        hgsave(HIST,['HIST_' filename]);
    end
    
    fgs.HIST = HIST;
    
end

end
