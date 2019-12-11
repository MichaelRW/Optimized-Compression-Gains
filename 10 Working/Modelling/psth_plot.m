% function psth_plot( psth_struct )
% 
% %% AUDIOGRAMS
% 
% figure;
% plot(psth_struct.audiogram_struct.F/1e3, -psth_struct.audiogram_struct.H,'rx', 'Markersize', 20, 'Linewidth', 4); ylim([-120 10]); set(gca, 'XScale', 'linear')
% 
% hold on;
% 
% line(3.5*ones(1,2),[10 -120],'color','k','Linestyle','--');
% line(4.5*ones(1,2),[10 -120],'color','k','Linestyle','--');
% line(5.5*ones(1,2),[10 -120],'color','k','Linestyle','--');
% line(6.5*ones(1,2),[10 -120],'color','k','Linestyle','--');
%   
% line([0 8000], 10*ones(1,2),'color','k');
% line([0 8000], -10*ones(1,2),'color','k');
% line([0 8000], -30*ones(1,2),'color','k');
% line([0 8000], -50*ones(1,2),'color','k');
% line([0 8000], -70*ones(1,2),'color','k');
% line([0 8000], -90*ones(1,2),'color','k');
% line([0 8000], -110*ones(1,2),'color','k');
% line([0 8000], -120*ones(1,2),'color','k');
%     
% hold off;
% set(gca,'XAxisLocation','top');
% set(gca, 'XTick', [1 2 3 4 5 6 7]);
% set(gca, 'XTickLabel', [125 250 500 1000 2000 4000 8000]);
% set(gca, 'YTick', fliplr([0 -20 -40 -60 -80 -100 -120]));
% set(gca, 'YTickLabel', fliplr([0 20 40 60 80 100 120]));
% xlim([1 7]);
% ylim([-120 10]);
% grid;
% set(gca, 'GridLineStyle', '-');
% xlabel('Frequency (Hz)', 'FontSize', 16);
% ylabel('Hearing Threshold Loss (dB)', 'FontSize', 16);
% title([psth_struct.audiogram_struct.type ' Hearing Loss'], 'FontSize', 16);
% 
% %% PLOTTINGS
% if strcmp(psth_struct.type, 'FINE')
% 
%     %% PSTH PLOT (Need to revise max/min)
%     figure;
% %     imagesc(psth_struct.psth_time, psth_struct.psth_freq, psth_struct.psth, psth_struct.psth_mnmx);
%     imagesc(psth_struct.psth_time, psth_struct.psth_freq, abs(log10(psth_struct.psth)), [0 log10(4318)]);
%     
%     axis xy;
%     set(gca, 'YScale', 'log');
%     set(gca, 'YTick', [250 500 750 1000 1500 2000 3000 4000 6000 8000]);
%     set(gca, 'YTickLabel', [250 500 750 1000 1500 2000 3000 4000 6000 8000]);
%     ylabel('Center Frequency (Hz)')
%     xlabel('Time (s)')
%     title('Neurogram')
%     xlim([0 0.03])
% 
% %     %% PWRR
% % 
% %     figure;
% %     plot(psth_struct.psth_freq,psth_struct.pwrr);
% %     xlabel('Centre Frequencies (Hz)')
% %     title(['Power Ratio'])
% %     ylim([0 1]);
% %     hold on;
% %     led = 'legend(gca';
% %     for i = 1:length(psth_struct.data_struct.approx_formants)
% %         num = psth_struct.data_struct.approx_formants(i);
% %         line([num num], [0 1], 'LineStyle', ':', 'Color', 'k');
% %         led = [led ',''' num2str(num) ' Hz'''];
% %     end
% %     led = [led ');'];
% %     eval(led)
% %     hold off;
% % 
% %     %% PHSR
% %     figure;
% %     % plot(psth_struct.phsr_freq,psth_struct.phsr);
% %     hold on;
% %     mx = max([psth_struct.phsr{:}]);
% %     mn = min([psth_struct.phsr{:}]);
% %     for i = 1:length(psth_struct.data_struct.approx_formants)
% %         num = psth_struct.data_struct.approx_formants(i);
% %         plot(psth_struct.phsr_freq{i},psth_struct.phsr{i})
% %         line([num num], [mn mx], 'LineStyle', ':', 'Color', 'k');
% %     end
% %     ylim([mn mx])
% %     ylabel('Relative Phase Shift (deg)')
% %     xlabel('Centre Frequencies (Hz)')
% %     title(['Phase Response'])
% % 
% %     %% BOXP
% %     h1 = figure;
% %     h3 = fill([0.1 10 10 0.1],[0.1*2^(1/2) 10*2^(1/2) 10/2^(1/2) 0.1/2^(1/2)],0.9*ones(1,3));
% %     set(h3,'edgecolor','none')
% %     set(gca,'xscale','log','yscale','log')
% %     title('Box Plot')
% %     xlabel('BF (kHz)')
% %     ylabel('Frequency (kHz)')
% %     axis([0.1 10 0.095 4.5])
% %     set(gca,'xticklabel','0.1|1.0|10')
% %     set(gca,'yticklabel','0.1|1.0')
% %     hold on
% %     plot([0.5 0.5],[0.095 4.5],'k--')
% %     plot([1.7 1.7],[0.095 4.5],'k--')
% %     plot([2.5 2.5],[0.095 4.5],'k--')
% %     plot([0.1 10],[0.5 0.5],'k--')
% %     plot([0.1 10],[1.7 1.7],'k--')
% %     plot([0.1 10],[2.5 2.5],'k--')
% %     wysiwyg
% % 
% %     figure(h1)
% %     for i = 1:length(psth_struct.psth_freq)
% %     h2 = loglog(psth_struct.psth_freq(i)/1e3,psth_struct.boxp_freq/1e3,'ks','markerfacecolor','k');
% %         for lp=1:length(h2)
% %             sz = psth_struct.boxp(i,lp);
% %             if sz < 1
% %                 set(h2(lp),'marker','none')
% %             else
% %                 set(h2(lp),'markersize',2*sz);
% %             end
% %         end
% %     end
% 
% elseif strcmp(psth_struct.type, 'AVG')
%     
%     %% HIST
%     figure;
%     imagesc(psth_struct.hist_bins, psth_struct.psth_freq, psth_struct.hist, [0 895]);
% %     imagesc(psth_struct.hist_bins, psth_struct.psth_freq, psth_struct.hist, psth_struct.hist_mnmx);
% %     imagesc(psth_struct.hist_bins, psth_struct.psth_freq, log10(psth_struct.hist), [0 log10(1101)]);
%     axis xy
%     axis square;
%     set(gca, 'YScale', 'log');
%     set(gca, 'YTick', [250 500 750 1000 1500 2000 3000 4000 6000 8000]);
%     set(gca, 'YTickLabel', [250 500 750 1000 1500 2000 3000 4000 6000 8000]);
%     title('Histogram Response')
%     ylabel('Center Frequency (Hz)')
%     xlabel('Rate Bins (# of spikes/sec/fiber)')
% end