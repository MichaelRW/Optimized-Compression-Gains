function [ fgs ] = error_plot( error )
% ERROR_PLOT plots the error data in the 'error' data structure.

if strcmp(error.type, 'FINE')

%     fgs.phsr = figure;
%     plot(error.ADJ, error.phsr) % legend 
%     title('Phase Response Error')
%     xlabel('Gain Adjustment')
%     ylabel('Error')
%     
%     led = 'legend(gca';
%     for i = 1:length(error.approx_formants)
%         num = error.approx_formants(i);
%         opt = error.phsr_opti(i);
%         led = [led ',''Form: ' num2str(num) 'Hz, Opt: ' num2str(opt) 'dB'''];
%     end
%     led = [led ');'];
%     eval(led);
%     
%     fgs.pwrr = figure;
%     plot(error.ADJ, error.pwrr)
%     title('Power Ratio Error')
%     xlabel('Gain Adjustment')
%     ylabel('Error')
% 
%     led = 'legend(gca';
%     for i = 1:length(error.approx_formants)
%         num = error.approx_formants(i);
%         opt = error.pwrr_opti(i);
%         led = [led ',''Form: ' num2str(num) 'Hz, Opt: ' num2str(opt) 'dB'''];
%     end
%     led = [led ');'];
%     eval(led);
%     
%     fgs.psth = figure;
%     plot(error.ADJ, error.psth)
%     title('Neurogram Error')
%     title(['Neurogram Error, Opt:' num2str(error.psth_opti) ])
%     xlabel('Gain Adjustment')
%     ylabel('Error')
%     
%     fgs.boxp = figure;
%     plot(error.ADJ, error.boxp)
%     title(['Box Error, Opt:' num2str(error.boxp_opti) ])
%     xlabel('Gain Adjustment')
%     ylabel('Error')
       
%     fgs.opti = figure;
%     hold on;
%     plot(error.SPL_uniq, error.psth_opti,'r');
%     plot(error.SPL_uniq, error.pwrr_opti(1,:),'g');
%     plot(error.SPL_uniq, error.boxp_opti,'b');
% %     plot(error.SPL_uniq, error.phsr_opti);
%     xlim([40 90]);
%     ylim([-40 40]);
%     title('Optimal Gains')
%     ylabel('Gain Adjustment (dB)')
%     xlabel('SPL (dB)')
%     legend(gca, 'Neurogram', 'F1 Pwr Ratio', 'Box Plot');
%     hold off;
    
elseif strcmp(error.type, 'AVG')
    
%     fgs.hist = figure;
%     plot(error.ADJ, error.hist)
%     title('Histogram Error')
%     title(['Histogram Error, Opt:' num2str(error.hist_opti) ])
%     xlabel('Gain Adjustment')
%     ylabel('Error')
    
%     fgs.opti = figure;
%     hold on;
%     plot(error.SPL_uniq, error.hist_opti);
%     xlim([40 90]);
%     ylim([-40 40]);
%     title('Optimal Gains')
%     ylabel('Gain Adjustment (dB)')
%     xlabel('SPL (dB)')
%     legend(gca, 'Histogram');
%     hold off;
    
end


end
