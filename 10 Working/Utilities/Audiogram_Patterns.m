function [FREQS, H_IMPAIRED, H_NORMAL] = Audiogram_Patterns(num, plt)
% [FREQS, H_IMPAIRED, H_NORMAL] = Audiogram_Patterns(num, plt)
% plt = 'yes' or 'no'
% num =
% 1: Mild, Gently Sloping
% 2: Moderate, Flat Loss
% 3: Moderate, Steeply Sloping
% 4: Profound, Gently Sloping

if nargin == 1
    plt = 'no';
end

switch num
    case 1
        % Mild, Gently Sloping
        t = 'Mild, Gently Sloping';
        F = [250 500 1000 2000 3000 4000 6000 8000];
        H_IMP = [20 20 30 40 45 50 50 50];
        H_NOR = [0 0 0 0 0 0 0 0];
    case 2
        % Moderate, Flat Loss
        t = 'Moderate, Flat Loss';
        F = [250 500 1000 2000 4000 6000 8000];
        H_IMP = [40 40 40 40 40 40 40];
        H_NOR = [0 0 0 0 0 0 0];
    case 3
        % Moderate, Steeply Sloping
        t = 'Moderate, Steeply Sloping';
        F = [250 500 1000 1500 2000 3000 4000 6000 8000];
        H_IMP = [25 30 55 65 80 85 90 90 90];
        H_NOR = [0 0 0 0 0 0 0 0 0];
    case 4
        % Profound, Gently Sloping
        t = 'Profound, Gently Sloping';
        F = [250 500 1000 1500 2000 3000 4000 6000 8000];
        H_IMP = [80 85 95 100 110 110 110 110 110];
        H_NOR = [0 0 0 0 0 0 0 0 0];
    otherwise
        disp('Not a Selection')
        F = [];
        H_IMP = [];
        H_NOR = [];
end

% FREQS = [250 500 750 1000 1500 2000 3000 4000 6000 8000];
% H_IMPAIRED = interp1(F,H_IMP,FREQS,'linear');
% H_NORMAL = interp1(F,H_NOR,FREQS,'linear');
% 
% if strncmpi(plt, 'yes',2)
%     plot(F, -H_IMP,'r+'); ylim([-120 0]); set(gca, 'XScale', 'log')
%     hold on;
%     plot(FREQS, -H_IMPAIRED,'bO'); ylim([-120 0]); set(gca, 'XScale', 'log')
%     hold off;
% end


FREQS = [250 500 750 1000 1500 2000 3000 4000 6000 8000];
H_IMPAIRED = interp1(F,H_IMP,FREQS,'linear');
H_NORMAL = interp1(F,H_NOR,FREQS,'linear');
F = floor(log2(F/125)/0.5)/2 + 1;
if nargout == 0
%     figure;
%     plot(F, -H_IMP,'rx', 'Markersize', 20, 'Linewidth', 4); ylim([-120 10]); set(gca, 'XScale', 'linear')
%     hold on;
% %     plot(FREQS, -H_IMPAIRED,'bO'); ylim([-120 5]); set(gca, 'XScale', 'linear')
%     
%     
%     line(3.5*ones(1,2),[10 -120],'color','k','Linestyle','--')
%     line(4.5*ones(1,2),[10 -120],'color','k','Linestyle','--')
%     line(5.5*ones(1,2),[10 -120],'color','k','Linestyle','--')
%     line(6.5*ones(1,2),[10 -120],'color','k','Linestyle','--')
%   
%     line([0 8000], 10*ones(1,2),'color','k')
%     line([0 8000], -10*ones(1,2),'color','k')
%     line([0 8000], -30*ones(1,2),'color','k')
%     line([0 8000], -50*ones(1,2),'color','k')
%     line([0 8000], -70*ones(1,2),'color','k')
%     line([0 8000], -90*ones(1,2),'color','k')
%     line([0 8000], -110*ones(1,2),'color','k')
%     line([0 8000], -120*ones(1,2),'color','k')   
%     
%     hold off;
%     set(gca,'XAxisLocation','top')
%     set(gca, 'XTick', [1 2 3 4 5 6 7]);
%     set(gca, 'XTickLabel', [125 250 500 1000 2000 4000 8000]);
%     set(gca, 'YTick', fliplr([0 -20 -40 -60 -80 -100 -120]));
%     set(gca, 'YTickLabel', fliplr([0 20 40 60 80 100 120]));
%     xlim([1 7])
%     ylim([-120 10])
%     grid;
%     set(gca, 'GridLineStyle', '-')
%     xlabel('Frequency (Hz)', 'FontSize', 16);
%     ylabel('Hearing Threshold Loss (dB)', 'FontSize', 16);
%     title([t ' Hearing Loss'], 'FontSize', 16)
end