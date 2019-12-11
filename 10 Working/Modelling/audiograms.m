function audiogram = audiograms(num, CF_Points, varargin)
%AUDIOGRAMS returns a sample audiogram profile.
% 
% Audiogram_Strucutre = audiograms(LookUpNumber, CF_Points) uses a 15 row
% look-up take to return a sample audiogram. LookUpNumber number 1 contains
% a normal audiogram profile with an overall zero dB loss. Each audiogram
% profile contains at least 10 logarithmically spaced hearing losses. The
% parameter CF_Points specifies the number of threshold points along the
% audiogram profile, with 10 being the least.
%  
% REF: Susan Scollie, Richard Seewald, Leonard Cornelisse, Sheila Moodie,
% Marlene Bagatto, Diana Laurnagaray, Steve Beaulac and John Pumford. The
% Desired Sensation Level Multistage Input/Output Algorithm. Trends Amplif.
% 2005;9(4):159-97.

AG = [[ 0  0  0  0  0  0  0  0  0  0];...
       [20 20 25 30 40 45 50 50 50 50];...
       [25 30 45 55 65 80 85 90 90 90];...  
       [30 30 30 30 30 30 30 30 30 30];...
       %[30 30 35 35 40 40 45 45 50 50];...
       [30 30 35 40 45 50 60 70 75 75];...
       [30 35 40 45 50 60 80 90 90 90];...
       [40 40 40 40 40 40 40 40 40 40];...
       [50 50 50 50 50 50 50 50 50 50];...
       [50 50 55 55 60 60 65 65 70 70];...
       [50 50 50 55 65 70 80 90 95 95];...
       [70 70 70 70 70 70 70 70 70 70];...
       [70 70 75 75 80 80 85 85 90 90];...
       [90 90 90 90 90 90 90 90 90 90];...
       [80 85 90 95 100 110 110 110 110 110]];

type = {'Normal', 'Mild', 'Moderate', 'Mild', 'Mild', 'Moderate', 'Moderate', 'Moderate', 'Moderate', 'Moderate', 'Moderate', 'Moderate', 'Severe', 'Severe', 'Profound',};

%type = {'Normal','Mild'};
if num < 1
    num = 1;
end

if num > 15
    num = 15;
end

%%
t = type{num};
F = [250 500 750 1000 1500 2000 3000 4000 6000 8000];

if exist('CF_Points') && ( CF_Points > 10 )
    FREQS = logspace( log10( F(1) ), log10( F(end) ), CF_Points );
    H = interp1(F,AG(num,:),FREQS,'pchip'); % 'linear' causes problems!
else
    FREQS = F;
    H = AG(num,:);
end

if ~isempty(varargin)
    IOHC_loss = varargin{1};
else
    IOHC_loss = 'Mixed';
end
    
if strncmpi(IOHC_loss, 'OHCL', 4)
    [Cohc,Cihc]=fitaudiogram2(FREQS,H,1);
    disp('Impairment due to Outer Hair Cell Loss');
elseif strncmpi(IOHC_loss, 'IHCL', 4)
    [Cohc,Cihc]=fitaudiogram2(FREQS,H,1);
    disp('Impairment due to Inner Hair Cell Loss');
else
    IOHC_loss = 'Mixed';
    [Cohc,Cihc]=fitaudiogram2(FREQS,H,1);
    disp('Mixed Inner & Outer Hair Cell Loss');
end

%%
audiogram.type = t;
audiogram.F = round(FREQS);
audiogram.H = H;
audiogram.Cohc = Cohc;
audiogram.Cihc = Cihc;
audiogram.orig_F = F;
audiogram.orig_H = AG(num,:);
audiogram.IOHC_loss = IOHC_loss;

%%
if nargout == 0
    F = floor(log2(F/125)/0.5)/2 + 1;
%     figure;
%     plot(F, -AG(num,:),'rx', 'Markersize', 20, 'Linewidth', 4); ylim([-120 10]); set(gca, 'XScale', 'linear')
%     hold on;
% 
%     line(3.5*ones(1,2),[10 -120],'color','k','Linestyle','--');
%     line(4.5*ones(1,2),[10 -120],'color','k','Linestyle','--');
%     line(5.5*ones(1,2),[10 -120],'color','k','Linestyle','--');
%     line(6.5*ones(1,2),[10 -120],'color','k','Linestyle','--');
%   
%     line([0 8000], 10*ones(1,2),'color','k');
%     line([0 8000], -10*ones(1,2),'color','k');
%     line([0 8000], -30*ones(1,2),'color','k');
%     line([0 8000], -50*ones(1,2),'color','k');
%     line([0 8000], -70*ones(1,2),'color','k');
%     line([0 8000], -90*ones(1,2),'color','k');
%     line([0 8000], -110*ones(1,2),'color','k');
%     line([0 8000], -120*ones(1,2),'color','k');
%     
%     hold off;
%     set(gca,'XAxisLocation','top');
%     set(gca, 'XTick', [1 2 3 4 5 6 7]);
%     set(gca, 'XTickLabel', [125 250 500 1000 2000 4000 8000]);
%     set(gca, 'YTick', fliplr([0 -20 -40 -60 -80 -100 -120]));
%     set(gca, 'YTickLabel', fliplr([0 20 40 60 80 100 120]));
%     xlim([1 7]);
%     ylim([-120 10]);
%     grid;
%     set(gca, 'GridLineStyle', '-');
%     xlabel('Frequency (Hz)', 'FontSize', 16);
%     ylabel('Hearing Threshold Loss (dB)', 'FontSize', 16);
%     title([t ' Hearing Loss'], 'FontSize', 16);
end