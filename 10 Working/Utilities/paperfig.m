function []=paperfig(target,fontsize,linewidth,markersize)
% PAPERFIG  Set fontsizes to suit article figures for 2 column format
%

if strcmp(get(target,'type'),'figure')
   trgts = get(target,'child');
   for lp=1:length(trgts)
      modify(trgts(lp),fontsize,linewidth,markersize)
   end
else
   modify(target,fontsize,linewidth,markersize)
end
   
function []=modify(trgt,fontsize,linewidth,markersize)

set(trgt,'fontsize',fontsize)
set(get(trgt,'title'),'fontsize',fontsize)
set(get(trgt,'zlabel'),'fontsize',fontsize)
set(get(trgt,'xlabel'),'fontsize',fontsize)
set(get(trgt,'ylabel'),'fontsize',fontsize)
ch = get(trgt,'children');
n = length(ch);
for lp=1:length(ch)
   if strcmp(get(ch(lp),'type'),'text')
      set(ch(lp),'fontsize',fontsize)
   elseif strcmp(get(ch(lp),'type'),'line')|strcmp(get(ch(lp),'type'),'patch')
      set(ch(lp),'linewidth',linewidth)
      set(ch(lp),'markersize',markersize)
%       linecolor(lp,:) = get(ch(lp),'color');
%       set(ch(lp),'markeredgecolor',linecolor(lp,:))
%       set(ch(lp),'markerfacecolor',linecolor(lp,:))
   end
end