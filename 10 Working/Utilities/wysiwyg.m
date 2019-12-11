function wysiwyg(figh)
% Based on a function from the MATDRAW toolbox, Keith Rogers 3/95
% Modified by Ian C. Bruce (ibruce@ieee.org) June 2007

if nargin==0
    figh=gcf;
end
units=get(figh,'units');
set(figh,'units',get(figh,'PaperUnits'));
set(figh,'Position',get(figh,'PaperPosition'));
set(figh,'Units',units)
set(figh,'Position',get(figh,'Position')-[0 50 0 0]) % use this line if figure window ends up partly off the screen