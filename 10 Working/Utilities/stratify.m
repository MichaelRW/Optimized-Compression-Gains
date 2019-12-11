function [ array ] = stratify( data, seg, type ,ovl );
% window = 256;
% shift = window*0.5;
% data = stratify(data, window, 'shift', shift);

if strcmp(type, 'shift')
    lip = ovl;
elseif strcmp(type, 'overlap')
    lip = (seg - ovl);
else
    disp('type help');
    array = [];
    return
end

len = length(data);
remainder = rem((len - seg),lip);

if remainder == 0;
    n = (len - seg)/lip;
else
    data = [data zeros(1,(lip - remainder))];
    len = length(data);
    n = (len - seg)/lip;
end

seg = 1:seg;
array = [];

for i = 0:n
    array(i+1,:) = data(seg + i*lip);
end