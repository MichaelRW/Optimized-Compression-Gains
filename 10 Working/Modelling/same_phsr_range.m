function [ x1o, x2o, y1o, y2o ] = same_phsr_range( psth_struct1, psth_struct2 )
% Two problems: if length of either x is <= 1, or if mutually exclusive.
% Check is 'FINE' and that same dimentions!
% Check is both length(psth_struct1&2.phsr_freq) are the same!

if (length(psth_struct1.phsr_freq) ~= length(psth_struct2.phsr_freq)) && ~strcmp(psth_struct.type, 'FINE')
    disp('Stuctures cannot be compared')
end

for i = 1:length(psth_struct1.phsr_freq)

    x1 = psth_struct1.phsr_freq{i};
    x2 = psth_struct2.phsr_freq{i};
    y1 = psth_struct1.phsr{i};
    y2 = psth_struct2.phsr{i};

    if ((length(x1) <=1) || (length(x2) <=1))
        x1o{i} = 0;
        x2o{i} = 0;
        y1o{i} = 0;
        y2o{i} = 0;
    else

        if x1(1) <= x2(1)
            xmin = x2(1);
        else
            xmin = x1(1);
        end

        if x1(end) <= x2(end)
            xmax = x1(end);
        else
            xmax = x2(end);
        end

        if xmax <= xmin
            x1o{i} = 0;
            x2o{i} = 0;
            y1o{i} = 0;
            y2o{i} = 0;
        else
            x1o{i} = x1(find(x1 == xmin):find(x1 == xmax));
            y1o{i} = y1(find(x1 == xmin):find(x1 == xmax));
            x2o{i} = x2(find(x2 == xmin):find(x2 == xmax));
            y2o{i} = y2(find(x2 == xmin):find(x2 == xmax));            
        end
    end
end