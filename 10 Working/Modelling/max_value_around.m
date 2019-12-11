function [ freq, val ] = max_value_around( data_struct, fi)


f0 = data_struct.fundamental;
tol = 0.05*f0;
fi = round(fi / f0) * f0;
f = data_struct.F;
X = data_struct.DATA;

for i = 1:length(fi)
    
    mask = (f < fi(i) + tol) & (f > fi(i) - tol);
    [val(i), cel] = max(X.*mask);
    freq(i) = f(cel);
    
end