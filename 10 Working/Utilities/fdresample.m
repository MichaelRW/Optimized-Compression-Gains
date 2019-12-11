function [ xr ] = fdresample( x,fs_res,fs_orig )

t = 0:1/fs_orig:length(x)/fs_orig; t(end) = [];
tr = 0:1/fs_res:length(x)/fs_orig; tr(end) = [];
xr = interp1(t,x,tr,'spline');