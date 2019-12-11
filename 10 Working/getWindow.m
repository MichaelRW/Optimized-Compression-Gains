
function [ window ] = getWindow( start, stop, FS, psth_time, psth_freq )

window_start = start/FS+(5e-3);
window_end = stop/FS+(5e-3);

start_index = find(abs((psth_time(1,:)-window_start))<1e-4,1);
finish_index = find(abs((psth_time(1,:)-window_end))<1e-4,1);


ramp_indices = round(10e-3*FS);
ramp = zeros(1,ramp_indices);

for ramp_count = 1:ramp_indices-1
    ramp(1,ramp_count) = (ramp_indices-ramp_count)/ramp_indices;
end

window = [zeros(1,start_index-1) ones(1,finish_index-start_index+1) ramp];

if length(window) > length(psth_time)
    window = window(:,1:length(psth_time));
else
    window = [window zeros(1,length(psth_time)-length(window))];
end

window = repmat(window,length(psth_freq),1);


