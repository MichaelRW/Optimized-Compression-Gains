## Copyright (C) 2018 Helen
## 
## This program is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see
## <https://www.gnu.org/licenses/>.

## -*- texinfo -*- 
## @deftypefn {} {@var{retval} =} getWindow (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Helen <Helen@DESKTOP-7BI49KH>
## Created: 2018-06-19

function window = getWindow(start,stop, FS, psth_time,psth_freq)
        window_start = start/FS+(5e-3); % start  
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
endfunction






