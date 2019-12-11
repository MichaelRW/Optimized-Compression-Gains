function imagefd( time, freq, data, clim )
imagesc(time, freq, data, clim);
axis xy;
set(gca, 'YScale', 'log');
set(gca, 'YTick', [250 500 750 1000 1500 2000 3000 4000 6000 8000]);
set(gca, 'YTickLabel', [250 500 750 1000 1500 2000 3000 4000 6000 8000]);

