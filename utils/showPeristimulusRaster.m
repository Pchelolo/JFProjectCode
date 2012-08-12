function showPeristimulusRaster(ts, r)
    lag = -20:20;
    av = [];
    ts = ceil(ts);
    while ts(1) <= -lag(1)
        ts = ts(2:end);
    end
    while ts(end) >= length(r) - lag(end)
        ts = ts(1:end-1);
    end  
    for t = 1:length(ts)
        av = [av; r(ts(t)+lag)];
    end
    f = figure;
    set(f, 'Color', 'w');
    %[~, perm] = sort(mean(av(:, 1:-lag(1)), 2));
    %size(perm)
    image(ceil(av * 60 / max(max(av))));
    colormap(hot);
end

