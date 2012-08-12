function [result] = STA(ts, beh, lags ,toplot, subplot_r, subplot_c)
if toplot == 1
    f = figure;
    set(f, 'Color', 'w');
end
ts = ceil(ts);
while ts(1) <= -lags(1)
    ts = ts(2:end);
end
while ts(end) >= size(beh,2) - lags(end)
    ts = ts(1:end-1);
end
averages = [];
upper_limits = [];
down_limits = [];
for var = 1 : size(beh,1)
dat = beh(var, :);
all = zeros(length(ts), length(lags));
for i = 1:length(ts)
    t = ts(i);
    disp(t);
    all(i,:) = dat(t+lags);
end
baseline = mean(all,1) - mean(mean(all(:, 1:-lags(1)), 1));
upper = baseline + 1.96*(std(all,[],1)./sqrt(length(ts)));
down = baseline - 1.96*(std(all,[],1)./sqrt(length(ts)));

averages = [averages;baseline];
upper_limits = [upper_limits; upper];
down_limits = [down_limits; down];
if toplot == 1
    subplot(subplot_r, subplot_c, var);
    plot(lags/20, baseline, 'r', lags/20, upper, 'b', lags/20, down, 'b');
end
end
result = struct('av', averages, 'up', upper_limits, 'down', down_limits);
end
