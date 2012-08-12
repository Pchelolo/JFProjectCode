function [r1, r2] = plotTimeXcorr( x, y, delta )
    timestamps = 0 : delta : max(x)+delta;
    disp(max(timestamps));
    %disp(max()
    r1 = zeros(1, length(timestamps));
    k = 1;
    for idx = 1 : length(x)
       while timestamps(k) < x(idx)
          k = k + 1;
       end
       r1(k) = r1(k) + 1;
    end
    timestamps = 0 : delta : max(y) + delta;
    r2 = zeros(1, length(timestamps));
    k = 1;
    for idx = 1 : length(y)
        while timestamps(k) < y(idx)
           k = k + 1;
        end
           r2(k) = r2(k) + 1;
    end
    figure;
    subplot(2,1,1);
    plot(-ceil(1000 / delta):ceil(1000 / delta), xcorr(r1, r2, ceil(1000 / delta)));
    subplot(2,1,2);
    hold on;
    plot(1:length(r1), r1*2);
    plot(1:length(r2), r2, 'r');
    hold off;
end

