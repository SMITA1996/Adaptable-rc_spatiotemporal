function theta = find_tipping_TGAM_exact(y, plotFlag)
% EXACT TGAM threshold detection matched to Python UnivariateSpline.
%
%   theta = find_tipping_TGAM_exact(y);
%   theta = find_tipping_TGAM_exact(y, 0);  % no plot
%
%   Uses spaps so smoothing tolerance matches Python's UnivariateSpline(..., s=len)

if nargin < 2
    plotFlag = 1;
end

y = y(:);
t = (1:length(y))';

% Candidate thresholds 
candidates = (round(0.1*length(y)) : 10 : round(0.9*length(y)))';

bestAIC = inf;
theta = NaN;

for th = candidates'

    idx1 = (t <= th);
    idx2 = (t > th);

    [sp1, ~] = spaps(t(idx1), y(idx1), length(y(idx1)));
    [sp2, ~] = spaps(t(idx2), y(idx2), length(y(idx2)));

    y_hat = zeros(size(y));
    y_hat(idx1) = fnval(sp1, t(idx1));
    y_hat(idx2) = fnval(sp2, t(idx2));

    % AIC 
    rss = sum((y - y_hat).^2);
    k = 6;  % constant DF 
    AIC = length(y)*log(rss/length(y)) + 2*k;

    if AIC < bestAIC
        bestAIC = AIC;
        theta = th;
        bestf1 = sp1;
        bestf2 = sp2;
    end
end

fprintf('\nEstimated tipping point = %d\n', theta);

%% Plot
if plotFlag
    figure; hold on;
    plot(t, y, 'k', 'LineWidth', 1.1);
    plot(t(t <= theta), fnval(bestf1, t(t <= theta)), 'r', 'LineWidth', 2);
    plot(t(t > theta), fnval(bestf2, t(t > theta)), 'b', 'LineWidth', 2);
    xline(theta, '--m', 'LineWidth', 2);

    legend('Data', 'Fit Left of θ', 'Fit Right of θ', 'Detected θ');
    xlabel('Time Index');
    ylabel('Value');
    title(sprintf('TGAM Tipping Point : t = %d', theta));
    grid on;
end

end
