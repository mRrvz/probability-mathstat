X = [-10.82,-9.27,-9.65,-9.36,-9.27,-11.25,-9.89,-9.26,-11.15,-8.90,-11.02,-8.28,-9.18,-10.16,-10.59,-10.82,-9.05,-9.47,-10.98,-11.50,-8.64,-10.86,-10.76,-11.49,-11.09,-9.33,-9.32,-9.66,-8.79,-10.54,-9.12,-10.40,-8.59,-10.22,-9.06,-10.59,-10.60,-10.25,-9.35,-11.44,-9.81,-9.32,-9.95,-9.33,-10.64,-9.45,-10.99,-10.15,-10.39,-10.36,-10.49,-11.67,-10.00,-10.87,-11.11,-9.68,-10.77,-9.13,-8.62,-10.33,-11.36,-10.24,-9.41,-11.05,-10.15,-9.35,-11.45,-9.87,-10.41,-10.11,-10.84,-11.48,-7.77,-10.79,-9.88,-10.70,-9.07,-9.47,-10.15,-9.93,-11.52,-9.04,-10.93,-10.13,-9.56,-11.39,-9.79,-9.19,-11.09,-9.86,-10.67,-10.26,-9.07,-10.53,-11.24,-10.16,-11.33,-8.76,-8.88,-10.53,-10.12,-8.98,-9.84,-9.90,-10.13,-9.32,-9.31,-9.99,-8.55,-11.64,-11.32,-10.51,-11.71,-10.50,-10.50,-12.20,-11.68,-10.45,-7.88,-10.84]
gamma = 0.9;

% 1-2
[muhat, muci] = my_normfit_mu(X, 1 - gamma);
[s2hat, s2ci] = my_normfit_s2(X, 1 - gamma);

% 3
process_mu(X, gamma, muhat);
process_s2(X, gamma, s2hat);


function [muhat, muci] = normfit_mu(X, alpha)
    [muhat, ~, muci, ~] = normfit(X, alpha);
end

function [s2hat, s2ci] = normfit_s2(X, alpha)
    [~, sigmahat, ~, sigmaci] = normfit(X, alpha);
    s2hat = sigmahat ^ 2;
    s2ci = sigmaci .^ 2;
end

function [muhat, muci] = my_normfit_mu(X, alpha)
    muhat = mean(X);
    s = std(X);
    gamma = 1 - alpha;
    n = length(X);
    mu_bottom = muhat + s * tinv((1 - gamma) / 2, n - 1) / sqrt(n);
    mu_top = muhat + s * tinv((1 + gamma) / 2, n - 1) / sqrt(n);
    muci = [mu_bottom, mu_top];
end

function [s2hat, s2ci] = my_normfit_s2(X, alpha)
    s2hat = var(X);
    gamma = 1 - alpha;
    n = length(X);
    s2_top = (n - 1) * s2hat / chi2inv((1 - gamma) / 2, n - 1);
    s2_bottom = (n - 1) * s2hat / chi2inv((1 + gamma) / 2, n - 1);
    s2ci = [s2_bottom, s2_top];
end

function process_parameter(X, gamma, est, fit, line_legend, est_legend, top_legend, bottom_legend)
    N = length(X);
    figure;
    hold on;
    grid on;
    plot([1, N], [est, est]);
    ests = [];
    cis_bottom = [];
    cis_top = [];
    for n = 1:N
        [est, cis] = fit(X(1:n), 1 - gamma);
        ests = [ests, est];
        cis_bottom = [cis_bottom, cis(1)];
        cis_top = [cis_top, cis(2)];
    end
    plot(1:N, ests);
    plot(1:N, cis_bottom);
    plot(1:N, cis_top);
    l = legend(line_legend, est_legend, top_legend, bottom_legend);
    set(l, 'Interpreter', 'latex', 'fontsize', 18);
    hold off;
end

function process_mu(X, gamma, muhat)
    process_parameter(X, gamma, muhat, @my_normfit_mu, '$\hat\mu(\vec x_N)$', '$\hat\mu(\vec x_n)$', '$\underline\mu(\vec x_n)$', '$\overline\mu(\vec x_n)$');
end

function process_s2(X, gamma, S2)
    process_parameter(X, gamma, S2, @my_normfit_s2, '$\hat\sigma^2(\vec x_N)$', '$\hat\sigma^2(\vec x_n)$', '$\underline\sigma^2(\vec x_n)$', '$\overline\sigma^2(\vec x_n)$');
end

