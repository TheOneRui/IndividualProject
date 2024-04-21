tiledlayout(2,1)

nexttile
hold on
for i = 1 : length(results_half) %Plot all the results in one Figure
    plot(results_half{4,i});
    plot(results_half{5,i});
end
yline(0.99, 'k', 'DisplayName', 'lower Statutory Limit');
yline(1.01, 'k', 'DisplayName', 'upper Statutory Limit');
plot_title = sprintf('Limit Violating Response Curves - Halved Injection');
title(plot_title);
xlabel('Time (s) since LoG Event')
ylabel('Frequency (pu)')
legend('show')
hold off

nexttile
hold on
for i = 1 : length(results) %Plot all the results in one Figure
    plot(results{4,i});
    plot(results{5,i});
end
yline(0.99, 'k', 'DisplayName', 'lower Statutory Limit');
yline(1.01, 'k', 'DisplayName', 'upper Statutory Limit');
plot_title = sprintf('Limit Violating Response Curves - Ideal Injection');
title(plot_title);
xlabel('Time (s) since LoG Event')
ylabel('Frequency (pu)')
legend('show')
hold off