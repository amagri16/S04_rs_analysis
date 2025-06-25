N = 15; % Show top 15 features
topN = min(N, numel(featuresSorted));
figure('Name','Top Features by Frequency','Color','w');
bar(countsSorted(1:topN), 'FaceColor', [0.2 0.6 0.8]);
hold on;
yyaxis right
plot(1:topN, meanAccSorted(1:topN), 'o-r', 'LineWidth', 2, 'MarkerSize', 6);
ylabel('Mean Accuracy (%)');
yyaxis left
set(gca, 'XTick', 1:topN, 'XTickLabel', featuresSorted(1:topN), 'XTickLabelRotation', 45);
ylabel('Number of Tasks Selected');
xlabel('Feature');
title(sprintf('Top %d Features: Frequency (bar) & Mean Accuracy (red line)', topN));
grid on;
hold off;