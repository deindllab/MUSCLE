%% Initialize
% select any cell array
% filePrefix = 'Mismatch_ident_vs_nonident_';
filePrefix = 'Example_seq_';
% filePrefix = 'Outlier_seq_after_3G_GATC_';
cell = mismatch;

sel_FRET = sel_FRET_B;
% cell = mismatch;
clear data
n_points = length(cell);
data_FRET{n_points} = [];
data_K{n_points} = [];
axesPosition = [0.2 0.2 0.6 0.6];

for i = 1:n_points
    temp = cell{i};
    data_FRET{i} = [];
    data_K{i} = [];
    for j = 1:size(temp,1)
        x_plot = temp(j,1);
        y_plot = temp(j,2);
        t = sel_FRET(x_plot,y_plot);
        t = t{1};
        t = t(t>-inf);
        data_FRET{i} = [data_FRET{i} t];
        data_K{i} = [data_K{i} htmap_K(x_plot,y_plot)];

    end
end
%% SEM

data_mean = zeros(n_points,1);

errors = zeros(n_points,1);
for i = 1:n_points
%     i
    temp = data_dG{i};
    errors(i) = std(temp)/sqrt(length(temp));
    mean_FRET = mean(temp); 
    data_mean(i) = mean_FRET;
end


%% Combining data for dG vs T

dG_vs_T_10C = data_mean;
dG_vs_T_10C_SEM = errors;
% dG_vs_T_15C = data_mean;
% dG_vs_T_15C_SEM = errors;
% dG_vs_T_20C = data_mean;
% dG_vs_T_20C_SEM = errors;

%% Plot
% errors = [(data_mean - lowerCIs)'; (upperCIs - data_mean)'];
fig = figure('Name','LinePlot','Units', 'pixels', 'Position', [100 300 240 180]);
axes('Position', axesPosition);
x_plot = 0:n_points-1;
% axes('OuterPosition', [0 0 1 1]); % Adjust the position of the axes
% plot(x_plot, data_mean, 'k-', 'LineWidth', 2,'Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
errorbar(x_plot, data_mean, errors, 'k-', 'LineWidth', 2,'Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');

% Adjusting plot appearance
set(gca, 'FontSize', 14, 'FontName', 'Arial', 'TickDir', 'out', 'LineWidth', 2);
set(gca, 'XColor', 'k', 'YColor', 'k');
set(gcf, 'color', 'w');
set(gca,'TickLength',[0.02, 0.01]);
box off;

% grid minor;
% xlabel('# of consecutive mismatches', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial', 'Color', 'k');
% ylabel('dG (kcal/mol)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial', 'Color', 'k');
ylabel('K', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial', 'Color', 'k');
figTitle = 'dG vs base pairs';
% figTitle = 'K vs base pairs';
title(strrep(filePrefix, '_', ' '), 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');

% Adjusting axis limits and tick marks
xlim([-0.5, 3.5]);
% ylim([0.3, 0.8]);
xticks(min(x_plot):1:max(x_plot));
% yticks(0.4:0.05:0.55);

if ~exist('saveFolder','var')
    saveFolder = uigetdir('','Select a folder to save figures in');
end
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
% Define base file name
baseFileName = [filePrefix, cleanFileName(figTitle),'_'];
% savefig([path1 '\average_plot']);
% Add timestamp to base file name
figFileName = fullfile(saveFolder, [baseFileName, timestamp, '.fig']);
pdfFileName = fullfile(saveFolder, [baseFileName, timestamp, '.pdf']);

% Save figure as .fig
saveas(fig, figFileName);

% Save figure as .pdf
exportgraphics(fig,pdfFileName,'ContentType','vector');
disp('Done!');

%% Plot histograms
fitting_flag = false;
if ~exist('saveFolder','var')
    saveFolder = uigetdir('','Select a folder to save figures in');
end
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
for i = 1:n_points
    figTitle = ['Hist_',num2str(i)];
    fig = figure('Name',figTitle,'Units', 'pixels', 'Position', [100 300 240 180]);
    axes('Position', axesPosition);
    histogram(data_FRET{i}, -0.2:0.05:1.2 , 'FaceColor', [0 174 239] / 255, 'EdgeColor', 'k', 'LineWidth', 1);
    % Adjusting plot appearance
    set(gca, 'FontSize', 14, 'FontName', 'Arial', 'TickDir', 'out', 'LineWidth', 2);
    set(gca, 'XColor', 'k', 'YColor', 'k');
    set(gcf, 'color', 'w');
    set(gca,'TickLength',[0.02, 0.01]);
    box off;
    hold on;
    if fitting_flag
        
        dataVector = data_FRET{i};
    %         dataVector = sel_FRET_B{i,j};
        dataVector = dataVector(dataVector>-0.5);
        scale = -0.2:0.05:1.2;
        [numCounts, edges] = histcounts(dataVector, scale);
        binCenters = (edges(1:end-1) + edges(2:end))/2; % Calculate bin centers for plotting and fitting

        % Double Gaussian function
        doubleGaussian = @(p, x) p(1)*exp(-((x-p(2))/p(3)).^2) + p(4)*exp(-((x-0.75)/closedFRET_sigma).^2);
    %         singleGaussian = @(p, x) p(1)*exp(-((x-0.75)/closedFRET_sigma).^2);
    %         singleGaussian = @(p, x) p(1)*exp(-((x-p(2))/p(3)).^2);

        % Initial guesses for the parameters: [Amplitude1, Mean1, Sigma1, Amplitude2, Mean2, Sigma2]
        % initialGuess = [1, 0.35, 0.1, 1, 0.75, 0.1];
        initialGuess_2 = [max(numCounts)/2, 0.4, 0.1, max(numCounts)/2];
    %         initialGuess_1 = [1, 0.75, 0.1];
    %         initialGuess_1 = 1;
        % Use lsqcurvefit to fit the model to the data
        options = optimoptions('lsqcurvefit', 'Algorithm', 'trust-region-reflective', 'MaxFunctionEvaluations', 100000, 'Display', 'off');
    %         [pFit_1, chi2_1] = lsqcurvefit(singleGaussian, initialGuess_1, binCenters, numCounts, [], [], options);
        [pFit_2, chi2_2] = lsqcurvefit(doubleGaussian, initialGuess_2, binCenters, numCounts, [0,0,0.02,0], [inf, 0.55, 0.3, inf], options);
        % Generate data from the fit for plotting
    %         fitY_1 = singleGaussian(pFit_1, binCenters);
        fitY_2 = doubleGaussian(pFit_2, binCenters);

    % %         Plot the histogram and the fit
    %     subplot(16,16, plot_num);
    %     histogram(dataVector, scale);
    %     box off
    %     set(gca,'YTickLabel',[]);
    %     set(gca,'XTickLabel',[]);
    %     set(gca,'XTick',[0 0.5 1]);
    %     set(gca,'YTick',[]);
    %     set(gca,'TickDir','out');
    %     ax=gca;
    %     axes('position',ax.Position,'box','on','ytick',[],'xtick',[],'color','none');
    %         figure;
    %         bar(binCenters, numCounts, 'FaceColor', [0.7, 0.7, 0.7]); % Histogram
    %         plot(binCenters, fitY_1, 'r-', 'LineWidth', 2); % Fit
        plot(binCenters, fitY_2, 'r-', 'LineWidth', 2); % Fit
        % Parameters to display
        str = sprintf('Low FRET = %g, sigma = %g, K = %g', pFit_2(2), pFit_2(3), (pFit_2(1)*pFit_2(3))/(pFit_2(4)*closedFRET_sigma));
        
        % Add text to the plot
        text(-0.1, max(fitY_2)*0.8, str, 'FontSize', 6, 'BackgroundColor', 'none');
    end
    % grid minor;
    xlabel('FRET', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial', 'Color', 'k');
    ylabel('Count', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial', 'Color', 'k');
    title([strrep(filePrefix, '_', ' '), num2str(i)], 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');

    % Adjusting axis limits and tick marks
    xlim([-0.2, 1.2]);
%     ylim([0.47, 0.68]);
    xticks(0:0.5:1);
    hold off
%     yticks(0.5:0.05:0.65);

%     savefig([path1 '\histogram_' num2str(i)]);

    baseFileName = [filePrefix, cleanFileName(figTitle),'_'];
    % savefig([path1 '\average_plot']);
    % Add timestamp to base file name
    figFileName = fullfile(saveFolder, [baseFileName, timestamp, '.fig']);
    pdfFileName = fullfile(saveFolder, [baseFileName, timestamp, '.pdf']);

    % Save figure as .fig
    saveas(fig, figFileName);

    % Save figure as .pdf
    exportgraphics(fig,pdfFileName,'ContentType','vector');
end
disp('Done!');

%% dG vs T

x = 10:5:20;

y = [dG_vs_T_10C, dG_vs_T_15C, dG_vs_T_20C];
y_err = [dG_vs_T_10C_SEM, dG_vs_T_15C_SEM, dG_vs_T_20C_SEM];

axesPosition = [0.2 0.2 0.6 0.6];
fig = figure('Name','LinePlot','Units', 'pixels', 'Position', [100 300 240 180]);
axes('Position', axesPosition);
% x_plot = 0:n_points-1;
% axes('OuterPosition', [0 0 1 1]); % Adjust the position of the axes
% plot(x_plot, data_mean, 'k-', 'LineWidth', 2,'Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
errorbar(x, y, y_err, 'k-', 'LineWidth', 2,'Marker', 'o', 'MarkerSize', 4);

% Adjusting plot appearance
set(gca, 'FontSize', 14, 'FontName', 'Arial', 'TickDir', 'out', 'LineWidth', 2);
set(gca, 'XColor', 'k', 'YColor', 'k');
set(gcf, 'color', 'w');
set(gca,'TickLength',[0.02, 0.01]);
box off;

% grid minor;
% xlabel('# of consecutive mismatches', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial', 'Color', 'k');
% ylabel('dG (kcal/mol)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial', 'Color', 'k');
ylabel('dG (kcal/mol)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial', 'Color', 'k');
% figTitle = 'dG vs base pairs';
figTitle = 'dG vs T';
% title(strrep(filePrefix, '_', ' '), 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');

% Adjusting axis limits and tick marks
xlim([7.5, 22.5]);
ylim([-0.5, 2]);
xticks(10:5:20);
yticks(-1:1:3);

if ~exist('saveFolder','var')
    saveFolder = uigetdir('','Select a folder to save figures in');
end
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
% Define base file name
baseFileName = [cleanFileName(figTitle),'_'];
% savefig([path1 '\average_plot']);
% Add timestamp to base file name
figFileName = fullfile(saveFolder, [baseFileName, timestamp, '.fig']);
pdfFileName = fullfile(saveFolder, [baseFileName, timestamp, '.pdf']);

% Save figure as .fig
saveas(fig, figFileName);

% Save figure as .pdf
exportgraphics(fig,pdfFileName,'ContentType','vector');
disp('Done!');