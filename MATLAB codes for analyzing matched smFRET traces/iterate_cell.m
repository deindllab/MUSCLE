%% Initialize 

% Dispersion of the closed state is fixed to the value obtained for a
% target sequence that is almost always in the closed state.
closedFRET_sigma = 0.1658;
% select any cell array
% filePrefix = 'Mismatch_ident_vs_nonident_';
% filePrefix = 'Mismatch_consec_';

filePrefix = 'Outlier_seq_after_3G_GATC_';
cell = mismatch;
sel_FRET = sel_FRET_A;
% sel_FRET = sel_FRET_B;
% cell = mismatch;
clear data
n_points = length(cell);
data{n_points} = [];
axesPosition = [0.2 0.2 0.6 0.6];
% path1 = uigetdir('D:\Experiments\MUSCLE DONUTS\','Choose a directory to save the figures');

for i = 1:n_points
    temp = cell{i};
    data{i} = [];
    for j = 1:size(temp,1)
        x_plot = temp(j,1);
        y_plot = temp(j,2);
        t = sel_FRET(x_plot,y_plot);
        t = t{1};
        t = t(t>-inf);
        data{i} = [data{i} transpose(t)];
    end
end
%% SEM
n_bootstraps = 100;
data_mean = zeros(n_points,1);
lowerCIs = zeros(n_points,1);
upperCIs = zeros(n_points,1);
errors = zeros(n_points,1);
for i = 1:n_points
%     i
    temp = data{i};
    temp = temp (temp>0&temp<1);
        % Compute bootstrapped confidence intervals
%     cis = bootci(n_bootstraps, {@mean, temp}, 'Alpha', 0.01);
    
    SEM = std(temp)/sqrt(length(temp));
    mean_FRET = mean(temp);
%     % Store the confidence intervals
%     lowerCIs(i) = cis(1);
%     upperCIs(i) = cis(2);
%     
    lowerCIs(i) = mean_FRET-SEM;
    upperCIs(i) = mean_FRET+SEM;
    
    data_mean(i) = mean_FRET;
end

%% Save ident
data_mean_ident = data_mean;
lowerCIs_ident = lowerCIs;
upperCIs_ident = upperCIs;


%% Save non-ident
data_mean_nonident = data_mean;
lowerCIs_nonident = lowerCIs;
upperCIs_nonident = upperCIs;


%% Plot histograms
fitting_flag = true;
if ~exist('saveFolder','var')
    saveFolder = uigetdir('','Select a folder to save figures in');
end
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
for i = 1:n_points
    figTitle = ['Hist_',num2str(i)];
    fig = figure('Name',figTitle,'Units', 'pixels', 'Position', [100 300 240 180]);
    axes('Position', axesPosition);
    histogram(data{i}, -0.2:0.05:1.2 , 'FaceColor', [0 174 239] / 255, 'EdgeColor', 'k', 'LineWidth', 1);
    % Adjusting plot appearance
    set(gca, 'FontSize', 14, 'FontName', 'Arial', 'TickDir', 'out', 'LineWidth', 2);
    set(gca, 'XColor', 'k', 'YColor', 'k');
    set(gcf, 'color', 'w');
    set(gca,'TickLength',[0.02, 0.01]);
    box off;
    hold on;
    if fitting_flag
        
        dataVector = data{i};
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
%% Equilibrium constant

% select any cell array
% filePrefix = 'Mismatch_num_ident_vs_nonident_K_';
filePrefix = 'Mismatch_consec_K_';
% filePrefix = 'Example_seq_before_';
cell = mismatch;
% sel_FRET = sel_FRET_A;
% sel_FRET = sel_FRET_B;
% cell = mismatch;
clear data
n_points = length(cell);
data{n_points} = [];
axesPosition = [0.2 0.2 0.6 0.6];
% path1 = uigetdir('D:\Experiments\MUSCLE DONUTS\','Choose a directory to save the figures');

for i = 1:n_points
    temp = cell{i};
    data{i} = [];
    for j = 1:size(temp,1)
        x_plot = temp(j,1);
        y_plot = temp(j,2);

        data{i} = [data{i} htmap_K(x_plot,y_plot)];
    end
end

%% SEM Equilibrium constant

n_bootstraps = 100;
data_mean = zeros(n_points,1);
lowerCIs = zeros(n_points,1);
upperCIs = zeros(n_points,1);
errors = zeros(n_points,1);
for i = 1:n_points
%     i
    temp = data{i};
%     temp = temp (temp>0&temp<1);
        % Compute bootstrapped confidence intervals
%     cis = bootci(n_bootstraps, {@mean, temp}, 'Alpha', 0.01);
    
    SEM = std(temp)/sqrt(length(temp));
    mean_FRET = mean(temp);
%     % Store the confidence intervals
%     lowerCIs(i) = cis(1);
%     upperCIs(i) = cis(2);
%     
    lowerCIs(i) = mean_FRET-SEM;
    upperCIs(i) = mean_FRET+SEM;
    
    data_mean(i) = mean_FRET;
end

%% Plot equilibrium constant

errors = [(data_mean - lowerCIs)'; (upperCIs - data_mean)'];
fig = figure('Name','LinePlot','Units', 'pixels', 'Position', [100 300 240 180]);
axes('Position', axesPosition);
x_plot = 0:n_points-1;
% axes('OuterPosition', [0 0 1 1]); % Adjust the position of the axes
% plot(x_plot, data_mean, 'k-', 'LineWidth', 2,'Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
errorbar(x_plot, data_mean, errors(1,:),errors(2,:), 'k-', 'LineWidth', 2,'Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');

% Adjusting plot appearance
set(gca, 'FontSize', 14, 'FontName', 'Arial', 'TickDir', 'out', 'LineWidth', 2);
set(gca, 'XColor', 'k', 'YColor', 'k');
set(gcf, 'color', 'w');
set(gca,'TickLength',[0.02, 0.01]);
box off;

% grid minor;
% xlabel('# of consecutive mismatches', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial', 'Color', 'k');
ylabel('K', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial', 'Color', 'k');
figTitle = 'K vs consec MM';
title(strrep(filePrefix, '_', ' '), 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');

% Adjusting axis limits and tick marks
xlim([-0.5, 4.5]);
ylim([0, 8]);
xticks(min(x_plot):1:max(x_plot));
yticks(0:2:8);

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

%%
for i = 1:n_points
    figure('Name',strcat(num2str(i),'_trace'));
    plot(data{i});
    ylim([0 1.2]);
%     savefig([path1 '\traces_' num2str(i)]);
end

%% HMM 

% select any cell array
filePrefix = 'Mismatch_num_';
cell = mismatch;
sel_FRET = HMM_kclose;
% sel_FRET = HMM_kopen;
% cell = mismatch;
clear data
n_points = length(cell);
data{n_points} = [];
axesPosition = [0.2 0.2 0.6 0.6];
% path1 = uigetdir('D:\Experiments\MUSCLE DONUTS\','Choose a directory to save the figures');

for i = 1:n_points
    temp = cell{i};
    data{i} = [];
    for j = 1:size(temp,1)
        x_plot = temp(j,1);
        y_plot = temp(j,2);
        t = sel_FRET(x_plot,y_plot);
%         t = t{1};
        data{i} = [data{i} t];
    end
end
% n_bootstraps = 100;
data_mean = zeros(n_points,1);
errors = zeros(n_points,1);
% upperCIs = zeros(n_points,1);
% errors = zeros(n_points,1);
for i = 1:n_points
%     i
    temp = data{i};
%     temp = temp (temp>0&temp<1);
        % Compute bootstrapped confidence intervals
%     cis = bootci(n_bootstraps, {@mean, temp}, 'Alpha', 0.01);
    
    % Store the confidence intervals
    errors(i) = nanstd(temp)/sqrt(numel(temp));
%     upperCIs(i) = cis(2);
    
%     data_mean(i) = nanmean(temp);
    data_mean(i) = mean(temp, "omitnan");
    
end

% errors = [(data_mean - lowerCIs)'; (upperCIs - data_mean)'];
fig = figure('Name','LinePlot');
axes('Position', axesPosition);
x_plot = 0:n_points-1;
% axes('OuterPosition', [0 0 1 1]); % Adjust the position of the axes
% plot(x_plot, data_mean, 'k-', 'LineWidth', 2,'Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
errorbar(x_plot, data_mean, errors,errors, 'k-', 'LineWidth', 2,'Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');

% Adjusting plot appearance
set(gca, 'FontSize', 14, 'FontName', 'Arial', 'TickDir', 'out', 'LineWidth', 2);
set(gca, 'XColor', 'k', 'YColor', 'k');
set(gcf, 'color', 'w');
set(gca,'TickLength',[0.02, 0.01]);
box off;

% grid minor;
% xlabel('# of consecutive mismatches', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial', 'Color', 'k');
ylabel('HMM krewind', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial', 'Color', 'k');
figTitle = 'HMM krewind';
title(strrep(filePrefix, '_', ' '), 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');

% Adjusting axis limits and tick marks
% xlim([-0.5, 4.5]);
% ylim([0.42, 0.63]);
xticks(min(x_plot):1:max(x_plot));
yticks(0.2:0.2:0.6);

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


%% Histogram for a selected sequence
% sequence = 'CTTG';
% sequence = 'ACAG';
sequence = 'GGCC';
% code = ['GCAT';'ATGC';'TGCA';'GCAT']; % Cas9
xcorr =4*(strfind(code(3,:),sequence(3))-1)+(strfind(code(1,:),sequence(1)));
ycorr =4*(strfind(code(4,:),sequence(4))-1)+(strfind(code(2,:),sequence(2)));
% xcorr = 16;
% ycorr = 2;

filePrefix = ['Hist_before_', sel_seq{xcorr,ycorr},'_'];
if ~exist('saveFolder','var')
    saveFolder = uigetdir('','Select a folder to save figures in');
end
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
% for i = 1:n_points
figTitle = ['Hist_',num2str(sel_seq{xcorr,ycorr})];
fig = figure('Name',figTitle,'Units', 'pixels', 'Position', [100 300 240 180]);
% axes('Position', axesPosition);
histogram(sel_FRET_B{xcorr,ycorr}, -0.2:0.05:1.2 , 'FaceColor', [0 174 239] / 255, 'EdgeColor', 'k', 'LineWidth', 1);
% Adjusting plot appearance
set(gca, 'FontSize', 14, 'FontName', 'Arial', 'TickDir', 'out', 'LineWidth', 2);
set(gca, 'XColor', 'k', 'YColor', 'k');
set(gcf, 'color', 'w');
set(gca,'TickLength',[0.02, 0.01]);
box off;

% grid minor;
xlabel('FRET', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial', 'Color', 'k');
ylabel('Count', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial', 'Color', 'k');
title(sel_seq{xcorr,ycorr}, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');

% Adjusting axis limits and tick marks
xlim([-0.2, 1.2]);
%     ylim([0.47, 0.68]);
xticks(0:0.5:1);
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
% end
disp('Done!');

%% Plot identical vs non-identical mismatches
% errors = [(data_mean - lowerCIs)'; (upperCIs - data_mean)'];
fig = figure('Name','LinePlot','Units', 'pixels', 'Position', [100 300 240 180]);
axes('Position', axesPosition);
x_plot = 0:n_points-1;
% axes('OuterPosition', [0 0 1 1]); % Adjust the position of the axes
hold on

errors_ident = [(data_mean_ident - lowerCIs_ident)'; (upperCIs_ident - data_mean_ident)'];
errors_nonident = [(data_mean_nonident - lowerCIs_nonident)'; (upperCIs_nonident - data_mean_nonident)'];
errorbar(x_plot, data_mean_ident, errors_ident(1,:),errors_ident(2,:), 'r-', 'LineWidth', 2,'Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
errorbar(x_plot, data_mean_nonident, errors_nonident(1,:),errors_nonident(2,:), 'k-', 'LineWidth', 2,'Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
% plot(x_plot, data_mean_ident, 'r-', 'LineWidth', 2,'Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
% plot(x_plot, data_mean_nonident, 'k-', 'LineWidth', 2,'Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
hold off
% errorbar(x_plot, data_mean, errors(1,:),errors(2,:), 'k-', 'LineWidth', 2,'Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');

% Adjusting plot appearance
set(gca, 'FontSize', 14, 'FontName', 'Arial', 'TickDir', 'out', 'LineWidth', 2);
set(gca, 'XColor', 'k', 'YColor', 'k');
set(gcf, 'color', 'w');
set(gca,'TickLength',[0.02, 0.01]);
box off;

% grid minor;
% xlabel('# of consecutive mismatches', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial', 'Color', 'k');
% ylabel('Mean FRET', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial', 'Color', 'k');
% figTitle = 'Mean FRET';
ylabel('K', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial', 'Color', 'k');
figTitle = 'K ident vs nonident';
title(strrep(filePrefix, '_', ' '), 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');

% Adjusting axis limits and tick marks
xlim([-0.5, 4.5]);
% ylim([0, 8]);
xticks(min(x_plot):1:max(x_plot));
% yticks(0:2:8);
ylim([0.35, 0.75]);
yticks(0.35:0.1:0.75);

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

%% Plot histograms Triple Gauss fit
fitting_flag = true;
if ~exist('saveFolder','var')
    saveFolder = uigetdir('','Select a folder to save figures in');
end
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
for i = 1:n_points
    figTitle = ['Hist_',num2str(i)];
    fig = figure('Name',figTitle,'Units', 'pixels', 'Position', [100 300 240 180]);
    axes('Position', axesPosition);
    histogram(data{i}, -0.2:0.05:1.2 , 'FaceColor', [0 174 239] / 255, 'EdgeColor', 'k', 'LineWidth', 1);
    % Adjusting plot appearance
    set(gca, 'FontSize', 14, 'FontName', 'Arial', 'TickDir', 'out', 'LineWidth', 2);
    set(gca, 'XColor', 'k', 'YColor', 'k');
    set(gcf, 'color', 'w');
    set(gca,'TickLength',[0.02, 0.01]);
    box off;
    hold on;
    if fitting_flag
        
        dataVector = data{i};
    %         dataVector = sel_FRET_B{i,j};
        dataVector = dataVector(dataVector>-0.5);
        scale = -0.2:0.05:1.2;
        [numCounts, edges] = histcounts(dataVector, scale);
        binCenters = (edges(1:end-1) + edges(2:end))/2; % Calculate bin centers for plotting and fitting

        % Double Gaussian function
        tripleGaussian = @(p, x) p(1)*exp(-((x-p(2))/p(3)).^2) +  p(4)*exp(-((x-p(5))/p(6)).^2) + p(7)*exp(-((x-p(8))/p(9)).^2);
    %         singleGaussian = @(p, x) p(1)*exp(-((x-0.75)/closedFRET_sigma).^2);
    %         singleGaussian = @(p, x) p(1)*exp(-((x-p(2))/p(3)).^2);

        % Initial guesses for the parameters: [Amplitude1, Mean1, Sigma1, Amplitude2, Mean2, Sigma2]
        % initialGuess = [1, 0.35, 0.1, 1, 0.75, 0.1];
        initialGuess_2 = [max(numCounts)/3, 0.3, 0.1, max(numCounts)/3, 0.5, 0.1, max(numCounts)/3, 0.75, closedFRET_sigma];
    %         initialGuess_1 = [1, 0.75, 0.1];
    %         initialGuess_1 = 1;
        % Use lsqcurvefit to fit the model to the data
        options = optimoptions('lsqcurvefit', 'Algorithm', 'trust-region-reflective', 'MaxFunctionEvaluations', 100000, 'Display', 'off');
    %         [pFit_1, chi2_1] = lsqcurvefit(singleGaussian, initialGuess_1, binCenters, numCounts, [], [], options);
        [pFit_2, chi2_2] = lsqcurvefit(tripleGaussian, initialGuess_2, binCenters, numCounts, [0,0,0.05,0,0.45,0.05,0,0.7, 0.05], [inf, 0.45, 0.3, inf, 0.65, 0.3, inf,0.9, 0.3], options);
        % Generate data from the fit for plotting
    %         fitY_1 = singleGaussian(pFit_1, binCenters);
        fitY_2 = tripleGaussian(pFit_2, binCenters);

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
        str = sprintf('A 1 = %g Low FRET 1 = %g, sigma 1 = %g, A 2 = %g Low FRET 2 = %g, sigma 2 = %g, A 3 = %g High FRET 3 %g sigma 3 %g' , pFit_2(1), pFit_2(2), pFit_2(3), pFit_2(4), pFit_2(5), pFit_2(6), pFit_2(7), pFit_2(8), pFit_2(9));
        
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

%% Plot Cas9 HMM
axesPosition = [0.2 0.2 0.6 0.6];
filePrefix = 'HMM MM kclose';
fig = figure('Name','LinePlot','Units', 'pixels', 'Position', [100 300 240 180]);
axes('Position', axesPosition);

% axes('OuterPosition', [0 0 1 1]); % Adjust the position of the axes
% plot(x_plot, data_mean, 'k-', 'LineWidth', 2,'Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
errorbar(x, kclose_MM, kclose_MM_SEM, 'k-', 'LineWidth', 2,'Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
hold on;
% errorbar(x, kopen_GC, kopen_GC_SEM, 'r-', 'LineWidth', 2,'Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');

% Adjusting plot appearance
set(gca, 'FontSize', 14, 'FontName', 'Arial', 'TickDir', 'out', 'LineWidth', 2);
set(gca, 'XColor', 'k', 'YColor', 'k');
set(gcf, 'color', 'w');
set(gca,'TickLength',[0.02, 0.01]);
box off;

% grid minor;
% xlabel('# of consecutive mismatches', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial', 'Color', 'k');
ylabel('Rate (s-1)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial', 'Color', 'k');
figTitle = filePrefix;
title(filePrefix, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');

% Adjusting axis limits and tick marks
xlim([-0.5, 4.5]);
ylim([0, 0.6]);
xticks(0:4);
yticks([0:0.2:0.6]);

if ~exist('saveFolder','var')
    saveFolder = uigetdir('','Select a folder to save figures in');
end
hold off;
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


%% HMM vs FRET scatter plots

axesPosition = [0.2 0.2 0.6 0.6];
filePrefix = 'FRET vs kopen';
fig = figure('Name','LinePlot','Units', 'pixels', 'Position', [100 300 240 180]);
axes('Position', axesPosition);

% axes('OuterPosition', [0 0 1 1]); % Adjust the position of the axes
scatter(FRET, kopen, 'k-','Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'k', 'MfarkerEdgeColor', 'k');
% errorbar(x, kclose_MM, kclose_MM_SEM, 'k-', 'LineWidth', 2,'Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
hold on;
% errorbar(x, kopen_GC, kopen_GC_SEM, 'r-', 'LineWidth', 2,'Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');

% Adjusting plot appearance
set(gca, 'FontSize', 14, 'FontName', 'Arial', 'TickDir', 'out', 'LineWidth', 2);
set(gca, 'XColor', 'k', 'YColor', 'k');
set(gcf, 'color', 'w');
set(gca,'TickLength',[0.02, 0.01]);
box off;

% grid minor;
% xlabel('# of consecutive mismatches', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial', 'Color', 'k');
ylabel('Rate (s-1)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial', 'Color', 'k');
figTitle = filePrefix;
title(filePrefix, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');

% Adjusting axis limits and tick marks
xlim([-0.5, 4.5]);
ylim([0, 0.6]);
xticks(0:4);
yticks([0:0.2:0.6]);

if ~exist('saveFolder','var')
    saveFolder = uigetdir('','Select a folder to save figures in');
end
hold off;
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

%% NN vs measurement scatter plots

axesPosition = [0.2 0.2 0.6 0.6];
filePrefix = 'dG measured vs NN';
fig = figure('Name','LinePlot','Units', 'pixels', 'Position', [100 300 240 180]);
axes('Position', axesPosition);

% axes('OuterPosition', [0 0 1 1]); % Adjust the position of the axes
scatter(htmap_dG(:), deltaG_NN_htmap(:), 'k-','Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'k', 'MfarkerEdgeColor', 'k');
% errorbar(x, kclose_MM, kclose_MM_SEM, 'k-', 'LineWidth', 2,'Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
hold on;
% errorbar(x, kopen_GC, kopen_GC_SEM, 'r-', 'LineWidth', 2,'Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');

% Adjusting plot appearance
set(gca, 'FontSize', 14, 'FontName', 'Arial', 'TickDir', 'out', 'LineWidth', 2);
set(gca, 'XColor', 'k', 'YColor', 'k');
set(gcf, 'color', 'w');
set(gca,'TickLength',[0.02, 0.01]);
box off;

% grid minor;
% xlabel('# of consecutive mismatches', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial', 'Color', 'k');
ylabel('Rate (s-1)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial', 'Color', 'k');
figTitle = filePrefix;
title(filePrefix, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');

% Adjusting axis limits and tick marks
xlim([-0.5, 1.5]);
ylim([-9, 3]);
xticks(-0.5:0.5:1.5);
yticks([-9:3:3]);

if ~exist('saveFolder','var')
    saveFolder = uigetdir('','Select a folder to save figures in');
end
hold off;
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