
%These scripts have been modified to make the histograms look nice

% Dispersion of the closed state is fixed to the value obtained for a
% target sequence that is almost always in the closed state.
closedFRET_sigma = 0.1658;
%%%%% COVERAGE %%%%%%%%%%%%
% Define the sequence codes
code = ['GCAT';'ATGC';'TGCA';'GCAT']; % Cas9

% Define the x and y axis ranges
x_axis = 1:16;
y_axis = 1:16;

% Initialize the heatmap matrix with zeros
htmap = zeros(16,16);

% Define the indices of the sequences to be analyzed
idx = 27:30; % Cas9

% Extract the selected sequences from the full dataset
Seq1 = Seq(:,idx);

% Define Threshold for the number of sequences

threshold = 5;

% Iterate over the selected sequences

for i = 1:length(Seq1)
   % Extract the current sequence
   s = Seq1(i,:);
   
   % Map the current sequence to x,y coordinates on the heatmap
   x = 4*(strfind(code(4,:),s(4))-1)+(strfind(code(3,:),s(3))); % Hairpin
   y = 4*(strfind(code(1,:),s(1))-1)+(strfind(code(2,:),s(2)));
    
%      x = 4*(strfind(code(3,:),s(3))-1)+(strfind(code(1,:),s(1))); % Cas9
%      y = 4*(strfind(code(4,:),s(4))-1)+(strfind(code(2,:),s(2)));

     htmap(x,y) = htmap(x,y) + 1;

   
end

numberoftraces = htmap;
%% Plotting an array of histograms

sel_FRET = sel_FRET_B;

scale = -0.2:0.05:1.2;
plot_num = 1;

figure('Position', [100, 100, 1000, 800]);

for j = 1:16
    for i = 1:16
        
        F_t = sel_FRET{i,j};
        F_t = F_t(~isnan(F_t));
        F_t = F_t(F_t>-0.2);
%         if numberoftraces(i,j) >= threshold
           subplot(16,16, plot_num);
           histogram(F_t, scale);
           box off
           set(gca,'YTickLabel',[]);
           set(gca,'XTickLabel',[]);
           set(gca,'XTick',[0 0.5 1]);
           set(gca,'YTick',[]);
           set(gca,'TickDir','out');
           ax=gca;
           axes('position',ax.Position,'box','on','ytick',[],'xtick',[],'color','none');
%         end
        plot_num = plot_num + 1;
    end
end


%% Saving the array of histograms

% Assuming your figure is the current figure
fig = gcf; % Get the handle to the current figure

% Export to PDF
% !Change the output path accordingly!
exportgraphics(fig, 'D:\hist_all_before.pdf', 'ContentType', 'vector');


%% Histogram fitting analysis
htmap_p = NaN(16);
htmap_lowFRETfraction = NaN(16);
htmap_Kfromfit = NaN(16);

% i = 6;
% j = 9;
for i = 1:16
    for j = 1:16
        dataVector = sel_FRET_A{i,j};
%         dataVector = sel_FRET_B{i,j};
        dataVector = dataVector(dataVector>-0.5);
        scale = -0.2:0.05:1.2;
        [numCounts, edges] = histcounts(dataVector, scale, 'Normalization', 'pdf');
        binCenters = (edges(1:end-1) + edges(2:end))/2; % Calculate bin centers for plotting and fitting

        % Double Gaussian function
        doubleGaussian = @(p, x) p(1)*exp(-((x-p(2))/p(3)).^2) + p(4)*exp(-((x-0.75)/closedFRET_sigma).^2);
        singleGaussian = @(p, x) p(1)*exp(-((x-0.75)/closedFRET_sigma).^2);
%         singleGaussian = @(p, x) p(1)*exp(-((x-p(2))/p(3)).^2);

        % Initial guesses for the parameters: [Amplitude1, Mean1, Sigma1, Amplitude2, Mean2, Sigma2]
        % initialGuess = [1, 0.35, 0.1, 1, 0.75, 0.1];
        initialGuess_2 = [1, 0.4, 0.1, 1];
%         initialGuess_1 = [1, 0.75, 0.1];
        initialGuess_1 = 1;
        % Use lsqcurvefit to fit the model to the data
        options = optimoptions('lsqcurvefit', 'Algorithm', 'trust-region-reflective', 'MaxFunctionEvaluations', 100000, 'Display', 'off');
        [pFit_1, chi2_1] = lsqcurvefit(singleGaussian, initialGuess_1, binCenters, numCounts, [], [], options);
        [pFit_2, chi2_2] = lsqcurvefit(doubleGaussian, initialGuess_2, binCenters, numCounts, [0,0,0.02,0], [inf, 0.55, 0.3, inf], options);
        % Generate data from the fit for plotting
        fitY_1 = singleGaussian(pFit_1, binCenters);
        fitY_2 = doubleGaussian(pFit_2, binCenters);

% %         Plot the histogram and the fit
%         figure;
%         bar(binCenters, numCounts, 'FaceColor', [0.7, 0.7, 0.7]); % Histogram
%         hold on;
%         plot(binCenters, fitY_1, 'r-', 'LineWidth', 2); % Fit
%         plot(binCenters, fitY_2, 'g-', 'LineWidth', 2); % Fit
%         xlabel('Data Values');
%         ylabel('Probability Density');
%         title('Histogram and Double Gaussian Fit');
%         legend('Histogram', 'Double Gaussian Fit');

        htmap_p(i,j) = ftest(length(binCenters),1,4,chi2_1,chi2_2);
        htmap_lowFRETvalue(i,j) = pFit_2(2);
        htmap_Kfromfit(i,j) = (pFit_2(1)*pFit_2(3))/(pFit_2(4)*closedFRET_sigma);
    end
end


%% Filter heatmaps
htmap_lowFRETvalue_filt = htmap_lowFRETvalue;
htmap_Kfromfit_filt = htmap_Kfromfit;

idx = htmap_p>1e-5;
htmap_lowFRETvalue_filt(idx) = NaN;
htmap_Kfromfit_filt(idx) = NaN;

idx = htmap_lowFRETvalue>0.5;
htmap_lowFRETvalue_filt(idx) = NaN;
htmap_Kfromfit_filt(idx) = NaN;

idx = htmap_Kfromfit<0.3;
htmap_lowFRETvalue_filt(idx) = NaN;
htmap_Kfromfit_filt(idx) = NaN;

%% Plotting heatmaps
if ~exist('saveFolder','var')
    saveFolder = uigetdir('','Select a folder to save figures in');
end
% MUSCLE_heatmap(htmap_Kfromfit,[0 8], code, 'hmTitle', 'Equilibrium constant Fit', 'boxLabelsFontSize', 7, 'saveFolder', saveFolder);
MUSCLE_heatmap(htmap_lowFRETvalue_filt,[0.3 0.55], code, 'hmTitle', 'Open FRET value Fit', 'boxLabelsFontSize', 7, 'saveFolder', saveFolder);
% MUSCLE_heatmap(htmap_Kfromfit_filt,[0 8], code, 'hmTitle', 'Equilibrium constant Fit', 'saveFolder', saveFolder);
MUSCLE_heatmap(htmap_lowFRETvalue_filt,[0.3 0.55], code, 'hmTitle', 'Open FRET value Fit', 'saveFolder', saveFolder);

%% Plotting an array of histograms with fits
sel_FRET = sel_FRET_A;

scale = -0.2:0.05:1.2;
plot_num = 1;

figure('Position', [100, 100, 1000, 800]);

for j = 1:16
    for i = 1:16
        
       
        dataVector = sel_FRET{i,j};
%         dataVector = sel_FRET_B{i,j};
        dataVector = dataVector(dataVector>-0.5);
        scale = -0.2:0.05:1.2;
        [numCounts, edges] = histcounts(dataVector, scale, 'Normalization', 'pdf');
        binCenters = (edges(1:end-1) + edges(2:end))/2; % Calculate bin centers for plotting and fitting

        % Double Gaussian function
        doubleGaussian = @(p, x) p(1)*exp(-((x-p(2))/p(3)).^2) + p(4)*exp(-((x-0.75)/closedFRET_sigma).^2);
%         singleGaussian = @(p, x) p(1)*exp(-((x-0.75)/closedFRET_sigma).^2);
%         singleGaussian = @(p, x) p(1)*exp(-((x-p(2))/p(3)).^2);

        % Initial guesses for the parameters: [Amplitude1, Mean1, Sigma1, Amplitude2, Mean2, Sigma2]
        % initialGuess = [1, 0.35, 0.1, 1, 0.75, 0.1];
        initialGuess_2 = [1, 0.4, 0.1, 1];
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
        subplot(16,16, plot_num);
        histogram(dataVector, scale);
        box off
        set(gca,'YTickLabel',[]);
        set(gca,'XTickLabel',[]);
        set(gca,'XTick',[0 0.5 1]);
        set(gca,'YTick',[]);
        set(gca,'TickDir','out');
        ax=gca;
        axes('position',ax.Position,'box','on','ytick',[],'xtick',[],'color','none');
%         figure;
%         bar(binCenters, numCounts, 'FaceColor', [0.7, 0.7, 0.7]); % Histogram
        hold on;
%         plot(binCenters, fitY_1, 'r-', 'LineWidth', 2); % Fit
        plot(binCenters, fitY_2, 'r-', 'LineWidth', 2); % Fit
%         xlabel('Data Values');
%         ylabel('Probability Density');
%         title('Histogram and Double Gaussian Fit');
%         legend('Histogram', 'Double Gaussian Fit');
        
%         if numberoftraces(i,j) >= threshold
           
%         end
        plot_num = plot_num + 1;
    end
end

%% Saving the array of histograms

% Assuming your figure is the current figure
fig = gcf; % Get the handle to the current figure

% Export to PDF
% !Change the output path accordingly!
exportgraphics(fig, 'D:\hist_all_after_fit.pdf', 'ContentType', 'vector');