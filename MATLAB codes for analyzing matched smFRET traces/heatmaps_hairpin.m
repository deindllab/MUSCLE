%% Plotting FRET after
% code = 'ATGC';
sel_FRET = sel_FRET_B;
% sel_FRET = sel_FRET_A;
% x_axis = 1:16;
% y_axis = 1:16;
htmap = NaN(16,16);
len = zeros(16,16);

% % idx = [18:21,23,24];
FRET_threshold = 0.5; % midpoint between 0.15 and 0.85 for the hairpin

min_length = 10;

for i = 1:16
    for j = 1:16
        F_t = sel_FRET{i,j};
        F_t = F_t(F_t>-0.2&F_t<1.2);
        temp_n = length(F_t);
        len(i,j) = temp_n;
        if temp_n>min_length
%             htmap(i,j) = nanmean(F_t);
%             temp_std = nanstd(F_t);
            htmap(i,j) = mean(F_t, "omitnan");
            temp_std = std(F_t, "omitnan");
        end
    end
end
%% Create custom htmap code if needed
code = code_htmap;
%% Plot coverage heatmap
% figure();
% heatmap(htmap,'ColorLimits',[0.4 0.8]);
% % figure();
% heatmap(len, 'ColorLimits',[0 2000]);
if ~exist('saveFolder','var')
    saveFolder = uigetdir('','Select a folder to save figures in');
end
max_num = prctile(len(:),95,'all');
% Determine the order of magnitude of the number
orderOfMagnitude = floor(log10(abs(max_num)));

% Calculate the rounding factor
roundingFactor = 10^orderOfMagnitude;

% Round the number down
htmap_scale =  floor(max_num / roundingFactor) * roundingFactor;

% MUSCLE_heatmap(len,[1400 htmap_scale], code, 'hmTitle', 'Coverage after (# time points)', 'saveFolder', saveFolder,'cmap', flipud(autumn));
MUSCLE_heatmap(len,[600 20000], code, 'hmTitle', 'Coverage (# time points)', 'saveFolder', saveFolder,'cmap', flipud(autumn));
% MUSCLE_heatmap(len,[0 htmap_scale], code, 'hmTitle', 'Coverage before (# time points)', 'boxLabelsFontSize', 6, 'saveFolder', saveFolder,'cmap', flipud(autumn));
%% Plot FRET heatmap
if ~exist('saveFolder','var')
    saveFolder = uigetdir('','Select a folder to save figures in');
end
% MUSCLE_heatmap(htmap,[0.15 0.5], code, 'hmTitle', 'Mean FRET with Cas9', 'boxLabelsFontSize', 7, 'saveFolder', saveFolder);
MUSCLE_heatmap(htmap,[0.15 0.75], code, 'hmTitle', 'Mean FRET', 'boxLabelsFontSize', 7, 'saveFolder', saveFolder);
MUSCLE_heatmap(htmap,[0.15 0.75], code, 'hmTitle', 'Mean FRET', 'saveFolder', saveFolder);

%% K and dG based on mean FRET

mean_FRET_all=NaN(16,16);
mean_FRET_all_err=NaN(16,16);

for k = 1:16
    for l = 1:16 
       mean_FRET_all(k,l) = mean(sel_FRET_B{k,l});
       mean_FRET_all_err(k,l)=std(sel_FRET_B{k,l})/sqrt(htmap_num_traces(k,l));        
    end     
end
%%

%Determining the FRET value detrmining the open state 
[min_val, min_idx] = min(mean_FRET_all(:));
[row, col] = ind2sub(size(mean_FRET_all), min_idx); % identifying the indicies of the minimal element on the matrix
dataVector = sel_FRET_B{row,col};
dataVector = dataVector(dataVector>-0.5);
scale = -0.2:0.05:1.2;
[numCounts, edges] = histcounts(dataVector, scale, 'Normalization', 'pdf');
binCenters = (edges(1:end-1) + edges(2:end))/2; % Calculate bin centers for plotting and fitting
singleGaussian = @(p, x) p(1)*exp(-((x-p(2))/p(3)).^2);

initialGuess_2 = [1, 0.15, 0.2];
%         initialGuess_1 = [1, 0.75, 0.1];
        % Use lsqcurvefit to fit the model to the data
options = optimoptions('lsqcurvefit', 'Algorithm', 'trust-region-reflective', 'MaxFunctionEvaluations', 100000, 'Display', 'off');
[pFit_1, chi2_1] = lsqcurvefit(singleGaussian, initialGuess_2, binCenters, numCounts, [0,0.02,0], [inf, 0.22, 0.3], options);
 % Generate data from the fit for plotting
fitY_1 = singleGaussian(pFit_1, binCenters);



figure;
bar(binCenters, numCounts, 'FaceColor', [0.7, 0.7, 0.7]); % Histogram
hold on;
plot(binCenters, fitY_1, 'r-', 'LineWidth', 2); % Fit
xlabel('Data Values');
ylabel('Probability Density');
title('Histogram and Single Gaussian Fit');
legend('Histogram', 'Single Gaussian Fit');

FRET_open = pFit_1(2);
sigma_FRET_open = pFit_1(3);

%%
%Determining the FRET value detrmining the closed state 

% [max_val, max_idx] = max(mean_FRET_all(:));
% [row_max, col_max] = ind2sub(size(mean_FRET_all), max_idx); % identifying the indicies of the minimal element on the matrix
row_max = 16;
col_max = 16;


dataVector = sel_FRET_B{row_max,col_max};
dataVector = dataVector(dataVector>-0.5);
scale = -0.2:0.05:1.2;
[numCounts, edges] = histcounts(dataVector, scale, 'Normalization', 'pdf');
binCenters = (edges(1:end-1) + edges(2:end))/2; % Calculate bin centers for plotting and fitting
      
% doubleGaussian = @(p, x) p(1)*exp(-((x-p(2))/p(3)).^2) + p(4)*exp(-((x-FRET_open)/sigma_FRET_open ).^2)
doubleGaussian = @(p, x) p(1)*exp(-((x-p(2))/p(3)).^2) + p(4)*exp(-((x-p(5))/p(6) ).^2)
% initialGuess_2 = [1, 0.7, 0.2, 1];
initialGuess_2 = [1, 0.7, 0.2, 1,0.1,0.1];
options = optimoptions('lsqcurvefit', 'Algorithm', 'trust-region-reflective', 'MaxFunctionEvaluations', 100000, 'Display', 'off');
% [pFit_2, chi2_2] = lsqcurvefit(doubleGaussian, initialGuess_2, binCenters, numCounts, [0,0.4,0.02,0], [inf, 0.9, 0.3, inf], options);
[pFit_2, chi2_2] = lsqcurvefit(doubleGaussian, initialGuess_2, binCenters, numCounts, [0,0.4,0.02,0,0.1, 0.02], [inf, 1, 0.3, inf, 0.4, 0.3], options);
fitY_2 = doubleGaussian(pFit_2, binCenters);

figure;
bar(binCenters, numCounts, 'FaceColor', [0.7, 0.7, 0.7]); % Histogram
hold on;
plot(binCenters, fitY_2, 'g-', 'LineWidth', 2); % Fit
xlabel('Data Values');
ylabel('Probability Density');
title('Histogram and Double Gaussian Fit');
legend('Histogram', 'Double Gaussian Fit');

FRET_close = pFit_2(2);
sigma_FRET_close = pFit_2(3);

%% FRET open and close based on 10deg
FRET_close = 0.9185;
FRET_open = 0.1609;

%%
if ~exist('saveFolder','var')
    saveFolder = uigetdir('','Select a folder to save figures in');
end

%Calculationg K based on the mean FRET
mean_FRET_all_norm = NaN(16,16);
mean_FRET_all_norm_err = NaN(16,16);
htmap_K = NaN(16,16);
htmap_K_err = NaN(16,16);

R = 1.987e-3; % Universal gas constant for dG calculations in kcal*K-1*mol-1
% T = 283; %10deg Celsius in K
T = 293; %20deg Celsius in K
% 
% figure('Name','Average FRET 10deg');
% heatmap(mean_FRET_all,'ColorLimits',[prctile(mean_FRET_all,5,'all') prctile(mean_FRET_all,95,'all')]);



for k = 1:16
    for l = 1:16 
           mean_FRET_all_norm(k,l) =  (mean_FRET_all(k,l) - FRET_open)/(FRET_close- FRET_open);
           mean_FRET_all_norm_err(k,l) = mean_FRET_all_err(k,l)/(FRET_close- FRET_open);
    end     
end

htmap_K = mean_FRET_all_norm ./ (1 - mean_FRET_all_norm);
htmap_K_err = mean_FRET_all_norm_err./((1-mean_FRET_all_norm).^2);


% MUSCLE_ap(htmap_K,[0 2.5] , code_htmap, 'hmTitle', 'Equilibrium constant 10deg from mean FRET', 'boxLabelsFontSize', 7, 'saveFolder', saveFolder);
% htmap_K(htmap_K == 0) = NaN; % to avoid inf dG
% htmap_dG = -R*T*log(htmap_K);
% MUSCLE_heatmap(htmap_dG,[-0.5 2] , code_htmap, 'hmTitle', 'dG 10deg from mean FRET', 'boxLabelsFontSize', 7, 'saveFolder', saveFolder);

MUSCLE_heatmap(htmap_K,[0 2.5] , code_htmap, 'hmTitle', 'Equilibrium constant 10deg from mean FRET', 'saveFolder', saveFolder);
htmap_K(htmap_K == 0) = NaN; % to avoid inf dG
htmap_dG = -R*T*log(htmap_K);
htmap_dG_err = -R*T*htmap_K_err./htmap_K;
MUSCLE_heatmap(htmap_dG,[-0.5 3] , code_htmap, 'hmTitle', 'dG 10deg from mean FRET', 'saveFolder', saveFolder);
% MUSCLE_heatmap(htmap_dG_err,[-0.5 2] , code_htmap, 'hmTitle', 'dG 10deg from mean FRET', 'boxLabelsFontSize', 7);


