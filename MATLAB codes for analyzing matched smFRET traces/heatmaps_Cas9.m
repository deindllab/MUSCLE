%% Plotting FRET before

sel_FRET = sel_FRET_B;
htmap = NaN(16,16);
htmap_K = NaN(16,16);
len = zeros(16,16);
FRET_threshold = 0.57; % midpoint between 0.4 and 0.74 for Cas9 is 0.57
min_length = 10;

for i = 1:16
    for j = 1:16
        F_t = sel_FRET{i,j};
        F_t = F_t(F_t>-0.2&F_t<1.2);
        temp_n = length(F_t);
        len(i,j) = temp_n;
        if temp_n>min_length
            %             all_FRET = [all_FRET F_t];
%             htmap(i,j) = nanmean(F_t);
            htmap(i,j) = mean(F_t, "omitnan");
%             temp_std = nanstd(F_t);
            temp_std = std(F_t, "omitnan");
            htmap_K(i,j) = sum(F_t<FRET_threshold)/sum(F_t>FRET_threshold);
%             htmap_K(i,j) = sum(F_t>FRET_threshold)/sum(F_t<FRET_threshold);
            % Calculate the t-value for a 99% confidence interval
        end
    end
end

%% Plot coverage heatmap
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


MUSCLE_heatmap(len,[1400 6000], code, 'hmTitle', 'Coverage (# time points)', 'saveFolder', saveFolder,'cmap', flipud(autumn));
% MUSCLE_heatmap(len,[0 htmap_scale], code, 'hmTitle', 'Coverage before (# time points)', 'boxLabelsFontSize', 6, 'saveFolder', saveFolder,'cmap', flipud(autumn));
%% Plot FRET heatmap Before
if ~exist('saveFolder','var')
    saveFolder = uigetdir('','Select a folder to save figures in');
end

% MUSCLE_heatmap(htmap,[0.4 0.75], code, 'hmTitle', 'Mean FRET', 'boxLabelsFontSize', 7, 'saveFolder', saveFolder);
MUSCLE_heatmap(htmap,[0.4 0.75], code, 'hmTitle', 'Mean FRET without Cas9', 'saveFolder', saveFolder);
MUSCLE_heatmap(htmap_K,[0 5], code, 'hmTitle', 'Equilibrium constant without Cas9', 'saveFolder', saveFolder);
R = 1.987e-3; % Universal gas constant for dG calculations in kcal*K-1*mol-1
T = 293; %20deg Celsius in K
htmap_dG = -R*T*log(htmap_K);
MUSCLE_heatmap(htmap_dG,[-1 3] , code, 'hmTitle', 'dG without Cas9', 'boxLabelsFontSize', 7, 'saveFolder', saveFolder);

%% Plotting FRET with Cas9

sel_FRET = sel_FRET_A;
htmap = NaN(16,16);
htmap_K = NaN(16,16);
len = zeros(16,16);
FRET_threshold = 0.57; % midpoint between 0.4 and 0.74 for Cas9 is 0.57
min_length = 10;

for i = 1:16
    for j = 1:16
        F_t = sel_FRET{i,j};
        F_t = F_t(F_t>-0.2&F_t<1.2);
        temp_n = length(F_t);
        len(i,j) = temp_n;
        if temp_n>min_length
            %             all_FRET = [all_FRET F_t];
            %             htmap(i,j) = nanmean(F_t);
            htmap(i,j) = mean(F_t, "omitnan");
%             temp_std = nanstd(F_t);
            temp_std = std(F_t, "omitnan");
            htmap_K(i,j) = sum(F_t<FRET_threshold)/sum(F_t>FRET_threshold);
%             htmap_K(i,j) = sum(F_t>FRET_threshold)/sum(F_t<FRET_threshold);
            % Calculate the t-value for a 99% confidence interval
        end
    end
end

%% Plot coverage heatmap
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


MUSCLE_heatmap(len,[1400 6000], code, 'hmTitle', 'Coverage (# time points)', 'saveFolder', saveFolder,'cmap', flipud(autumn));
% MUSCLE_heatmap(len,[0 htmap_scale], code, 'hmTitle', 'Coverage before (# time points)', 'boxLabelsFontSize', 6, 'saveFolder', saveFolder,'cmap', flipud(autumn));

%% Plot FRET heatmap with Cas9
if ~exist('saveFolder','var')
    saveFolder = uigetdir('','Select a folder to save figures in');
end

% MUSCLE_heatmap(htmap,[0.4 0.75], code, 'hmTitle', 'Mean FRET', 'boxLabelsFontSize', 7, 'saveFolder', saveFolder);
MUSCLE_heatmap(htmap,[0.4 0.75], code, 'hmTitle', 'Mean FRET with Cas9', 'saveFolder', saveFolder);
MUSCLE_heatmap(htmap_K,[0 5], code, 'hmTitle', 'Equilibrium constant with Cas9', 'saveFolder', saveFolder);
R = 1.987e-3; % Universal gas constant for dG calculations in kcal*K-1*mol-1
T = 293; %20deg Celsius in K
htmap_dG = -R*T*log(htmap_K);
MUSCLE_heatmap(htmap_dG,[-1 3] , code, 'hmTitle', 'dG with Cas9', 'boxLabelsFontSize', 7, 'saveFolder', saveFolder);
