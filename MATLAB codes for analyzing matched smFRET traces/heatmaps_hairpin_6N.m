%% Heatmap number of traces

code = 'ATGC';
code_htmap = ['ATGC';'ATGC';'ATGC';'TACG';'TACG';'TACG']; % Hairpin
htmap_num_traces = zeros(64);
idx = 39:41;
idx_2 = 2:4;

for i = 1: length(Seq1)
    s = Seq(i,idx);
    s_2 = Seq_2(i,idx_2);
    invalidChars = [~ismember(s, code) ~ismember(s_2, code)];
    if any(invalidChars)
        continue
    end
    x =16*(strfind(code,s(1))-1)+4*(strfind(code,s(2))-1)+(strfind(code,s(3)));
    y =16*(strfind(code,s_2(1))-1)+4*(strfind(code,s_2(2))-1)+(strfind(code,s_2(3)));
    htmap_num_traces(x,y) = htmap_num_traces(x,y) + 1;
end

%% Plot number of traces
if ~exist('saveFolder','var')
    saveFolder = uigetdir('','Select a folder to save figures in');
end
max_num = prctile(htmap_num_traces(:),98)/2;
% Determine the order of magnitude of the number
orderOfMagnitude = floor(log10(abs(max_num)));

% Calculate the rounding factor
roundingFactor = 10^orderOfMagnitude;

% Round the number down
htmap_scale =  floor(max_num / roundingFactor) * roundingFactor;

% MUSCLE_heatmap_6N(htmap_num_traces,[0 htmap_scale] , code_htmap, 'hmTitle', 'Number of traces', 'saveFolder', saveFolder);
MUSCLE_heatmap_6N(htmap_num_traces,[10 300] , code_htmap, 'hmTitle', 'Number of traces', 'saveFolder', saveFolder, 'cmap', flipud(autumn));
%% Initialization
% code = 'ATGC';
R = 1.987e-3; % Universal gas constant for dG calculations in kcal*K-1*mol-1
T = 293; %10deg Celsius in K
FRET_threshold = 0.5;
htmap = NaN(64,64);
htmap_K = NaN(64,64);
htmap_dG = NaN(64,64);
len = zeros(64,64);
min_length = 10;

%% Plotting FRET
sel_FRET = sel_FRET_B;
%%
for i = 1:64
    for j = 1:64
        len(i,j) = length(sel_FRET{i,j});
%         if length(sel_FRET{i,j})>20
        F_t = sel_FRET{i,j};
        if iscell(F_t)
            F_t = cat(2, F_t{:});
        end
        F_t = F_t(F_t>0.2&F_t<1.2);
%         all_FRET_B = [all_FRET_B F_t];
%         htmap(i,j) = nanmean(F_t);
        htmap(i,j) = mean(F_t, "omitnan");
        htmap_K(i,j) = sum(F_t>FRET_threshold)/sum(F_t<FRET_threshold);
%         end
    end
end
%% Plot FRET heatmaps
if ~exist('saveFolder','var')
    saveFolder = uigetdir('','Select a folder to save figures in');
end
max_num = prctile(len(:),98)/2;
% Determine the order of magnitude of the number
orderOfMagnitude = floor(log10(abs(max_num)));

% Calculate the rounding factor
roundingFactor = 10^orderOfMagnitude;

% Round the number down
htmap_scale =  floor(max_num / roundingFactor) * roundingFactor;

% MUSCLE_heatmap_6N(htmap_num_traces,[0 htmap_scale] , code_htmap, 'hmTitle', 'Number of traces', 'saveFolder', saveFolder);
MUSCLE_heatmap_6N(len,[100 3000] , code_htmap, 'hmTitle', 'Number of time points', 'saveFolder', saveFolder, 'cmap', flipud(autumn));

MUSCLE_heatmap_6N(htmap,[0.35 0.65] , code_htmap, 'hmTitle', 'Mean FRET', 'saveFolder', saveFolder);

MUSCLE_heatmap_6N(htmap_K,[0 2] , code_htmap, 'hmTitle', 'Equilibrium constant', 'saveFolder', saveFolder);

htmap_dG = -R*T*log(htmap_K);
MUSCLE_heatmap_6N(htmap_dG,[-1 1.5] , code_htmap, 'hmTitle', 'dG', 'saveFolder', saveFolder);