window = 1;
threshold = 0.02; % Threshold for Cy5 direct relative to total intensity with green excitation
clear sel_FRET sel_FRET_A sel_FRET_B sel_Cy3 sel_Cy5
code = ['GCAT';'ATGC';'TGCA';'GCAT']; % Cas9
idx = 27:30; % Ha
reject_norm = [];
n_traces = length(Seq);
max_int = zeros(n_traces,1);
max_int_all = zeros(n_traces,1);
min_int_all = zeros(n_traces,1);
sel_FRET_A{16,16} = [];
sel_Cy3_A{16,16} = [];
sel_Cy5_A{16,16} = [];
sel_FRET_B{16,16} = [];
sel_Cy3_B{16,16} = [];
sel_Cy5_B{16,16} = [];
sel_delta_FRET{16,16} = [];
%start_frame = 41;
start_frame = 10;

% preallocating the sel_FRET

max_length_B = start_frame-1;
max_length_A = size(Cy3,2)-max_length_B;
number_htmap = zeros(16,16);
for i = 1: length(Seq)
   s = Seq(i,idx);
   x =4*(strfind(code(3,:),s(3))-1)+(strfind(code(1,:),s(1)));
   y =4*(strfind(code(4,:),s(4))-1)+(strfind(code(2,:),s(2)));
   number_htmap(x,y) = number_htmap(x,y) + 1;
%    if (x == 18) & (y == 21)
%        weird = [weird; Seq(i,:)];
%    end
end

for i = 1:16
    for j = 1:16
        len_B = max_length_B*number_htmap(i,j);
        len_A = max_length_A*number_htmap(i,j);
        sel_FRET_A{i,j} = NaN(len_A,1);
        sel_Cy3_A{i,j} = NaN(len_A,1);
        sel_Cy5_A{i,j} = NaN(len_A,1);
        sel_FRET_B{i,j} = NaN(len_B,1);
        sel_Cy3_B{i,j} = NaN(len_B,1);
        sel_Cy5_B{i,j} = NaN(len_B,1);
    end
end

htmap_counter = zeros(16);
%% plotting number heatmap

if ~exist('saveFolder','var')
    saveFolder = uigetdir('','Select a folder to save figures in');
end
max_num = max(number_htmap(:))/2;
% Determine the order of magnitude of the number
orderOfMagnitude = floor(log10(abs(max_num)));

% Calculate the rounding factor
roundingFactor = 10^orderOfMagnitude;

% Round the number down
htmap_scale =  floor(max_num / roundingFactor) * roundingFactor;

% MUSCLE_heatmap(number_htmap,[0 htmap_scale] , code, 'hmTitle', 'Number of traces', 'boxLabelsFontSize', 7, 'saveFolder', saveFolder, 'cmap', flipud(autumn));
MUSCLE_heatmap(number_htmap,[150 800] , code, 'hmTitle', 'Number of traces', 'boxLabelsFontSize', 7, 'saveFolder', saveFolder, 'cmap', flipud(autumn));

%% Plotting intensity histograms
Int_hist_G = mean((Cy3(:,2:2:2*start_frame)+Cy5(:,2:2:2*start_frame)),2);

figure('Name','Green intensity');
histogram(Int_hist_G,100, 'BinLimits', [prctile(Int_hist_G,2,'all'), prctile(Int_hist_G,98,'all')]);

Int_hist_R = mean((Cy5(:,1:2:2*start_frame)),2);

figure('Name','Red intensity');
histogram(Int_hist_R,100, 'BinLimits', [prctile(Int_hist_R,2,'all'), prctile(Int_hist_R,98,'all')]);

FRET_hist = mean(Cy5(:,2:2:2*start_frame),2)./Int_hist_G;

figure('Name','FRET');
histogram(FRET_hist,-0.2:0.02:1.2);

Int_hist_G = mean((Cy3(:,2*(start_frame+1))+Cy5(:,2*(start_frame+1))),2);

figure('Name','Green intensity after');
histogram(Int_hist_G,100, 'BinLimits', [prctile(Int_hist_G,2,'all'), prctile(Int_hist_G,98,'all')]);

Int_hist_R = mean((Cy5(:,2*start_frame+1)),2);

figure('Name','Red intensity after');
histogram(Int_hist_R,100, 'BinLimits', [prctile(Int_hist_R,2,'all'), prctile(Int_hist_R,98,'all')]);

%%
Red_thr = 1000;
Green_thr = 5000;
Red_thr_after = 300;
Green_thr_after = 5000;
FRET_threshold_min = 0.6;
FRET_threshold_max = 0.9;

%% Processing traces


for i =1: n_traces
    if rem(i,1000) == 0
        i
    end
    s = Seq(i,idx);

    x =4*(strfind(code(3,:),s(3))-1)+(strfind(code(1,:),s(1)));
    y =4*(strfind(code(4,:),s(4))-1)+(strfind(code(2,:),s(2)));
    htmap_counter(x,y) = htmap_counter(x,y)+1;

    Cy5t = Cy5(i,2:2:end);
    Cy3t = Cy3(i,2:2:end);
    
    C5dir = Cy5(i,1:2:end);
    

    Int = Cy5t + Cy3t;
    FRET = Cy5t./Int;
    Max = mean(Int(1:start_frame));
    D = mean(C5dir(1:start_frame));

    if (mean(FRET(1:start_frame))< FRET_threshold_min)||(mean(FRET(1:start_frame))> FRET_threshold_max)
        continue
    end
    
    if Max<Green_thr 
        continue
    end
    
    if ~isnan(Red_thr)
        if D<Red_thr 
            continue
        end
    end
    Max = prctile(Int(start_frame+1:end),95);
    D = prctile(C5dir(start_frame+1:end),95);
    
    if Max<Green_thr_after 
        continue
    end
    
    if ~isnan(Red_thr_after)
        if D<Red_thr_after 
            continue
        end
    end
    

    Int_ave = smooth(Int,window,'moving');
    Cy3_ave = smooth(Cy3t,window,'moving');
    Cy5_ave = smooth(Cy5t,window,'moving');

    max_int_all(i) = Max;
    [Min,I] = min(Int_ave);
    min_int_all(i) = Min;
    Length = length(Int);
    Int_ave = smooth(Int,window,'moving');
    Cy5t_ave = smooth(Cy5t,window,'moving');
    FRET_ave = Cy5t_ave./Int_ave;
    Max = mean(Int(1:10));
    max_int(i) = Max;

    idx_t = (FRET>0)&(FRET<1.2)&(Int_ave'>0.6*Max)&(C5dir>0.6*D);
    idx_t_B = zeros(1,length(idx_t));
    idx_t_A = zeros(1,length(idx_t));
    idx_t_B = logical(idx_t_B);
    idx_t_A = logical(idx_t_A);
    idx_t_B(1:start_frame-1) = idx_t(1:start_frame-1);
    idx_t_A(start_frame:end) = idx_t(start_frame:end);
    
    
    Cy3_B_t = NaN(max_length_B,1);
    Cy3_B_t(idx_t_B) = Cy3t(idx_t_B);
    Cy3_A_t = NaN(max_length_A,1);
    Cy3_A_t(idx_t_A) = Cy3t(idx_t_A);
    
    Cy5_B_t = NaN(max_length_B,1);
    Cy5_B_t(idx_t_B) = Cy5t(idx_t_B);
    Cy5_A_t = NaN(max_length_A,1);
    Cy5_A_t(idx_t_A) = Cy5t(idx_t_A);
    
    FRET_B_t = NaN(max_length_B,1);
    FRET_B_t(idx_t_B) = FRET(idx_t_B);
    FRET_A_t = NaN(max_length_A,1);
    FRET_A_t(idx_t_A) = FRET(idx_t_A);
    
    
    sel_FRET_A{x,y}(((htmap_counter(x,y)-1)*max_length_A+1):(htmap_counter(x,y)*max_length_A)) =  FRET_A_t;
    sel_FRET_B{x,y}(((htmap_counter(x,y)-1)*max_length_B+1):(htmap_counter(x,y)*max_length_B)) =  FRET_B_t;
    sel_Cy3_A{x,y}(((htmap_counter(x,y)-1)*max_length_A+1):(htmap_counter(x,y)*max_length_A)) =  Cy3_A_t;
    sel_Cy3_B{x,y}(((htmap_counter(x,y)-1)*max_length_B+1):(htmap_counter(x,y)*max_length_B)) =  Cy3_B_t;
    sel_Cy5_A{x,y}(((htmap_counter(x,y)-1)*max_length_A+1):(htmap_counter(x,y)*max_length_A)) =  Cy5_A_t;
    sel_Cy5_B{x,y}(((htmap_counter(x,y)-1)*max_length_B+1):(htmap_counter(x,y)*max_length_B)) =  Cy5_B_t;
   
end

%% Splitting traces for HMM with minimal length
minLength = 5;
for a = 16:-1:1
    for b = 16:-1:1
%         sel_FRET_B_HMM{a,b} = {};
        
        % Ensure the data starts and ends with NaN for consistent processing
        dataExtended = [NaN; sel_FRET_A{a,b}; NaN];

        % Find indices where the data transitions from or to NaN
        nanIndices = find(isnan(dataExtended));

        % Calculate the start and end indices of each stretch
        startIndices = nanIndices(1:end-1) + 1;
        endIndices = nanIndices(2:end) - 1;

        % Calculate lengths of each stretch
        lengths = endIndices - startIndices + 1;

        % Filter stretches longer than minLength
        longStretchesIdx = find(lengths >minLength);

        % Preallocate cell array for stretches longer than N points
        stretches = cell(length(longStretchesIdx), 1);

        % Extract and store each long stretch into the cell array
        for i = 1:length(longStretchesIdx)
            idx = longStretchesIdx(i);
            stretches{i} = dataExtended(startIndices(idx):endIndices(idx));
        end
        sel_FRET_A_HMM{a,b} = stretches;
    end
end
