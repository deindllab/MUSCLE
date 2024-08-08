% window = 1;
threshold = 0.02; % Threshold for Cy5 direct relative to total intensity with green excitation
clear sel_FRET sel_FRET_A sel_FRET_B sel_Cy3 sel_Cy5
% code = 'ATGC';
% code = ['TGCA';'GCAT';'ATGC';'CATG']; % CCR5
% code = ['GCAT';'ATGC';'TGCA';'GCAT']; % Ha
code = 'ATGC';
idx = 39:41;
idx_2 = 2:4;
% idx = 17:20; % CCR5
% idx = 27:30; % Ha
% FRET_threshold = 0.5;
% Seq1 = Seq(:,idx);
% Cy3_accepted = [];
% Cy5_accepted = [];
% Seq_accepted = [];
% FRET_accepted = [];
% FRET_all = zeros(size(Cy3,2)/2);
reject_norm = [];
n_traces = length(Seq);
max_int = zeros(n_traces,1);
max_int_all = zeros(n_traces,1);
min_int_all = zeros(n_traces,1);
% sel_FRET_A{16,16} = [];
% sel_Cy3_A{16,16} = [];
% sel_Cy5_A{16,16} = [];
sel_FRET_B{64,64} = [];
sel_Cy3_B{64,64} = [];
sel_Cy5_B{64,64} = [];
% sel_delta_FRET{16,16} = [];
%start_frame = 41;
% start_frame = 6;
window = 6;
n_traces = length(Seq);

%% Plotting intensity histograms
Int_hist_G = mean((Cy3(:,2:2:window)+Cy5(:,2:2:window)),2);

figure('Name','Green intensity');
histogram(Int_hist_G);

Int_hist_R = mean((Cy5(:,1:2:window)),2);

figure('Name','Red intensity');
histogram(Int_hist_R);

%%
Red_thr = 2500;
Green_thr = 5000;
%%

for i =1: n_traces
    if rem(i,5000) == 0
        i
    end
    s = Seq(i,idx);
    s_2 = Seq_2(i,idx_2);
    invalidChars = [~ismember(s, code) ~ismember(s_2, code)];
    if any(invalidChars)
        continue
    end
    
%     x = 16*(strfind(code,s(5))-1)+4*(strfind(code,s(3))-1)+(strfind(code,s(1)));
%     y = 16*(strfind(code,s(6))-1)+4*(strfind(code,s(4))-1)+(strfind(code,s(2)));Cy5t = Cy5(i,:);
%     x =4*(strfind(code(3,:),s(3))-1)+(strfind(code(1,:),s(1)));
%     y =4*(strfind(code(4,:),s(4))-1)+(strfind(code(2,:),s(2)));
    x =16*(strfind(code,s(1))-1)+4*(strfind(code,s(2))-1)+(strfind(code,s(3)));
    y =16*(strfind(code,s_2(1))-1)+4*(strfind(code,s_2(2))-1)+(strfind(code,s_2(3)));
    
    
%     y =4*(strfind(code(1,:),s(1))-1)+(strfind(code(2,:),s(2)));
    Cy5t = Cy5(i,2:2:end);
    Cy3t = Cy3(i,2:2:end);
    
    C5dir = Cy5(i,1:2:end);
    

    Int = Cy5t + Cy3t;
    FRET = Cy5t./Int;
    Max = mean(Int(1:window/2));
    D = mean(C5dir(1:window/2));
%     if D<0.01*M
%         continue % the molecule probably lacks Cy5 altogether
%     end  
%     if mean(FRET(5:start_frame-1))< FRET_threshold
%         continue
%     end
    
    if Max<Green_thr % Used to be 10000, changed for the 3-round experiment
        continue
    end
    
    if D<Red_thr % Used to be 10000, changed for the 3-round experiment
        continue
    end
    
% Estimating the maximum intensity as average intensity in the first
% 10 frames

%     Int_ave = smooth(Int,window,'moving');
%     Cy3_ave = smooth(Cy3t,window,'moving');
%     Cy5_ave = smooth(Cy5t,window,'moving');
%         Max = quantile(Int_ave,0.95);
       
%     Is the trace normalizable? Intensity below 20% of maximum for at
%     least 20 frames. If so, normalize the trace. If not - move on.
    max_int_all(i) = Max;
    [Min,I] = min(Int);
    min_int_all(i) = Min;
    Length = length(Int);
%     FRET_all(i,:) = FRET;
%         Cy3_accepted = [Cy3_accepted; Cy3t];
%         Cy5_accepted = [Cy5_accepted; Cy5t];
%         Seq_accepted = [Seq_accepted; Seq(i,:)];
%         FRET_accepted = [FRET_accepted; FRET];
    idx_t = (FRET>0)&(FRET<1.2)&(Int>0.6*Max)&(C5dir>0.6*D);
%     idx_t_B = zeros(1,length(idx_t));
%     idx_t_A = zeros(1,length(idx_t));
%     idx_t_B = logical(idx_t_B);
%     idx_t_A = logical(idx_t_A);
%     idx_t_B(1:start_frame-1) = idx_t(1:start_frame-1);
%     idx_t_A(start_frame:end) = idx_t(start_frame:end);
%     sel_FRET_A{x,y} =  [sel_FRET_A{x,y} FRET(idx_t_A)];
%     sel_Cy3_A{x,y} =  [sel_Cy3_A{x,y} Cy3t(idx_t_A)];
%     sel_Cy5_A{x,y} =  [sel_Cy5_A{x,y} Cy5t(idx_t_A)];
    sel_FRET_B{x,y} =  [sel_FRET_B{x,y} FRET(idx_t)];
    sel_Cy3_B{x,y} =  [sel_Cy3_B{x,y} Cy3t(idx_t)];
    sel_Cy5_B{x,y} =  [sel_Cy5_B{x,y} Cy5t(idx_t)];
%     sel_delta_FRET{x,y} = [sel_delta_FRET{x,y} mean(FRET(idx_t_A))-mean(FRET(idx_t_B))];

end