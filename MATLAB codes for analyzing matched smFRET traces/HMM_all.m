% choose the sequence to plot
% seq_to_plot = 'TTAA';
% code = ['ATGC';'ATGC';'TACG';'TACG']; % Hairpin
code = ['GCAT';'ATGC';'TGCA';'GCAT']; % Cas9
% x = 4*(strfind(code(4,:),seq_to_plot(4))-1)+(strfind(code(3,:),seq_to_plot(3)));
% y = 4*(strfind(code(1,:),seq_to_plot(1))-1)+(strfind(code(2,:),seq_to_plot(2)));
minFRET = -0.2;
maxFRET = 1.2;
step = 0.05;
edges = [-Inf minFRET:step:maxFRET Inf]; %for discretization
trans = [0.98,0.02; 0.1,0.9]; % Initial guess of transition probabilities
sigma1 = 0.1165; % Width of emission peaks
sigma2 = 0.1185; % Width of emission peaks
center1 = 0.354; % Low FRET peak guess
center2 = 0.757; % High FRET peak guess


emis = [normpdf(edges,center1,sigma);
   normpdf(edges,center2,sigma)];
emis(emis==0) = 0.1;
emis(1,:) = emis(1,:)/sum(emis(1,:));
emis(2,:) = emis(2,:)/sum(emis(2,:));
% x = 9;
% y = 9;
clear HMM_estTR HMM_estE
HMM_estTR{16,16} = [];
HMM_estE{16,16} = [];
for i = 1:16
    for j = 1:16
        [i j]
        data = sel_FRET_A_HMM{i,j};

        clear data_discr;
        n_traces = length(data);
        if n_traces>0
            data_discr{n_traces} = [];
            % data_concat = [];
            % trace_edges = [];
            for k = 1:n_traces
%                data_discr{k} = transpose(discretize(data{k},edges));
                data_discr{k} = transpose(discretize(medfilt1(data{k},2),edges));
            %    data_concat = [data_concat data{i}];
            %    trace_edges = [trace_edges length(data_concat)];
            end
            % trace_edges = trace_edges(1:end-2);

            % data_concat = seq2;

            [estTR,estE] = hmmtrain_fixEm(data_discr,trans,emis);
            HMM_estTR{i,j} = estTR;
            HMM_estE{i,j} = estE;
            % [estTR,estE] = hmmtrain(data_concat,trans,emis,'Verbose',true);
        else
            HMM_estTR{i,j} = NaN(2,2);
            HMM_estE{i,j} = NaN(size(emis,1),size(emis,2));        
        end
    end
end
%% Kinetic constant heatmaps
frame_rate = 2.5;
HMM_kopen = NaN(16,16);
HMM_kclose = NaN(16,16);
for i = 1:16
    for j = 1:16
        temp = HMM_estTR{i,j}*frame_rate; % to convert to s-1
        HMM_kopen(i,j) = temp(2,1);
        HMM_kclose(i,j) = temp(1,2);
    end
end
if ~exist('saveFolder','var')
    saveFolder = uigetdir('','Select a folder to save figures in');
end
MUSCLE_heatmap(HMM_kclose,[0 0.7], code, 'hmTitle', 'krewind HMM', 'boxLabelsFontSize', 7, 'saveFolder', saveFolder);
MUSCLE_heatmap(HMM_kopen,[0 0.2], code, 'hmTitle', 'kunwind HMM', 'boxLabelsFontSize', 7, 'saveFolder', saveFolder);
% equilibrium constant heatmap
HMM_Kequil = HMM_kopen./HMM_kclose;
MUSCLE_heatmap(HMM_Kequil,[0 4], code, 'hmTitle', 'kunwind/krewind HMM', 'boxLabelsFontSize', 7, 'saveFolder', saveFolder);
R = 1.987e-3; % Universal gas constant for dG calculations in kcal*K-1*mol-1
T = 293; %20deg Celsius in K
HMM_dG = -R*T*log(HMM_Kequil);
MUSCLE_heatmap(HMM_dG,[-1 3] , code, 'hmTitle', 'dG HMM', 'boxLabelsFontSize', 7, 'saveFolder', saveFolder);
