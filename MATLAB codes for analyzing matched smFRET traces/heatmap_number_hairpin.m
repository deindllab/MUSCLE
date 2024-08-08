
% code = 'ATGC';
% code = ['TGCA';'GCAT';'ATGC';'CATG'];
code = ['ATGC';'ATGC';'TACG';'TACG'];
x_axis = 1:16;
y_axis = 1:16;
htmap_num_traces = zeros(16,16);
idx = [39,40,59,60];
% idx = [1:6];
Seq1 = Seq(:,idx);

% s1  = 'TATTTA';
% weird = [];
% Complementary sequences should be at the diagonal
for i = 1: length(Seq1)
   s = Seq1(i,:);
   x =4*(strfind(code(4,:),s(4))-1)+(strfind(code(3,:),s(3)));
   y =4*(strfind(code(1,:),s(1))-1)+(strfind(code(2,:),s(2)));
   htmap_num_traces(x,y) = htmap_num_traces(x,y) + 1;
%    if (x == 18) & (y == 21)
%        weird = [weird; Seq(i,:)];
%    end
end

%% create code for hairpin heatmaps
code_htmap = code;
code_htmap(1,:) = code(3,:);
code_htmap(3,:) = code(4,:);
code_htmap(4,:) = code(1,:);

%% Plot FRET heatmap Before
if ~exist('saveFolder','var')
    saveFolder = uigetdir('','Select a folder to save figures in');
end

MUSCLE_heatmap(htmap_num_traces,[80 1000] , code, 'hmTitle', 'Number of traces', 'boxLabelsFontSize', 7, 'saveFolder', saveFolder, 'cmap', flipud(autumn));
MUSCLE_heatmap(htmap_num_traces,[80 1000] , code_htmap, 'hmTitle', 'Number of traces', 'saveFolder', saveFolder, 'cmap', flipud(autumn));