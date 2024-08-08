% 1. Start by importing the Bioinformatics Toolbox.
import bioinformatics.*

% 2. Read the FASTQ file using the fastqread function.
% fastqStruct = fastqread('D:\20231110_Misha_+4_20sec_purified\20231110_20sec_+4_purified.fastq');
% fastqStruct = fastqread('D:\240318_Misha_GEMIN5\0s_test1_S1_L001_R1_001.fastq');
fastqStruct = fastqread('D:\240318_Misha_GEMIN5\20s_test1_S1_L001_R1_001.fastq');

% fastqStruct = fastqread('D:\20230921_FC_Full_Ha4N\test1_S1_L001_R1_001.fastq','blockread', [5 10]);
% 3. Initialize a cell array to store the sequences.
seqCellArray = cell(length(fastqStruct), 1);

% 4. Loop over the fastqStruct and extract the Sequence field, which contains the sequences.
for i = 1:length(fastqStruct)
    seqCellArray{i} = fastqStruct(i).Sequence;
end

% 5. Convert the cell array of sequences into a char array. This will only work if all the sequences are the same length.
% If the sequences are not the same length, they cannot be directly stored in a 2D char array.
seqCharArray = char(seqCellArray);
seqCharArray = seqCharArray(:,1:62);
clear fastqStruct seqCellArray

%% Select sequences
standard = 'GGTCTCGTCCAATCTAGGTCATCAATTATTATACATCGGAGCCCTGCCAATCTCGTATGCCG'; % Normalization standard CCR5
% standard = 'GGTCTCGCACAGCAGAAATCTCTGCTAGAATATAAAGATGAGACGCTGGAGTACAAACGTCA'; % Standard Cas9 library
library = 'GGTCTCGCACAGCAGAAATCTCTGCTGCTGTGAGGCTACGAGCGGTTGGAGTACAAACGTCA'; % GEMIN5

% Counting the number of CCR5 reads that serves as a normalization standard 
mask = bsxfun(@eq, seqCharArray, standard); % Create a logical matrix where each element is 1 if it matches the standard, 0 otherwise
rowSums = sum(mask, 2); % Sum along the rows
idx = find(rowSums > 40); % Find the indices of rows where the sum is greater than 40
% idx = find(rowSums > 20); % Find the indices of rows where the sum is greater than 40
CCR5_20s = sum(idx(:));

mask = bsxfun(@eq, seqCharArray, library); % Create a logical matrix where each element is 1 if it matches the standard, 0 otherwise
rowSums = sum(mask, 2); % Sum along the rows
idx = find(rowSums > 40); % Find the indices of rows where the sum is greater than 40
% idx = find(rowSums > 20); % Find the indices of rows where the sum is greater than 40
Seq = seqCharArray(idx,:);

%% Quantifying the number of traces for each library member

code = ['GCAT';'ATGC';'TGCA';'GCAT']; % GATG target

x_axis = 1:16;
y_axis = 1:16;
htmap_num_traces = zeros(16,16);
idx = 27:30; 
Seq1 = Seq(:,idx);


for i = 1: length(Seq1)
   s = Seq1(i,:);
   x = 4*(strfind(code(3,:),s(3))-1)+(strfind(code(1,:),s(1)));
   y = 4*(strfind(code(4,:),s(4))-1)+(strfind(code(2,:),s(2)));
   htmap_num_traces(x,y) = htmap_num_traces(x,y) + 1;
%    if (x == 18) & (y == 21)
%        weird = [weird; Seq(i,:)];
%    end
end

%% Normalizing read counts

htmap_20s = htmap_num_traces;
htmap_20s_norm = htmap_20s/CCR5_20s;

%% Repeating the same steps for 0 s
% 1. Start by importing the Bioinformatics Toolbox.
import bioinformatics.*

% 2. Read the FASTQ file using the fastqread function.
% fastqStruct = fastqread('D:\20231110_Misha_+4_20sec_purified\20231110_20sec_+4_purified.fastq');
% fastqStruct = fastqread('D:\240318_Misha_GEMIN5\0s_test1_S1_L001_R1_001.fastq');
fastqStruct = fastqread('D:\240318_Misha_GEMIN5\20s_test1_S1_L001_R1_001.fastq');

% fastqStruct = fastqread('D:\20230921_FC_Full_Ha4N\test1_S1_L001_R1_001.fastq','blockread', [5 10]);
% 3. Initialize a cell array to store the sequences.
seqCellArray = cell(length(fastqStruct), 1);

% 4. Loop over the fastqStruct and extract the Sequence field, which contains the sequences.
for i = 1:length(fastqStruct)
    seqCellArray{i} = fastqStruct(i).Sequence;
end

% 5. Convert the cell array of sequences into a char array. This will only work if all the sequences are the same length.
% If the sequences are not the same length, they cannot be directly stored in a 2D char array.
seqCharArray = char(seqCellArray);
seqCharArray = seqCharArray(:,1:62);
clear fastqStruct seqCellArray

%% Select sequences
standard = 'GGTCTCGTCCAATCTAGGTCATCAATTATTATACATCGGAGCCCTGCCAATCTCGTATGCCG'; % Normalization standard CCR5
% standard = 'GGTCTCGCACAGCAGAAATCTCTGCTAGAATATAAAGATGAGACGCTGGAGTACAAACGTCA'; % Standard Cas9 library
library = 'GGTCTCGCACAGCAGAAATCTCTGCTGCTGTGAGGCTACGAGCGGTTGGAGTACAAACGTCA'; % GEMIN5

% Counting the number of CCR5 reads that serves as a normalization standard 
mask = bsxfun(@eq, seqCharArray, standard); % Create a logical matrix where each element is 1 if it matches the standard, 0 otherwise
rowSums = sum(mask, 2); % Sum along the rows
idx = find(rowSums > 40); % Find the indices of rows where the sum is greater than 40
% idx = find(rowSums > 20); % Find the indices of rows where the sum is greater than 40
CCR5_0s = sum(idx(:));

mask = bsxfun(@eq, seqCharArray, library); % Create a logical matrix where each element is 1 if it matches the standard, 0 otherwise
rowSums = sum(mask, 2); % Sum along the rows
idx = find(rowSums > 40); % Find the indices of rows where the sum is greater than 40
% idx = find(rowSums > 20); % Find the indices of rows where the sum is greater than 40
Seq = seqCharArray(idx,:);

%% Quantifying the number of traces for each library member

code = ['GCAT';'ATGC';'TGCA';'GCAT']; % GATG target

x_axis = 1:16;
y_axis = 1:16;
htmap_num_traces = zeros(16,16);
idx = 27:30; 
Seq1 = Seq(:,idx);


for i = 1: length(Seq1)
   s = Seq1(i,:);
   x = 4*(strfind(code(3,:),s(3))-1)+(strfind(code(1,:),s(1)));
   y = 4*(strfind(code(4,:),s(4))-1)+(strfind(code(2,:),s(2)));
   htmap_num_traces(x,y) = htmap_num_traces(x,y) + 1;
%    if (x == 18) & (y == 21)
%        weird = [weird; Seq(i,:)];
%    end
end

%% Normalizing read counts

htmap_0s = htmap_num_traces;
htmap_0s_norm = htmap_0s/CCR5_0s;

%% Final normalization to the uncleaved sample 0s

htmap_20s_norm2 = htmap_20s_norm./htmap_0s_norm;

htmap_GEMIN = htmap_20s_norm2;
flattened_GEMIN = htmap_GEMIN(:);

%% Plotting a heatmap
if ~exist('saveFolder','var')
    saveFolder = uigetdir('','Select a folder to save figures in');
end
MUSCLE_heatmap(htmap_GEMIN,[0.15 0.6] , code, 'hmTitle', 'Number of traces', 'saveFolder', saveFolder, 'cmap', flipud(autumn));