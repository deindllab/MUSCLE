%% Listing all possible hairpin constructs

% Define the constant sequences
constSeq1 = 'TGGAGG';
constSeq2 = 'AAAAAAAAAAAAAAAAAA';
constSeq3 = 'CCTCCA';

% Define the nucleotides
nucleotides = ['A', 'C', 'G', 'T'];

% Initialize the array to store the sequences
sequences = cell(256, 1);
index = 1;

% Generate all combinations
for i = 1:length(nucleotides)
    for j = 1:length(nucleotides)
        for k = 1:length(nucleotides)
            for l = 1:length(nucleotides)
                % Construct the sequence
                seq = [constSeq1, nucleotides(i), nucleotides(j), ...
                       constSeq2, nucleotides(k), nucleotides(l), ...
                       constSeq3];
                % Store the sequence
                sequences{index} = seq;
                index = index + 1;
            end
        end
    end
end



% sequences now contains all 256 possible library sequences
%%
fileName = 'D:\MUSCLE paper\UNAfold 20240417 10C.txt';


% Read the file without using headers, assuming tab-delimiter
data = readtable(fileName, 'Delimiter', '\t', 'ReadVariableNames', false);

% Manually assign column names if the first row is not headers
data.Properties.VariableNames = {'Sequence', 'DeltaG', 'DeltaH', 'DeltaS', 'Tm', 'ImageLink', 'ThermodynamicDetails'};

% Display the first few rows to verify
% head(data);

%%

% Convert to string array if not already (for MATLAB R2016b and later)
data.DeltaG = string(data.DeltaG);

% Use regexprep to keep only the numeric part (including decimal point and minus sign)
numericStrings = regexprep(data.DeltaG, '[^\d.-]', '');

% Convert the cleaned strings to numeric values
numericValues = str2double(numericStrings);

%%
code = ['ATGC';'ATGC';'TACG';'TACG'];
idx = [7,8,27,28];
deltaG_UNAfold_htmap = NaN(16);
htmap_seq = cell(16);
for i = 1:256
    s = data{i,1};
    s = s{1}; 
    s = s(idx);
    x =4*(strfind(code(4,:),s(4))-1)+(strfind(code(3,:),s(3)));
    y =4*(strfind(code(1,:),s(1))-1)+(strfind(code(2,:),s(2)));
    htmap_seq{x,y} = s;
    deltaG_UNAfold_htmap(x,y) = numericValues(i);
end

min(deltaG_UNAfold_htmap(:))
max(deltaG_UNAfold_htmap(:))

%%
if ~exist('saveFolder','var')
    saveFolder = uigetdir('','Select a folder to save figures in');
end

MUSCLE_heatmap(deltaG_UNAfold_htmap,[-12 -6] , code_htmap, 'hmTitle', 'dG 10deg UNAfold', 'saveFolder', saveFolder);
