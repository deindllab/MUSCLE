%% Matching pattern
% 1: stem--
% 2: stem-+
% 3: stem+-
% 4: stem++
code = 'ATGC';
code_comp = 'TACG';
% code = ['ATGC';'ATGC';'TACG';'TACG'];
clear mismatch
mismatch{4} = [];

for i = 1:4
    for j = 1:4
%         seq_tmp = [code(i) code(j) code_comp(j) code_comp(i)];
%         seq_to_save = [seq_to_save; seq_tmp];
        x =4*(i-1)+j;
        y =4*(i-1)+j;
        mismatch{4} =[mismatch{4}; x y];
%         FRET_8bp_15deg = [FRET_8bp_15deg sel_FRET_A1{x,y}];
%         FRET_8bp_20deg = [FRET_8bp_20deg sel_FRET_A2{x,y}];
    end
end

for i = 1:4
    for j = 1:4
        for k = 1:4
                if (k~=j)
%                     seq_tmp = [code(i) code(j) code_comp(k) code_comp(i)];
                    x =4*(i-1)+k;
                    y =4*(i-1)+j;
                    mismatch{3} =[mismatch{3}; x y];
%                     seq_tmp = [code(j) code(i) code_comp(i) code_comp(k)];
                    x =4*(k-1)+i;
                    y =4*(j-1)+i;
                    mismatch{2} =[mismatch{2}; x y];
                end
            
        end
    end
end

for i = 1:4
    for j = 1:4
        for k = 1:4
            for l = 1:4
                if (k~=j)&&(l~=i)
%                     seq_tmp = [code(i) code(j) code_comp(k) code_comp(l)];
%                     seq_to_save = [seq_to_save; seq_tmp];    
                    x =4*(l-1)+k;
                    y =4*(i-1)+j;
                    mismatch{1} =[mismatch{1}; x y];
%                     FRET_6bp_15deg = [FRET_6bp_15deg sel_FRET_A1{x,y}];
%                     FRET_6bp_20deg = [FRET_6bp_20deg sel_FRET_A2{x,y}];
                end
            end
        end
    end
end

%% Example sequences
code_1 = ['ATGC';'ATGC';'TACG';'TACG'];


% example_seq = ['CATT';'GATC';'ACAG';'AAGG'];
example_seq = ['AAAA';'TTAA'];
% example_seq = ['GATG';'CATG';'TTTG';'ATGG';'TCAA'];
clear mismatch
mismatch{size(example_seq,1)} = [];
for i = 1:size(example_seq,1)
    s = example_seq(i,:);
    x =4*(strfind(code_1(4,:),s(4))-1)+(strfind(code_1(3,:),s(3)));
    y =4*(strfind(code_1(1,:),s(1))-1)+(strfind(code_1(2,:),s(2)));
    mismatch{i} = [x y];
end