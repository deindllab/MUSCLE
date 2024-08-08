%% Grouping sequences based on the number of consecutive mismatches

s1 = 'GATG';
clear mismatch mismatch_seq s
mismatch{5} = [];
mismatch_seq{5} = [];
% Pooling the sequences that have 0 to 4 mismatches starting from the
% PAM-distal side (e.g. 2: XXMM, 4:XXXX, where X is a mismatch, M is
% match
code_old = 'ATGC';
code = ['GCAT';'ATGC';'TGCA';'GCAT']; % Cas9

% weird = [];
s(4) = 'A';
% s = string(s);
for a = 1:4
    s(1) = code_old(a);
    for b = 1:4
        s(2) = code_old(b);
        for c = 1:4
            s(3) = code_old(c);
            for d = 1:4
                s(4) = code_old(d);

                x = 4*(strfind(code(3,:),s(3))-1)+(strfind(code(1,:),s(1)));
                y = 4*(strfind(code(4,:),s(4))-1)+(strfind(code(2,:),s(2)));
                
                switch num2str(s == s1)
                   case '1  1  1  1'
                        k = 1;
                        mismatch{k} = [mismatch{k}; x y];
                        mismatch_seq{k} = [mismatch_seq{k}; s];
                   case '0  1  1  1'
                        k = 2;
                        mismatch{k} = [mismatch{k}; x y];
                        mismatch_seq{k} = [mismatch_seq{k}; s];
                   case '0  0  1  1'
                        k = 3;
                        mismatch{k} = [mismatch{k}; x y];
                        mismatch_seq{k} = [mismatch_seq{k}; s];
                   case '0  0  0  1'
                        k = 4;
                        mismatch{k} = [mismatch{k}; x y];
                        mismatch_seq{k} = [mismatch_seq{k}; s];
                   case '0  0  0  0'
                        k = 5;
                        mismatch{k} = [mismatch{k}; x y];
                        mismatch_seq{k} = [mismatch_seq{k}; s];
                end                                                   
            end
        end
    end
end

%% Grouping sequences based on the total number of mismatches
s1 = 'GATG';
clear mismatch mismatch_seq s
mismatch{5} = [];
mismatch_seq{5} = [];
% Pooling the sequences that have 0 to 4 mismatches in total
code_old = 'ATGC';
code = ['GCAT';'ATGC';'TGCA';'GCAT']; % Cas9

s(4) = 'A';
% s = string(s);
for a = 1:4
    s(1) = code_old(a);
    for b = 1:4
        s(2) = code_old(b);
        for c = 1:4
            s(3) = code_old(c);
            for d = 1:4
                s(4) = code_old(d);

                x = 4*(strfind(code(3,:),s(3))-1)+(strfind(code(1,:),s(1)));
                y = 4*(strfind(code(4,:),s(4))-1)+(strfind(code(2,:),s(2)));
                
                switch sum(s == s1)
                   case 4
                        k = 1;
                        mismatch{k} = [mismatch{k}; x y];
                        mismatch_seq{k} = [mismatch_seq{k}; s];
                   case 3
                        k = 2;
                        mismatch{k} = [mismatch{k}; x y];
                        mismatch_seq{k} = [mismatch_seq{k}; s];
                   case 2
                        k = 3;
                        mismatch{k} = [mismatch{k}; x y];
                        mismatch_seq{k} = [mismatch_seq{k}; s];
                   case 1
                        k = 4;
                        mismatch{k} = [mismatch{k}; x y];
                        mismatch_seq{k} = [mismatch_seq{k}; s];
                   case 0
                        k = 5;
                        mismatch{k} = [mismatch{k}; x y];
                        mismatch_seq{k} = [mismatch_seq{k}; s];
                end                                                   
            end
        end
    end
end
%% Grouping sequences based on the total number of GC pairs
s1 = 'GATG';
clear mismatch mismatch_seq s
mismatch{5} = [];
mismatch_seq{5} = [];
% Pooling the sequences that have 0 to 4 mismatches in total
code_old = 'ATGC';
code = ['GCAT';'ATGC';'TGCA';'GCAT']; % Cas9

s(4) = 'A';
% s = string(s);
for a = 1:4
    s(1) = code_old(a);
    for b = 1:4
        s(2) = code_old(b);
        for c = 1:4
            s(3) = code_old(c);
            for d = 1:4
                s(4) = code_old(d);

                x = 4*(strfind(code(3,:),s(3))-1)+(strfind(code(1,:),s(1)));
                y = 4*(strfind(code(4,:),s(4))-1)+(strfind(code(2,:),s(2)));
                 % count the number of G or C instances
                numG = count(s, 'G');
                numC = count(s, 'C');

                % total count
                GC = numG + numC;
                switch GC
                   case 4
                        k = 5;
                        mismatch{k} = [mismatch{k}; x y];
                        mismatch_seq{k} = [mismatch_seq{k}; s];
                   case 3
                        k = 4;
                        mismatch{k} = [mismatch{k}; x y];
                        mismatch_seq{k} = [mismatch_seq{k}; s];
                   case 2
                        k = 3;
                        mismatch{k} = [mismatch{k}; x y];
                        mismatch_seq{k} = [mismatch_seq{k}; s];
                   case 1
                        k = 2;
                        mismatch{k} = [mismatch{k}; x y];
                        mismatch_seq{k} = [mismatch_seq{k}; s];
                   case 0
                        k = 1;
                        mismatch{k} = [mismatch{k}; x y];
                        mismatch_seq{k} = [mismatch_seq{k}; s];
                end                                                   
            end
        end
    end
end



%% Grouping sequences based on the position of a single mismatch
s1 = 'GATG';
clear mismatch mismatch_seq s
mismatch{4} = [];
mismatch_seq{4} = [];
% Pooling the sequences that have a single mismatch in positions 1 to 4 (e.g. 2: MXMM, 4:MMMX, where X is a mismatch, M is
% match
code_old = 'ATGC';

code = ['GCAT';'ATGC';'TGCA';'GCAT']; % Cas9

s(4) = 'A';
% s = string(s);
for a = 1:4
    s(1) = code_old(a);
    for b = 1:4
        s(2) = code_old(b);
        for c = 1:4
            s(3) = code_old(c);
            for d = 1:4
                s(4) = code_old(d);

                x = 4*(strfind(code(3,:),s(3))-1)+(strfind(code(1,:),s(1)));
                y = 4*(strfind(code(4,:),s(4))-1)+(strfind(code(2,:),s(2)));
                
                switch num2str(s == s1)
                           case '0  1  1  1'
                                k = 1;
                                mismatch{k} = [mismatch{k}; x y];
                                mismatch_seq{k} = [mismatch_seq{k}; s];
                           case '1  0  1  1'
                                k = 2;
                                mismatch{k} = [mismatch{k}; x y];
                                mismatch_seq{k} = [mismatch_seq{k}; s];
                           case '1  1  0  1'
                                k = 3;
                                mismatch{k} = [mismatch{k}; x y];
                                mismatch_seq{k} = [mismatch_seq{k}; s];
                           case '1  1  1  0'
                                k = 4;
                                mismatch{k} = [mismatch{k}; x y];
                                mismatch_seq{k} = [mismatch_seq{k}; s];   
                end
            end
        end
    end
end
%% Grouping sequences based on the number of identical mismatches
s1 = 'GATG';
s1_rc = 'CTAC';
clear mismatch mismatch_seq s
mismatch{5} = [];
mismatch_seq{5} = [];
code_old = 'ATGC';

code = ['GCAT';'ATGC';'TGCA';'GCAT']; % Cas9
s(4) = 'A';
% s = string(s);
for a = 1:4
    s(1) = code_old(a);
    for b = 1:4
        s(2) = code_old(b);
        for c = 1:4
            s(3) = code_old(c);
            for d = 1:4
                s(4) = code_old(d);

                x = 4*(strfind(code(3,:),s(3))-1)+(strfind(code(1,:),s(1)));
                y = 4*(strfind(code(4,:),s(4))-1)+(strfind(code(2,:),s(2)));
                
                switch sum(s == s1)
                   case 4
                        k = 1;
                        mismatch{k} = [mismatch{k}; x y];
                        mismatch_seq{k} = [mismatch_seq{k}; s];
                   case 3
                        if sum (s == s1_rc) == 1
                            k = 2;
                            mismatch{k} = [mismatch{k}; x y];
                            mismatch_seq{k} = [mismatch_seq{k}; s];
                        end
                        
                   case 2
                        if sum (s == s1_rc) == 2
                            k = 3;
                            mismatch{k} = [mismatch{k}; x y];
                            mismatch_seq{k} = [mismatch_seq{k}; s];
                        end
                        
                   case 1
                        if sum (s == s1_rc) == 3
                            k = 4;
                            mismatch{k} = [mismatch{k}; x y];
                            mismatch_seq{k} = [mismatch_seq{k}; s];
                        end
                        
                   case 0
                        if sum (s == s1_rc) == 4
                            k = 5;
                            mismatch{k} = [mismatch{k}; x y];
                            mismatch_seq{k} = [mismatch_seq{k}; s];
                        end
                        
                end                                                   
            end
        end
    end
end

%% Grouping sequences based on the number of non-identical mismatches
s1 = 'GATG';
s1_rc = 'CTAC';
clear mismatch mismatch_seq s
mismatch{5} = [];
mismatch_seq{5} = [];
code_old = 'ATGC';

code = ['GCAT';'ATGC';'TGCA';'GCAT']; % Cas9
s(4) = 'A';
% s = string(s);
for a = 1:4
    s(1) = code_old(a);
    for b = 1:4
        s(2) = code_old(b);
        for c = 1:4
            s(3) = code_old(c);
            for d = 1:4
                s(4) = code_old(d);

                x = 4*(strfind(code(3,:),s(3))-1)+(strfind(code(1,:),s(1)));
                y = 4*(strfind(code(4,:),s(4))-1)+(strfind(code(2,:),s(2)));
                
                switch sum(s == s1)
                   case 4
                        k = 1;
                        mismatch{k} = [mismatch{k}; x y];
                        mismatch_seq{k} = [mismatch_seq{k}; s];
                   case 3
                        if sum (s == s1_rc) == 0
                            k = 2;
                            mismatch{k} = [mismatch{k}; x y];
                            mismatch_seq{k} = [mismatch_seq{k}; s];
                        end
                        
                   case 2
                        if sum (s == s1_rc) == 0
                            k = 3;
                            mismatch{k} = [mismatch{k}; x y];
                            mismatch_seq{k} = [mismatch_seq{k}; s];
                        end
                        
                   case 1
                        if sum (s == s1_rc) == 0
                            k = 4;
                            mismatch{k} = [mismatch{k}; x y];
                            mismatch_seq{k} = [mismatch_seq{k}; s];
                        end
                        
                   case 0
                        if sum (s == s1_rc) == 0
                            k = 5;
                            mismatch{k} = [mismatch{k}; x y];
                            mismatch_seq{k} = [mismatch_seq{k}; s];
                        end
                        
                end                                                   
            end
        end
    end
end

%% Example sequences

% example_seq = ['CATT';'GATC';'ACAG';'AAGG'];
example_seq = ['GATC'];
% example_seq = ['GATG';'CATG';'TTTG';'ATGG';'TCAA'];
clear mismatch
mismatch{size(example_seq,1)} = [];
for i = 1:size(example_seq,1)
    s = example_seq(i,:);
    x = 4*(strfind(code(3,:),s(3))-1)+(strfind(code(1,:),s(1)));
    y = 4*(strfind(code(4,:),s(4))-1)+(strfind(code(2,:),s(2)));
    mismatch{i} = [x y];
end
