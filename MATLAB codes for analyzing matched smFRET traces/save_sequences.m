%% Sequence list
pathname = uigetdir('','Select a folder to save traces in');
% sequence_list = ['CCCG';'AAGG';'CGGG';'AGTG';'GGTG';'GGGT';'CCCC';'GCCC';'CCGA';'GATC';'GCGG';'CCCT';'ATGA';'ATGT';'ATGC';'CCGG';'CGAC';'GCCG';'GGGC';'TTTG';'TGTG';'GATG';'AATG';'CATT';'TTTA';'TTGC'];
% sequence_list = ['AGTT';'GCCT';'GTGA'];
% sequence_list = ['AACC';'CTTG'];
% sequence_list = ['GATG';'AATG';'CTGG';'TTTA'];
% sequence_list = ['AAAA';'CAAG';'CCGG'];
% sequence_list = ['GAAC'];
% sequence_list = ['GAAC'];
% sequence_list = ['GATG';'CATG';'TTTG';'ATGG';'TCAA'];
sequence_list = ['CATT';'GATC';'ACAG';'AAGG'];
% sequence_list = ['TGGG';'AGGC';'CCGG';'GGCG';'TACC';'CGAC';'AAGG';'ACAG';'ACTC';'CAGC';'GTGC';'AGCC';'CGGG';'AGCT';'GCAC'];
nSeq = size(sequence_list,1);
index = 27:30; % Ha
% index = [39,40,59,60];  % Hairpin
Seq1 = Seq(:,index);
for i = 1:nSeq
    sequence = sequence_list(i,:);
    filename = ['Traces_', sequence, '.mat'];
    idx = sum(Seq1 == sequence,2) == 4;
    % Apply index
    clear S
    S.Dist = Dist(idx);
    S.Cy3 = Cy3(idx,:);
    S.Cy5 = Cy5(idx,:);
    % FRET = FRET(idx);
    S.Seq = Seq(idx,:);
    if exist('click_points','var') == 1
        S.click_points = click_points(idx);
    end
    if exist('ALEX','var') == 1
        S.ALEX = ALEX(idx);
    end
    if exist('ID','var') == 1
        S.ID = ID(idx);
    end
    if exist('exp_name','var') == 1
        S.exp_name = exp_name(idx);
    end
%     if exist('x','var') == 1
%         S.x = x(idx);
%         S.y = y(idx);
%     end
    if exist('x_FQ','var') == 1
        S.x_FQ = x_FQ(idx);
        S.y_FQ = y_FQ(idx);
    end
    if exist('FOV','var') == 1
        S.FOV = FOV(idx);
    end
    % Save


    % If the user did not cancel the dialog
    save(fullfile(pathname, filename), '-struct', 'S');
end
    