%% Select sequences
standard = '  GGTCTCGTCCAATCTAGGTCATCAATTATTATACATCGGAGCCCTGCCAATCTCGTATGCCG'; % CCR5 control
% standard =
% 'GGTCTCGCACAGCAGAAATCTCTGCTAGAATATAAAGATGAGACGCTGGAGTACAAACGTCA'; % Cas9

idx = [];
for i = 1: length(Seq)
    if (sum(Seq(i,:) == standard)>40)
        idx = [idx i];
    end
end

%% Apply index
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
% if exist('x','var') == 1
%     S.x = x(idx);
%     S.y = y(idx);
% end
% if exist('x_FQ','var') == 1
%     S.x_FQ = x_FQ(idx);
%     S.y_FQ = y_FQ(idx);
% end
if exist('FOV','var') == 1
    S.FOV = FOV(idx);
end
if exist('Seq_2','var') == 1
    S.Seq_2 = Seq_2(idx,:);
end
% clearvars -except Cy3 Cy5 Seq time Dist x y x_FQ y_FQ exp_name ALEX FOV ID click_points
%% Save
% Open a save file dialog
[filename, pathname] = uiputfile('*.mat', 'Save as');

% If the user did not cancel the dialog
if ischar(filename)
    % Save the structure to the chosen file
    save(fullfile(pathname, filename), '-struct', 'S');
end
