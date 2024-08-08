% Sebastian Deindl 2012-2014
% Goes through existing .mat file and organizes all the indicated variables into
% a structure called "U"
%
%clear all
% %warning off MATLAB:divideByZero
%
%#######################################################
% Modifications made by Joakim Laksman May 2016 for 3 color FRET
% Takes all the *.mat files in one folder and combines them into
% a  *_struc.mat file which is a structure.
%#######################################################

warning off all
clear
fileName = 'Combined.mat'

path1 = uigetdir('D:\Experiments\MUSCLE DONUTS\','Read normalized Files');
files = fuf([path1 '\*.mat'],'detail');
Cy3_all = [];
Cy5_all = [];
Seq_all = [];
Seq_2_all = [];
Dist_all = [];
x_all = [];
y_all = [];
x_FQ_all = [];
y_FQ_all = [];
minlength = 140;
ID_all = strings(0);
FOV_all = strings(0);
click_points_all = cell(0);
exp_name_all = strings(0);
ALEX_all = strings(0);

ALEX = strings(0);
click_points = cell(0);
ID = strings(0);
FOV = strings(0);
exp_name = strings(0);
Dist = [];
x = [];
y = [];
x_FQ = [];
y_FQ = [];

prompt      = {'Experiment name:','ALEX [1: GREEN. & 2: RED.]:','Override frame rate:'};
dlg_title   = 'Parameters';
num_lines   = 1;
def         = {'','',''};
input_ans   = inputdlg(prompt, dlg_title, num_lines, def);
exp_prefix     = string(input_ans{1});
ALEX_user     = string(input_ans{2});
frame_rate = str2num(input_ans{3});

if ~isempty(files)
    for j = 1:size(files,1)
        FID = fopen(files{j});
        load(files{j});
        [~,name,~] = fileparts(files{j});
        name = erase(name,'_traces');
        n_traces_file = size(Cy3,1);
        if  ~isempty(Dist)
            Dist_all = [Dist_all Dist];
            Dist = [];
        end
%         clear x y
%         x = [];
%         y = [];
        if  isempty(x)
            x = NaN(1,n_traces_file);
            y = NaN(1,n_traces_file);        
        end
        x_all = [x_all x];
        y_all = [y_all y];
        x = [];
        y = [];
        
        if  isempty(click_points)
            click_points = cell(1, n_traces_file);
        end
        
        click_points_all = [click_points_all click_points];
        click_points = cell(0);
        
        if  isempty(x_FQ)            
            x_FQ = NaN(1,n_traces_file);
            y_FQ = NaN(1,n_traces_file);
        end
        x_FQ_all = [x_FQ_all x_FQ];
        y_FQ_all = [y_FQ_all y_FQ];
        x_FQ = [];
        y_FQ = [];

        if  isempty(ID)
            ID = strings(n_traces_file,1);
            for i = 1:n_traces_file
                [~,ID(i)] = fileparts(tempname);
            end
        end
        ID_all = [ID_all; ID];
        ID = strings(0);
        
        if  isempty(ALEX)
            ALEX = ALEX_user + strings(n_traces_file,1);
        end
        ALEX_all = [ALEX_all; ALEX];
        ALEX = strings(0);
        
        if isempty(exp_name)
            exp_name_all = [exp_name_all; exp_prefix+strings(n_traces_file,1)];
        else
            exp_name_all = [exp_name_all; exp_prefix+"_"+exp_name];
            exp_name = strings(0);
        end
        
        if isempty(FOV)
            FOV_all = [FOV_all; name+strings(n_traces_file,1)];
        else
            FOV_all = [FOV_all; name+"_"+FOV];
            FOV = strings(0);
        end
        
%         Seq = Seq(:,1:49); % temporary fix for a weird bug
        Seq = Seq(:,1:62);
        Seq_2 = Seq_2(:,1:62);
        n_traces = size(Cy3,1);
        length = min([minlength size(Cy3,2)]);
        Z = zeros(n_traces, minlength);
        Cy3t = Z;
        Cy3t(:,1:length) = Cy3(:,1:length);
        Cy5t = Z;
        Cy5t(:,1:length) = Cy5(:,1:length);
        
        
        Cy3_all = [Cy3_all; Cy3t];
        Cy5_all = [Cy5_all; Cy5t];
        Seq_all = [Seq_all; Seq];
        Seq_2_all = [Seq_2_all; Seq_2];
    end
    
    Cy3 = Cy3_all;
    Cy5 = Cy5_all;
    Seq = Seq_all;
    Seq_2 = Seq_2_all;
    click_points = click_points_all;
    Dist = Dist_all;
    ALEX = ALEX_all;
    ID = ID_all;
    exp_name = exp_name_all;
    FOV = FOV_all;
    n_traces = size(Cy3,1);
    if isempty(Dist)
        Dist = NaN(n_traces,1);
    end
    x = x_all;
    y = y_all;
    x_FQ = x_FQ_all;
    y_FQ = y_FQ_all;
    if ~isempty(frame_rate)
        time = 0:(size(Cy3,2)-1);
        time = time/frame_rate;
    end

    if ~isempty(exp_prefix)
        fileName = exp_prefix+'.mat';
    end
    
    save(fullfile(path1,fileName),'Cy3', 'Cy5', 'Seq','Seq_2','time','Dist','x','y','x_FQ','y_FQ','exp_name','ALEX','FOV','ID','click_points');
else
    disp('No files found');
end
fclose('all');
clearvars -except Cy3 Cy5 Seq Seq_2 time Dist x y x_FQ y_FQ exp_name ALEX FOV ID click_points
