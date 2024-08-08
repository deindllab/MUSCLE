% A version of the readtraces code to go through the MUSCLE data and
% manually curate and optionally normalize traces, while keeping track of
% the sequences
% Anton Jan 2023

% Modified by Joakim Laksman October 2016.


function readtraces_keyboard_seq_CCR5_ALEX



    clear
%     idx1 = 17:20; % for CCR5
    idx1 = 27:30; % for Ha
%     idx1 = [39,40,59,60];
    save_accepted_only = true; % Save all traces or only the accepted ones. 
    % Use false to generate training sets for neural networks
    B=2;

    C=2;
    % A=1;
    Cy3 = [];
    Cy5 = [];
    time = [];
    Seq = [];
    Seq_2 = [];
    Dist = [];
    idx = [];
    ALEX = strings(0);

    ID = strings(0);
    FOV = strings(0);
    exp_name = strings(0);
    x = [];
    y = [];
    x_FQ = [];
    y_FQ = [];

    % sort_wrap = 0;
    % n_thr = 50;
    % Cy5_thr = 10000;

    %function for processing keyboard input when choosing what to do with a
    %molecule
    function key_pressed_fcn(~,eventDat)
            switch eventDat.Key
                case 'f' 
                    B = 1;
                    delete(gcf);
                case 'd' 
                    B = 2;
                    delete(gcf);
                case 's' 
                    B = 3;
                    delete(gcf);
                case 'a' 
                    B = 4;
                    delete(gcf);
                case 'r'
                    B = 6;
                    delete(gcf);
            end

    end
    %function for processing keyboard input when choosing whether to accept a
    %trace
    function key_pressed_fcn2(~,eventDat)
            switch eventDat.Key
                case 'f' 
                    C = 1;
                    delete(gcf);
                case 'd' 
                    C = 2;
                    delete(gcf);

            end

    end

    % function key_pressed_fcn3(~,eventDat)
    %         switch eventDat.Key
    %             case 'f' 
    %                 A = 2;
    %                 delete(gcf);
    %             case 'd' 
    %                 A = 1;
    %                 delete(gcf);
    %             case 's' 
    %                 A = 2;
    %                 delete(gcf);
    %             case 'a' 
    %                 A = 4;
    %                 delete(gcf);
    %             
    %         end
    %            
    % end
        %1 = Analyze molecule
        %2 = Next molecule
        %3 = Previous molecule
        %4 = Accept molecule
        %functions for button clicks when choosing what to do with a molecule
    function callback1(~,~)
        B=1;
        delete(gcf);
    end

    function callback2(~,~)
        B=2;
        delete(gcf);
    end

    function callback3(~,~)
        B=3;
        delete(gcf);
    end

    function callback4(~,~)
        B=4;
        delete(gcf);
    end

    function callback5(~,~)
        B=5;
        delete(gcf);
    end

    function callback6(~,~)
        B=6;
        delete(gcf);
    end
    %funcitons for button clicks when choosing wether to save a trace
    function callback1_1(~,~)
        C=1;
        delete(gcf);
    end

    function callback1_2(~,~)
        C=2;
        delete(gcf);
    end

%     function callback2_1(~,~)
%         A=2;
%         delete(gcf);
%     end
%     
%     function callback2_2(fig_obj,eventDat)
%         A=1;
%         delete(gcf);
%     end
%     function callback2_3(fig_obj,eventDat)
%         A=2;
%         delete(gcf);
%     end
%     function callback2_4(fig_obj,eventDat)
%         A=4;
%         delete(gcf);
%     end
    function scale = make_scale(values)
        
        scale_min = prctile(values,2,'all');
        scale_max = prctile(values,98,'all');
        scale_span = scale_max - scale_min;
        scale_min = scale_min - 0.2*scale_span;
        scale_max = scale_max + 0.2*scale_span;
        scale = [scale_min scale_max];
    end
   
% fretthreshold = 10000; % this is a threshold value to prevent crazy fret traces
%%%%%%%%% First read in file from directory - either .traces or .mat
%WD = cd;

fileID = fopen('lastFolder.txt');
A = fread(fileID,'*char')';
fclose(fileID);
% disp(A);

%[fileName,path1] = uigetfile({'*.traces;*.mat'},'Read Traces File');
[fileName,path1] = uigetfile({[A '\*.mat']},'Read Traces File');
% fileName_temp=fileName;
% fileName_temp=fileName_temp(1:length(fileName_temp)-4);

%if the user changed the directory, save it here...
lastFolderFile = fopen('lastFolder.txt','w');
fprintf(lastFolderFile,'%s',path1);
fclose(lastFolderFile);
% disp(path1);
% disp(length(path));
addpath(path1);  %add path to path list each time to ensure proper file access

fopen(fileName);
load(fileName);
% fretthreshold = 100;
% FRETmin = -2;
% FRETmax = 2;
n_peaks = size(Cy3,1);
% n_peaks = n_peaks(1);
length_1 = size(Cy5,2);
idx_G = 2:2:length_1;
idx_R = 1:2:length_1;
% idx_G = 1:2:length_1; %because for the experiment with sucrose the order in ALEX was different
% idx_R = 2:2:length_1;
frame_rate = 5;
time = (0:(size(Cy3,2)-1))/frame_rate;
time_G = time(idx_G);
time_R = time(idx_R);
% frame_rate = 1/(time(2)-time(1));

if exist('click_points','var') == 0
    click_points = cell(1, n_peaks);
end

h = figure('Position', [200 100 600 500]);

%this while loop will continue until last trace in file is viewed
first = 1;
index = [];
% i = 1; % Change accordingly to start from a different molecule

for i = n_peaks:-1:1
    if ~isempty(click_points{i})
        break;
    end
end

flow_frame = 1;
flow_time = time(flow_frame);
click_lines = cell(0);
temp_save_path = strcat(path1,fileName(1:length(fileName)-4),'_weeded.mat');
while i < n_peaks+1
    if rem(i,30) == 0
        save(temp_save_path,'Cy3', 'Cy5', 'Seq','time','Dist','x','y','x_FQ','y_FQ','exp_name','ALEX','FOV','ID','click_points');
    end
    if not(isempty(Seq))
        set(h,'Name',strcat('Molecule #',num2str(i),'/', num2str(n_peaks),':::  ', string(Seq(i,idx1)), ' (GATG) Picked molecules #',num2str(size(index,2))));
    else
        set(h,'Name',strcat('Molecule #',num2str(i),'/', num2str(n_peaks),' Picked molecules #',num2str(size(index,2))));
    end
   
%     if vsel_1 ~= 0
        if first
            first = 0;
            axis1 = subplot(3,1,1);
            Cy5_plot = plot(time_G, Cy5(i,idx_G),'r');
            hold on
            xlabel('Time, seconds');
            ylabel('Intensity');
            title('Green excitation');
            grid on
            vline(flow_time);
            
            clicks = click_points{i};
            if ~isempty(clicks)
                nclicks = length(clicks);
                click_lines = cell(1,nclicks);
                for j = 1:nclicks
                    click_lines{j} = vline(time(clicks(j)),'g-');
                end
            end
            Cy3_plot = plot(time_G, Cy3(i,idx_G),'g');
            legend('Cy5','Cy3');
            axis1.YLim = make_scale([Cy3(i,idx_G) Cy5(i,idx_G)]);
            %     ,'Location','NorthEastOutside'        
               
            hold off
            
            axis2 = subplot(3,1,2);
            Cy5dir_plot = plot(time_R, Cy5(i,idx_R),'r');
            hold on
            xlabel('Time, seconds')
            ylabel('Intensity')
            title('Red excitation')
            grid on
            vline(flow_time);
            axis2.YLim = make_scale(Cy5(i,idx_R));
            hold off
            axis3 = subplot(3,1,3);
            dyesum_11 = Cy5(i,idx_G)+Cy3(i,idx_G);
%             [rows] = find(dyesum_11 <fretthreshold);
            FRET_11 = Cy5(i,idx_G)./dyesum_11;
%             FRET_11(rows) = 0;

            FRET_plot = plot(time_G, FRET_11,'k');
            hold on
            legend('Cy5/(Cy5+Cy3)');
            xlabel('Time, seconds');
            ylabel('FRET');
           
            set(gca,'ytick',-.4:.2:1.4);
            grid on
             vline(flow_time);
            hold off
            axis1.XLim = [0 time(end)];
            axis2.XLim = [0 time(end)];
            axis3.XLim = [0 time(end)];
            axis3.YLim = [-.4 1.4];
        else
            if ~isempty(click_lines)
                for j = 1:length(click_lines)
                    delete(click_lines{j});
                end
                click_lines = cell(0);
            end
            clicks = click_points{i};
            if ~isempty(clicks)
                nclicks = length(clicks);
                click_lines = cell(1,nclicks);
                for j = 1:nclicks
                    click_lines{j} = vline(time(clicks(j)),'g-');
                end
            end
            Cy5_plot.YData = Cy5(i,idx_G);
            Cy5dir_plot.YData = Cy5(i,idx_R);
            axis1.YLim = make_scale([Cy3(i,idx_G) Cy5(i,idx_G)]);
            axis2.YLim = make_scale(Cy5(i,idx_R));
            Cy3_plot.YData = Cy3(i,idx_G);
            
            dyesum_11 = Cy5(i,idx_G)+Cy3(i,idx_G);
%             [rows] = find(dyesum_11 <fretthreshold);
            FRET_11 = Cy5(i,idx_G)./dyesum_11;
%             FRET_11(rows) = 0;
            FRET_plot.YData = FRET_11;
            axis3.YLim = [-.4 1.4];
            
        end
        
%     end
    B = 2;
            
    %1 = Analyze molecule
    %2 = Next molecule
    %3 = Previous molecule
    %4 = Accept molecule
    
    d = dialog('Position',[100 800 130 240],'Name','My Dialog');

        uicontrol('Parent',d,...
           'Style','text',...
           'Position',[5 200 120 30],...
           'String','Analyze this molecule?');
        uicontrol('Parent',d,...
           'Position',[5 170 120 30],...
           'String','Yes (f)',...
           'Callback',@callback1);
        uicontrol('Parent',d,...
           'Position',[5 140 120 30],...
           'String','Next (d)',...
           'Callback',@callback2);
        uicontrol('Parent',d,...
           'Position',[5 110 120 30],...
           'String','Previous (s)',...
           'Callback',@callback3);
        uicontrol('Parent',d,...
           'Position',[5 80 120 30],...
           'String','Accept (a)',...
           'Callback',@callback4);
        uicontrol('Parent',d,...
           'Position',[5 50 120 30],...
           'String','Save fig (r)',...
           'Callback',@callback6);
        uicontrol('Parent',d,...
           'Position',[5 20 120 30],...
           'String','End program',...
           'Callback',@callback5);
    
    
    d.KeyPressFcn=@key_pressed_fcn;
    uiwait(d);
    
    
%     B = menu('analyze this molecule?', 'yes', 'next' , 'previous', 'accept.', 'end program');
    
    if B == 1
        
        
%         A = menu('Which dye molecule bleached', 'Cy7', 'Cy5', 'Cy5 & Cy7', 'Cy3, Cy5 & Cy7', 'Neither');
       
        [xt,yt] = ginput;

        %%% prevent program from crashing if area outside graph is clicked
        if isempty(xt)
            xt = -10
        end
        
        while (max(xt) > time(end)) || (min(xt) < 0)
            %                     'Try again - reclick both points'
            [xt,yt] = ginput;
            if isempty(xt)
                xt = -10
            end
        end


%                %round after multiplying to get closer to actual frames selected
%             x2_1 = round(x(2)*rate_fr);
%             x1_2 = round(x(1)*rate_fr);
%             x2_2 = round(x(2)*rate_fr);

        %             x1_1 = round(x(1)*rate_fr*vsel_1);
        %             x2_1 = round(x(2)*rate_fr*vsel_1);
        %             x1_2 = round(x(1)*rate_fr*vsel_2);
        %             x2_2 = round(x(2)*rate_fr*vsel_2);

%             
%             Cy5_norm_1      = Cy5(i,idx_G);
%             norm_1 = mean(Cy5_norm_1(idx_norm));
%             Cy5_norm_1 = Cy5_norm_1 -norm_1;
%             
%             Cy3_norm_1      = Cy3(i,idx_G);
%             norm_1 = mean(Cy3_norm_1(idx_norm));
%             Cy3_norm_1 = Cy3_norm_1 -norm_1;


%             Cy5_plot.YData = Cy5_norm_1;
%             Cy3_plot.YData = Cy3_norm_1;
%             dyesum_11 = Cy5_norm_1+Cy3_norm_1;
%             FRET_11 = Cy5_norm_1./dyesum_11;
%             FRET_plot.YData = FRET_11;

%             C = menu('Accept this trace', 'yes', 'no');
        C=2;
        if ~isempty(click_lines)
            for j = 1:length(click_lines)
                delete(click_lines{j});
            end
            click_lines = cell(0);
        end
        clicks = round(xt*frame_rate)+1;
        if mod(length(clicks),2) == 1
           clicks = [1; clicks];
        end
        if ~isempty(clicks)
            nclicks = length(clicks);
            click_lines = cell(1,nclicks);
            for j = 1:nclicks
                click_lines{j} = vline(time(clicks(j)),'g-');
            end
        end
        d = dialog('Position',[100 800 130 100],'Name','My Dialog');

            uicontrol('Parent',d,...
                'Style','text',...
                'Position',[5 70 120 20],...
                'String','Save this trace?');
            uicontrol('Parent',d,...
                'Position',[5 40 120 30],...
                'String','Yes (f)',...
                'Callback',@callback1_1);
            uicontrol('Parent',d,...
                'Position',[5 10 120 30],...
                'String','No (d)',...
                'Callback',@callback1_2);
        d.KeyPressFcn=@key_pressed_fcn2;
        uiwait(d);

        if C == 1
                index = [index i];
                click_points{i} = clicks;
%                     Cy5(i,:) = Cy5_norm_1;
%                     Cy3(i,:) = Cy3_norm_1;

        end
            

        
        i = i+1;
        
    elseif B == 2
        i = i+1;
        
    elseif B == 3
        i = i-1;
        
    elseif B == 4
       
        index = [index i];
        click_points{i} = [1 length(time)];
        i = i+1;
    elseif B == 5
        break
    
    elseif B == 6
        saveas(gcf, strcat(path1,fileName, '_trace_', num2str(i), '.fig'));
%         i = i+1;
    end
end

% if ~COMMON_START_FILTER
%     close(h);
%     % wait_times
% end


%% save all data in workspace for later manipulation

r = menu('Would you like to save this workspace?','yes','no');


close(h)


if r == 1
    %     uy=path1
    %     cd(path1);
    [file,path2,filterindex] = uiputfile(strcat(path1,fileName(1:length(fileName)-4),'_weeded.mat'),fileName);
    %filterindex is set to zero if cancel button is clicked (or an error occurs), otherwise it's 1
    if filterindex
        if save_accepted_only
            Cy3 = Cy3(index,:);
            Cy5 = Cy5(index,:);
            if not(isempty(Seq))
                Seq = Seq(index,:);
            end
            if not(isempty(Dist))
                Dist = Dist(index);
            end
            if not(isempty(x))
                x = x(index);
            end
            if not(isempty(y))
                y = y(index);
            end
            if not(isempty(x_FQ))
                x_FQ = x_FQ(index);
            end
            if not(isempty(y_FQ))
                y_FQ = y_FQ(index);
            end
            if not(isempty(exp_name))
               exp_name = exp_name(index);
            end
            if not(isempty(ALEX))
               ALEX = ALEX(index);
            end
            if not(isempty(FOV))
               FOV = FOV(index);
            end
            if not(isempty(ID))
               ID = ID(index);
            end
            if not(isempty(click_points))
               click_points = click_points(index);
            end
        end

        save_file_name = strcat(path2,file);
        % Ugly try-catch solution to ugly problem...
        saved = 0;
        while saved == 0
            try
                disp(['Attempting to save file:' save_file_name]);
                save(save_file_name,'Cy3', 'Cy5', 'Seq','time','Dist','x','y','x_FQ','y_FQ','exp_name','ALEX','FOV','ID','click_points');
                saved = 1;
            catch
                warning('Problem saving, retry saving...')
                pause(.5);
            end
        end

    end
    
    %     cd(WD);
end
fclose('all');
clearvars -except Cy3 Cy5 Seq time Dist x y x_FQ y_FQ exp_name ALEX FOV ID click_points
end