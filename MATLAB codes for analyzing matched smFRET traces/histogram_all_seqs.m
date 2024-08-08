
%These scripts have been modified to make the histograms look nice

%%%%% COVERAGE %%%%%%%%%%%%
% Define the sequence codes
code = ['ATGC';'ATGC';'TACG';'TACG'];

% Define the x and y axis ranges
x_axis = 1:16;
y_axis = 1:16;

% Initialize the heatmap matrix with zeros
htmap = zeros(16,16);

% Define the indices of the sequences to be analyzed
idx = [39,40,59,60];
% idx = 27:30; % Ha

% Extract the selected sequences from the full dataset
Seq1 = Seq(:,idx);

% Define Threshold for the number of sequences

threshold = 5;

% Iterate over the selected sequences

for i = 1:length(Seq1)
   % Extract the current sequence
   s = Seq1(i,:);
   
   % Map the current sequence to x,y coordinates on the heatmap
   x = 4*(strfind(code(4,:),s(4))-1)+(strfind(code(3,:),s(3))); % Hairpin
   y = 4*(strfind(code(1,:),s(1))-1)+(strfind(code(2,:),s(2)));
    
%      x = 4*(strfind(code(3,:),s(3))-1)+(strfind(code(1,:),s(1))); % Cas9
%      y = 4*(strfind(code(4,:),s(4))-1)+(strfind(code(2,:),s(2)));

     htmap(x,y) = htmap(x,y) + 1;

   
end

numberoftraces = htmap;
%% Plotting an array of histograms

sel_FRET = sel_FRET_B; % 10 deg data, change accordingly
% sel_FRET = sel_FRET_A;

scale = -0.2:0.05:1.2;
plot_num = 1;

figure('Position', [100, 100, 1000, 800]);

for j = 1:16
    for i = 1:16
        
        F_t = sel_FRET{i,j};
        F_t = F_t(~isnan(F_t));
        F_t = F_t(F_t>-0.2);
%         if numberoftraces(i,j) >= threshold
           subplot(16,16, plot_num);
           histogram(F_t, scale);
           box off
           set(gca,'YTickLabel',[]);
           set(gca,'XTickLabel',[]);
           set(gca,'XTick',[0 0.5 1]);
           set(gca,'YTick',[]);
           set(gca,'TickDir','out');
           ax=gca;
           axes('position',ax.Position,'box','on','ytick',[],'xtick',[],'color','none');
%         end
        plot_num = plot_num + 1;
    end
end

%% Saving the array of histograms

% Assuming your figure is the current figure
fig = gcf; % Get the handle to the current figure

% Export to PDF. 
%!Change the output path accordingly!
exportgraphics(fig, 'D:\hist_all.pdf', 'ContentType', 'vector');