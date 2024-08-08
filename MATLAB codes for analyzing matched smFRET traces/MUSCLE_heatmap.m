function MUSCLE_heatmap(htmap, scale, code, varargin)
    


    % Create an input parser object
    p = inputParser;

    % Add required inputs
    
    addRequired(p, 'htmap',  @(x) isnumeric(x) && ismatrix(x) && ~isvector(x));  % Check if htmap is a 2D numeric matrix
    addRequired(p, 'scale', @(x) isnumeric(x) && isvector(x) && numel(x)==2); % Check if scale is a numeric vector with 2 elements
    addRequired(p, 'code', @(x) ischar(x) && isequal(size(x),[4 4])); % Check if scale is a numeric vector with 2 elements
    
    % Add optional parameters
    addParameter(p, 'cmap', flipud(summer),@(x) isnumeric(x) && size(x,2) == 3 && all(x(:) >= 0) && all(x(:) <= 1)); % Check if colormap is valid
    addParameter(p, 'hmTitle', '', @ischar); % check if hmTitle is char
    addParameter(p, 'boxLabelsFontSize', 0, @(x) isnumeric(x) && numel(x)==1); % Optional flag to choose whether to add text labels to heatmap boxes
    addParameter(p, 'saveFolder', '', @(x) ischar(x) && isfolder(x)); % Check if input is a valid directory path

    % Parse the inputs
    parse(p, htmap, scale, code, varargin{:});

    % Access the inputs
    htmap = transpose(p.Results.htmap);
    scale = p.Results.scale;
    code = p.Results.code;
    cmap = p.Results.cmap;
    hmTitle = p.Results.hmTitle;
    boxLabelsFontSize = p.Results.boxLabelsFontSize;
    saveFolder = p.Results.saveFolder;
    
    %% Plot heatmap
    
    if ~isempty(hmTitle)
        fig = figure('Name',hmTitle);
    else
        fig = figure();
    end
    
    h = imagesc(htmap,scale);
    %% Format heatmap
    
    if ~isempty(hmTitle)
        title(hmTitle, 'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Arial');
    end
    colormap(cmap);
    c = colorbar('Ticks', scale, 'Location', 'eastoutside');
    % Specify the color of the tick labels
    c.Color = 'Black'; % Make last tick label red
    c.LineWidth = 1;
    alpha_data = ~isnan(htmap);
    axis square
    set(h, 'AlphaData', alpha_data);
    ylabels_pos2 = {code(2,1), code(2,2), code(2,3), code(2,4)};
    ylabels_pos2 = [ylabels_pos2 ylabels_pos2  ylabels_pos2  ylabels_pos2];
    xlabels_pos1 = {code(1,1), code(1,2), code(1,3), code(1,4)};
    xlabels_pos1 = [xlabels_pos1 xlabels_pos1 xlabels_pos1 xlabels_pos1];
    % xlabels_pos3 =[code(3,1),'   ', code(3,2),'   ', code(3,3),'   ', code(3,4)];
    % ylabels_pos4 =[code(4,1),'   ', code(4,2),'   ',code(4,3),'   ', code(4,4)];

    ylabels_pos4 ={code(4,1), code(4,2), code(4,3), code(4,4)};
    xlabels_pos3 ={code(3,1), code(3,2),code(3,3), code(3,4)};
    ax = gca;
    set(ax, 'TickLength', [0, 0]);
    set(ax, 'LineWidth', 1);
    
    % Set first level of axis labels

    set(ax, 'XTick', 1:16, 'XTickLabel', xlabels_pos1, 'TickLabelInterpreter', 'none','color', 'k', 'FontSize', 12, 'FontName', 'Arial');
    set(ax, 'YTick', 1:16, 'YTickLabel', ylabels_pos2, 'TickLabelInterpreter', 'none','color', 'k', 'FontSize', 12, 'FontName', 'Arial');
    % Change the color of the x and y tick labels to red
    set(ax, 'XColor', 'k', 'YColor', 'k');
    
    label_offset = 2.25;
    text(2.5:4:16,[16+label_offset 16+label_offset 16+label_offset 16+label_offset],xlabels_pos3, 'color', 'k','HorizontalAlignment' , 'center', 'FontSize', 14, 'FontName', 'Arial', 'FontWeight', 'bold');
    label_offset = 1;
    text([-label_offset -label_offset -label_offset -label_offset], 2.5:4:16,ylabels_pos4, 'color', 'k','HorizontalAlignment' , 'center', 'FontSize', 14, 'FontName', 'Arial', 'FontWeight', 'bold');

    ax.Title.FontSize = 18;
    if boxLabelsFontSize>0
        % Add values to boxes
        [R, C] = ndgrid(1:16, 1:16);
        % check if htmap contains any non-integer values
        if any(htmap ~= fix(htmap))
            boxlabels = arrayfun(@(x) sprintf('%.2f', x), htmap, 'UniformOutput', false);
            boxlabels = string(boxlabels);
        else
            boxlabels = cellstr(string(htmap));
        end

        text(C(:), R(:), boxlabels(:), 'color', 'k', 'HorizontalAlignment' , 'center', 'VerticalAlignment', 'middle', 'FontSize', boxLabelsFontSize, 'FontName', 'Arial');
    end
    
    %% Saving the files
    if ~isempty(saveFolder)
            % Create a timestamp string
            timestamp = datestr(now, 'yyyymmdd_HHMMSS');
            % Define base file name
            baseFileName = 'HM_';
        if ~isempty(hmTitle)
            baseFileName = [baseFileName, cleanFileName(hmTitle),'_'];
        end    
            
        % Add timestamp to base file name
        figFileName = fullfile(saveFolder, [baseFileName, timestamp, '.fig']);
        pdfFileName = fullfile(saveFolder, [baseFileName, timestamp, '.pdf']);

        % Save figure as .fig
        saveas(fig, figFileName);

        % Save figure as .pdf
        exportgraphics(fig,pdfFileName,'ContentType','vector');
            
    end

end