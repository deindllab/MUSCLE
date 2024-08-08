    function cleanFileName = cleanFileName(rawFileName)
        % Replace spaces with underscores
        cleanFileName = strrep(rawFileName, ' ', '_');

        % Remove characters not allowed in filenames
        illegalChars = {'.', '/', '\', ':', '*', '?', '"', '<', '>', '|'}; % You can add more if needed
        for k = 1:length(illegalChars)
            cleanFileName = erase(cleanFileName, illegalChars{k});
        end
    end