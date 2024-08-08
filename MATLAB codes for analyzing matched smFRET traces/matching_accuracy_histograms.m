%% Plotting an average FRET histogram for the dataset
stop = 10;
Cy3_B = mean(Cy3(:,2:2:stop),2);
Cy5_B = mean(Cy5(:,2:2:stop),2);
Cy5dir_B = mean(Cy5(:,1:2:stop),2);


idx = Cy5dir_B > 200; % filtering out traces without Cy5
Cy3_B = Cy3_B(idx);
Cy5_B = Cy5_B(idx);
FRET = Cy5_B./(Cy3_B + Cy5_B);

histogram(FRET, -0.2:0.1:1.2);