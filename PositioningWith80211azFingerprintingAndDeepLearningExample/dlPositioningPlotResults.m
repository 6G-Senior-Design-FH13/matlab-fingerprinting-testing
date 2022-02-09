function metric = dlPositioningPlotResults(mapFileName, labels, YPred, task)
%dlPositioningPlotResults Visualize indoor location estimation predictions
%   METRIC = dlPositioningPlotResults(MAPFILENAME, LABELS, YPRED, TASK)
%   visualises results of location estimation in a 3D environment modeled
%   in an stl file, MAPFILENAME. Plot the STAs true position according to
%   LABELS and display the predicted location, YPRED. These will be a
%   position vector of length 3 or a categorical location label depending
%   on the value of TASK.
%
%   METRIC is the mean position error when TASK = "positioning" and the
%   mean accuracy when TASK = "localization".

%   Copyright 2020 The MathWorks, Inc. 


YTrue = labels.classification;
locs = labels.regression;
if ~isa(mapFileName, 'triangulation')
    tri = stlread(mapFileName);
else
    tri = map;
end
figure
trisurf(tri, ...
    'FaceAlpha', 0.3, ...
    'FaceColor', [.5, .5, .5], ...
    'EdgeColor', 'none');
view(60, 30);
hold on; axis equal; grid off;
xlabel('x'); ylabel('y'); zlabel('z');
view([84.75 56.38])
% Plot edges
fe = featureEdges(tri,pi/20);
numEdges = size(fe, 1);
pts = tri.Points;
a = pts(fe(:,1),:); 
b = pts(fe(:,2),:); 
fePts = cat(1, reshape(a, 1, numEdges, 3), ...
    reshape(b, 1, numEdges, 3), nan(1, numEdges, 3));
fePts = reshape(fePts, [], 3);
plot3(fePts(:, 1), fePts(:, 2), fePts(:, 3), 'k', 'LineWidth', .5); 

if task =="positioning"
    positionError = plotRegression(locs, YPred);
    metric = positionError;
else
    classAccuracy = plotClassification(YTrue, YPred, locs);
    metric = classAccuracy;
end

end

function accuracy = plotClassification( YTrue, YPred, locs)
%PLOTCLASSIFICATION Visualize indoor localization results
%   plotClassification(YTRUE, YPRED, LOCS) produces a 3D scatter plot of
%   true STA positions, LOCS, on a map. These are colored corresponding
%   to their predicted class label, (YPRED). The true labels, YTRUE are 
%   then compared to the predictions and a confusion chart and percent 
%   accuracy are displayed.

% Convert class to "one-hot" logical vector
    numClasses = length(unique(YTrue));
    title({'{\bf{CNN Location Prediction}}';'\fontsize{10}True STA positions coloured by Predicted Class'},'FontWeight','Normal');

    % Visualize the class areas
    % Conference Room
    v1= [0 2.75 0.05; 0 0 0.05; 5 0 0.05; 5 2.75 0.05];
    f1 = [1 2 3 4];
    patch('Faces', f1, 'Vertices', v1, 'FaceColor', 'B', 'FaceAlpha', 0.3);

    % Storage
    v2= [0 2.75 0.05; 0 8.0 0.05; 0.5 8 0.05; 0.5 2.75 0.05];
    f2 = [1 2 3 4];
    patch('Faces', f2, 'Vertices', v2, 'FaceColor', 'Y', 'FaceAlpha', 0.3);

    % Desk 1 
    v3= [0.5 6.5 0.05; 0.5 8.0 0.05; 2.5 8 0.05; 2.5 6.5 0.05];
    f3 = [1 2 3 4];
    patch('Faces', f3, 'Vertices', v3, 'FaceColor', 'G', 'FaceAlpha', 0.3);

    % Desk 2
    v4= [2.5 6.5 0.05; 2.5 8.0 0.05; 5 8 0.05; 5 6.5 0.05];
    f4 = [1 2 3 4];
    patch('Faces', f4, 'Vertices', v4, 'FaceColor', 'M', 'FaceAlpha', 0.3);

    % Desk 3
    v5= [3.5 6.5 0.05; 5 6.5 0.05; 5 4.5 0.05; 3.5 4.5 0.05];
    f5 = [1 2 3 4];
    patch('Faces', f5, 'Vertices', v5, 'FaceColor', 'C', 'FaceAlpha', 0.3);

    % Desk 4
    v6= [3.5 4.5 0.05; 5 4.5 0.05; 5 2.75 0.05; 3.5 2.75 0.05];
    f6 = [1 2 3 4];
    patch('Faces', f6, 'Vertices', v6, 'FaceColor', 'R', 'FaceAlpha', 0.3);

    % Office (white) not plotted - remaining space
    locations = sort(categories(YTrue));
    
    % Visualise the predicted receiver classes
    colors = {'b','g','m', 'c' ,'r','w','y'};

    for p=1:size(YTrue,2)
        x = colors((YPred(p)==locations));
        scatter3(locs(1,p), locs(2,p), locs(3,p), char(x), 'filled');
    end

    l = zeros([numClasses 1]);
    for q = 1:numClasses
        l(q) = scatter(NaN, NaN, char(colors(q)), "filled", 'MarkerEdgeColor', 'k');
    end
    locations{1}="conference room"; %removes _ from label name for display
    lgd = legend(l, locations, "location", "bestoutside");
    lgd.Title.String = "Predicted Class";

    % Plot the confusion chart and Accuracy
    figure
    cm = confusionchart(cellstr(YTrue),cellstr(YPred));
    accuracy = sum(diag(cm.NormalizedValues))/sum(cm.NormalizedValues(:))*100;
    cm.title(['Accuracy: ',num2str(accuracy),'%']);
end

function meanErr = plotRegression(YTrue, YPred)
%PLOTREGRESSION Visualize indoor localization results
%   plotRegression(YTRUE, YPRED) produces a 3D scatter plot of
%   true STA positions, YTRUE, on a map. These are colored corresponding
%   to magnitude of distance error of predicted positions, (YPRED). The
%   calculated error is used to produce a CDF plot.

    YTrue = YTrue';
    mErr = zeros([1 size(YTrue, 1)]);
    
    % Compute the distance error
    for i=1:size(YTrue,1)
        mErr(i) = double(norm(YTrue(i,:) - YPred(i,:)));
    end

    % Set the color bar properties
    minErr = floor(min(mErr))*10;
    maxErr = ceil(max(mErr))*10;
    numColors =  (maxErr - minErr)/5;
    cm = colormap(jet(numColors));

    % Plot the true receiver locations - coloured by magnitude of error
    for i = 1:size(YTrue,1)
        cmIdx = find((mErr(i)*10 - (minErr:5:maxErr))<0, 1) - 1;
        if cmIdx >numColors
            cmIdx = numColors;
        end
        scatter3(YTrue(i, 1), YTrue(i, 2), YTrue(i, 3), 'MarkerEdgeColor', cm(cmIdx,:), 'MarkerFaceColor', cm(cmIdx,:), 'MarkerFaceAlpha', 1.0);
    end

    % Create colorbar
    cb = colorbar; % ('direction', 'reverse');
    cb.Label.String = 'Distance Error (m)';
    cbLim = cb.Limits;
    cb.Ticks = cbLim(1) + diff(cbLim)/(2*numColors) + ...
        (0:numColors-1)*diff(cbLim)/numColors;
    cb.TickLabels = (minErr/10):.5:(maxErr-5/10);
    title({'{\bf{CNN Position Prediction}}';'\fontsize{10}True STA Positions coloured by Distance Error'},'FontWeight','Normal');
    
    % Plot CDF and display mean distance error
    figure
    stairs(sort(mErr),(1:length(mErr))/length(mErr));
    meanErr = mean(mErr);
    grid on;
    xlabel('Distance error (m)') 
    ylabel('Cumulative Probability') 
    title(['CDF - Positioning Error  ', '(Mean = ',num2str(round(meanErr,2)),'m)']);

end