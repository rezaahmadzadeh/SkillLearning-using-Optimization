%compute aggregate stats for the LASA dataset

LASAResultsFolderName = 'results\LASA_dataset';
ext = 'mat';
doPlot = 0; 

[learnedModels,nModels] = getAllMatFiles(LASAResultsFolderName, ext);

LASA_stats = cell(1,5);

LASA_stats{1}.algoName = 'position-only';
LASA_stats{2}.algoName = 'tangent-only';
LASA_stats{3}.algoName = 'laplace-only';
LASA_stats{4}.algoName = 'uniform weighting';
LASA_stats{5}.algoName = 'optimal weighting';

for j = 1:5
    LASA_stats{j}.performanceMeasures.SEA.list = [];
    LASA_stats{j}.performanceMeasures.SSE.list = [];
    LASA_stats{j}.performanceMeasures.DTWD.list = [];

end

for i = 1:nModels
    disp(['computing results for ' learnedModels{i} '...'])
    [evaluationResults,repros] = compare_with_baslines(learnedModels{i},doPlot);
    save([LASAResultsFolderName '\' erase(learnedModels{i},'trained.mat') 'evaluation.mat'],'evaluationResults','repros')
    
    for j = 1:length(evaluationResults)
        LASA_stats{j}.performanceMeasures.SEA.list = [LASA_stats{j}.performanceMeasures.SEA.list evaluationResults{j}.performanceMeasures.SEA.list];
        LASA_stats{j}.performanceMeasures.SSE.list = [LASA_stats{j}.performanceMeasures.SSE.list evaluationResults{j}.performanceMeasures.SSE.list];
        LASA_stats{j}.performanceMeasures.DTWD.list = [LASA_stats{j}.performanceMeasures.DTWD.list evaluationResults{j}.performanceMeasures.DTWD.list];
    end
    
    clear evaluationResults
end

for j = 1:5
    LASA_stats{j}.performanceMeasures.SEA.mean = mean(LASA_stats{j}.performanceMeasures.SEA.list);
    LASA_stats{j}.performanceMeasures.SEA.std = std(LASA_stats{j}.performanceMeasures.SEA.list);
    LASA_stats{j}.performanceMeasures.SSE.mean = mean(LASA_stats{j}.performanceMeasures.SSE.list);
    LASA_stats{j}.performanceMeasures.SSE.std = std(LASA_stats{j}.performanceMeasures.SSE.list);
    LASA_stats{j}.performanceMeasures.DTWD.mean = mean(LASA_stats{j}.performanceMeasures.DTWD.list);
    LASA_stats{j}.performanceMeasures.DTWD.std = std(LASA_stats{j}.performanceMeasures.DTWD.list);
end

save([LASAResultsFolderName '\aggregate_results\LASA_aggregregate_evaluation.mat'],'LASA_stats');