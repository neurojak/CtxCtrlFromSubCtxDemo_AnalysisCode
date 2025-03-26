% Jason Keller
% 2020-2025
% 
% top level script for analyses in "Cortical control of innate behavior from subcortical demonstration":
% -----------------------------------------------------------------------------------------------------
% this script includes all analyses for creating figures, but other scripts are required for
% preprocessing raw data; the full analysis script loads raw data based on a data hierarchical
% structure: project->mouse->experiment, and some of the functions that support this functionality 
% are included here for clarity although  never called, as well as many functions for exploratory 
% data analyses
%
% this script loads a MAT file (which can be 10s of GB, so needs a decent amount of RAM >60GB), to illustrates
% the analyses methods used

%% program control:
behaviorData = 0; %set to 0/false to skip behavior data, etc.
kinematicData = 1;
ephysData = 0;

%% other constants
blueCmap = [0 0.4470 0.7410];
orangeCmap = [0.8500 0.3250 0.0980];	
yellowCmap = [0.9290 0.6940 0.1250];
medGreyCmap = [0.4 0.4 0.4];
lightGreyCmap = [0.7 0.7 0.7];
linesCmap = lines;
fontSz = 16;
maxTrialsPlotted = 50; %standard learning curve ~150 (must load multiple sessions)

rootFolder = 'D:\matFiles\';
projectFilename = 'exampleData.mat';

%% load data & plot
% load([rootFolder projectFilename]);

if behaviorData % set behaviorData to 0 or comment out plots below if desired
%% BEHAVIOR: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plotBasicVhexDataNorm(tsData_test, idxData_test); %normalized to see timing relationships better
%         plotBasicVhexData(tsData_test, idxData_test);

% LOG latency vs. trials combined no opto %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cmapi = lines(length(mouseList.ID));
figure; hold on;
maxJumps = 1;
%maxTrialsPlotted = 150; %standard learning curve ~150
shiftSec = 3; %add delay from platform to cold air to ensure no negative log
trialWindow = 5; %for Gaussian smoothing
for i = 1:length(mouseList.ID) %loop over mice for individual traces
   dataToUse = latencyVsJumps; % latencyVsJumpsNonOpto option to switch to non opto only

   currMouseDataToUse = dataToUse{i};
   currentOptoTrials = optoIndVsJumps{i};
%    currMouseDataToUse(currentOptoTrials) = []; %option to remove opto trials

   currentMaxJumps = length(currMouseDataToUse);
   if currentMaxJumps > maxJumps
       maxJumps = currentMaxJumps; %keep track of max for padding below
   end 

   latencyVsJumpsShifted{i} = smoothdata(currMouseDataToUse + shiftSec, 'gaussian', trialWindow); %smooth to get rid of NaN gaps
    
   %also loop over experiments to mark session boundaries:
   if isempty(mouseList.expList{i})
       currentExpList = 1:length(find(numJumpsAcrossExp(i,:))); 
   else
       currentExpList = mouseList.expList{i};
   end
   for j = currentExpList
       currentExpJumps = numJumpsAcrossExp(i,j);

       if j == currentExpList(1)
           startJump = 1;
       end
       if startJump < maxTrialsPlotted
%         semilogy(startJump, latencyVsJumpsShifted{i}(startJump), 'LineStyle', 'none', 'Marker', 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'm', 'MarkerSize', 7) %mark start of day/experiment with diff color
       end
       startJump = startJump+currentExpJumps; %update start jump for next exp
   end
end

%pad latency data with NaNs to make matrix for group statistics for bounded line:
latencyVsJumpsNonOptoPadded = latencyVsJumpsShifted;
for i = 1:length(mouseList.ID)
   currMouseDataToUse = dataToUse{i};
   currentOptoTrials = optoIndVsJumps{i};
%    currMouseDataToUse(currentOptoTrials) = []; %option to remove opto trials
    currentNumJumps = length(currMouseDataToUse);
%     currentNumJumps = length(latencyVsJumps{i}(~avoidIndVsJumps{i})); %option to compare only react (see change above as well)
    if currentNumJumps <= maxJumps
       padLength =  maxJumps - currentNumJumps + 1;
       latencyVsJumpsNonOptoPadded{i}(end:end+padLength-1) = NaN(padLength,1)';
    end
end
latencyVsJumpsNonOptoMat = cell2mat(latencyVsJumpsNonOptoPadded');

[meanLat,varLat] = grpstats(latencyVsJumpsNonOptoMat,[],{'mean' 'sem'});
x = 1:1:length(meanLat);
% semilogy(x,meanLat, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 2 )
[hlA,hpA] = boundedline(x(1:maxTrialsPlotted), meanLat(1:maxTrialsPlotted), varLat(1:maxTrialsPlotted), 'cmap', [0.1 0.1 0.1], 'alpha','y', 'transparency', 0.5, 'LineWidth', 2); 
set(gca,'yscale','log')

ylabel('extension latency (sec)', 'FontSize', fontSz);
xLims = get(gca, 'XLim');
semilogy(xLims, [shiftSec,shiftSec],'Color', 'r', 'LineStyle', '--', 'LineWidth', 1); %dotted line for latency=0 (cold air on)
set(gcf,'color','w');
xlabel('trials', 'FontSize', fontSz); 
hold off;
ylim([0.1 40])

% cumulative sum of avoids over trials %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; hold on;
%maxTrialsPlotted = 150; %standard learning curve ~150
maxJumps = 1; %keep track for padding below
for i = 1:length(mouseList.ID) %loop over mice

   currentData = latencyVsJumps{i};
   currentOptoTrials = optoIndVsJumps{i};
%    currentData(currentOptoTrials) = []; %option to remove opto trials to match above data
   cumulativeAvoids{i} = cumsum(currentData<0); % cumulative sum of avoid trials where latency < 0;

   currentMaxJumps = length(cumulativeAvoids{i});
   if currentMaxJumps > maxJumps
       maxJumps = currentMaxJumps; %keep track of max across mice for padding below
   end 

   if currentMaxJumps >= maxTrialsPlotted
        plot(1:1:maxTrialsPlotted, cumulativeAvoids{i}(1:maxTrialsPlotted), 'Color', lightGreyCmap); %plot individual mouse traces
   else
        plot(1:1:currentMaxJumps, cumulativeAvoids{i}, 'Color', lightGreyCmap); %plot individual mouse traces
   end
end

%pad latency data with NaNs to make matrix for group statistics for bounded line:
cumulativeAvoidsPadded = cumulativeAvoids;
for i = 1:length(mouseList.ID)
    currentNumJumps = length(cumulativeAvoids{i});
    if currentNumJumps <= maxJumps
       padLength =  maxJumps - currentNumJumps + 1;
       cumulativeAvoidsPadded{i}(end:end+padLength-1) = NaN(padLength,1)';
    end
end

%pad latency data with NaNs to make matrix for group statistics for bounded line:
cumulativeAvoidsPadded = cumulativeAvoids;
for i = 1:length(mouseList.ID)
    currentNumJumps = length(cumulativeAvoids{i});
    if currentNumJumps <= maxJumps
       padLength =  maxJumps - currentNumJumps + 1;
       cumulativeAvoidsPadded{i}(end:end+padLength-1) = NaN(padLength,1)';
    end
end
cumulativeAvoidsMat = cell2mat(cumulativeAvoidsPadded');
% cumulativeAvoidsMax = max(cumulativeAvoidsMat,[],2); %optionally save this variable for later comparing different groups

[meanCumulativeAvoids,varCumulativeAvoids] = grpstats(cumulativeAvoidsMat,[],{'mean' 'sem'});
x = 1:1:length(meanCumulativeAvoids);
[hlA,hpA] = boundedline(x(1:maxTrialsPlotted), meanCumulativeAvoids(1:maxTrialsPlotted), varCumulativeAvoids(1:maxTrialsPlotted), 'cmap', orangeCmap, 'alpha','y', 'LineWidth', 1); %#ok<*ASGLU>
xLims = get(gca, 'XLim');
% yLims = get(gca, 'YLim');
axis([xLims(1) maxTrialsPlotted 0 50]);
ylabel('cumulative # avoids', 'FontSize', fontSz);
xlabel('trials', 'FontSize', fontSz); 
xticks(0:50:150)
yticks(0:10:30)
hold off;
set(gcf,'color','w');

% TOTAL opto vs not latency distribution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; 
hold on;
latenciesOpto = [];
latenciesNoOpto = [];
shiftSec = 3; % for taking log

for i = 1:length(mouseList.ID)  
    currentOptoIndeces = optoIndVsJumps{i};
    currentLatenciesOpto = latencyVsJumps{i}(currentOptoIndeces);
    currentLatenciesNonOpto = latencyVsJumps{i}(~currentOptoIndeces);
    latenciesOpto = [latenciesOpto currentLatenciesOpto]; %#ok<*SAGROW>
    latenciesNoOpto = [latenciesNoOpto currentLatenciesNonOpto];
end

s1 = swarmchart(ones(1,length(latenciesNoOpto))*-0.5, latenciesNoOpto+shiftSec, 200, medGreyCmap, 'filled','MarkerFaceAlpha',0.4);
s2 = swarmchart(ones(1,length(latenciesOpto))*0.5, latenciesOpto+shiftSec,  200, blueCmap, 'filled','MarkerFaceAlpha',0.4);
xLims = get(gca, 'XLim');
plot(xLims, [shiftSec,shiftSec],'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 1); %dotted line for latency=0 (cold air on)
hold off;
set(gca,'yscale','log')
ylabel('latencies (sec)')
xlim([-1 1]);
ylim([0.4 10]);
xticks([-0.5 0.5])
xticklabels({'no opto','opto'})
makepretty;

pWilcoxRank = ranksum(latenciesOpto,latenciesNoOpto); %Wilcoxon rank sum test (Mann-Whitney U Test) for distributions
p = pWilcoxRank;
yLims = get(gca, 'YLim');
y = [yLims(1) yLims(2)];
if p < 0.001
   text(0, yLims(2)*0.9, '***', 'FontSize', 20);
elseif p < 0.01
   text(0, yLims(2)*0.9, '**', 'FontSize', 20);
elseif p < 0.05
   text(0, yLims(2)*0.9, '*', 'FontSize', 20);
else %p > 0.05
   text(0, yLims(2)*0.9, 'n.s.', 'FontSize', 20);
end

% TOTAL opto vs not avoid PROBABILITY distribution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; 
hold on;
probAvoidOpto = [];
probAvoidNoOpto = [];

for i = 1:length(mouseList.ID)  
    currentOptoIndeces = optoIndVsJumps{i};
    currentProbAvoidOpto = length(find(avoidIndVsJumps{i}(currentOptoIndeces))) / length(find(currentOptoIndeces)); %total # avoid extensions / total # extensions
    currentNonOptoIndeces = ~optoIndVsJumps{i};
    currentProbAvoidNonOpto = length(find(avoidIndVsJumps{i}(currentNonOptoIndeces))) / length(find(currentNonOptoIndeces)); 
    probAvoidOpto = [probAvoidOpto currentProbAvoidOpto]; %#ok<*SAGROW>
    probAvoidNoOpto = [probAvoidNoOpto currentProbAvoidNonOpto];
end

x1 = ones(1,length(probAvoidNoOpto))*-0.5;
x2 = ones(1,length(probAvoidOpto))*0.5;
scatter(x1, probAvoidNoOpto, 200, medGreyCmap, 'filled');
scatter(x2, probAvoidOpto, 200, blueCmap, 'filled');
plot([x1; x2], [probAvoidNoOpto; probAvoidOpto], 'LineStyle', '-', 'Color', medGreyCmap)
xlim([-1 1]);
xticks([-0.5 0.5])
xticklabels({'no opto','opto'})
ylabel('avoid probability')
set(gca, 'YLim', [0 1]);
makepretty;

pWilcoxSign = signrank(probAvoidOpto,probAvoidNoOpto); %Wilcoxon signed rank test for paired measurements
p = pWilcoxSign;
yLims = get(gca, 'YLim');
y = [yLims(1) yLims(2)];
if p < 0.001
   text(0, yLims(2)*0.9, '***', 'FontSize', 20);
elseif p < 0.01
   text(0, yLims(2)*0.9, '**', 'FontSize', 20);
elseif p < 0.05
   text(0, yLims(2)*0.9, '*', 'FontSize', 20);
else %p > 0.05
   text(0, yLims(2)*0.9, 'n.s.', 'FontSize', 20);
end

% opto vs not acceleration distribution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hAccelFig = figure; 
hold on;
accelOpto = [];
accelNoOpto = [];

for i = 1:length(mouseList.ID)  
    currentOptoIndeces = optoIndVsJumps{i}; % add this back in to test during high probability periods: & avoidProbVsJumps{i};
    currentAccelOpto = accelVsJumps{i}(currentOptoIndeces);
    currentNonOptoIndeces = ~currentOptoIndeces; % & ~avoidIndVsJumps{i}; % & avoidProbVsJumps{i}; %optionally only look at reactive for fair comparison
    currentAccelNonOpto = accelVsJumps{i}(currentNonOptoIndeces);
    accelOpto = [accelOpto currentAccelOpto]; %#ok<*SAGROW>
    accelNoOpto = [accelNoOpto currentAccelNonOpto];
end

data = {accelNoOpto; accelOpto};
plotSpread(data, 'xValues', [-0.6, 0.6], 'xNames', {'no opto','opto'}, 'spreadWidth', 1, 'spreadFcn', {'lin',[25]}, 'distributionColors', {medGreyCmap blueCmap});
ylabel('extension peak acceleration (cm/sec^2)')
set(gcf,'color','w');
makepretty;

pWilcox = ranksum(accelOpto,accelNoOpto); %Wilcoxon rank sum test (Mann-Whitney U Test)
p = pWilcox;
yLims = get(gca, 'YLim');
xLims = get(gca, 'XLim');
y = [yLims(1) yLims(2)];
if p < 0.001
   text(0, yLims(2)*0.9, '***', 'FontSize', 20);
elseif p < 0.01
   text(0, yLims(2)*0.9, '**', 'FontSize', 20);
elseif p < 0.05
   text(0, yLims(2)*0.9, '*', 'FontSize', 20);
else %p > 0.05
   text(0, yLims(2)*0.9, 'n.s.', 'FontSize', 20);
end
text(xLims(1)+0.1, yLims(2)*0.99, ['p = ' num2str(p,3)], 'FontSize', 20);

% avoid vs react acceleration distribution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hAccelFig = figure; 
hold on;
accelAvoid = [];
accelReact = [];

for i = 1:length(mouseList.ID)  
    currentAvoidIndeces = avoidIndVsJumps{i} & ~optoIndVsJumps{i}; 
    currentAccelAvoid = accelVsJumps{i}(currentAvoidIndeces);
    currentReactIndeces = ~avoidIndVsJumps{i} & ~optoIndVsJumps{i}; %
    currentAccelReact = accelVsJumps{i}(currentReactIndeces);
    accelAvoid = [accelAvoid currentAccelAvoid]; %#ok<*SAGROW>
    accelReact = [accelReact currentAccelReact];
end

data = {accelAvoid; accelReact};
plotSpread(data, 'xValues', [-0.6, 0.6], 'xNames', {'avoid','react'}, 'spreadWidth', 1, 'spreadFcn', {'lin',[25]}, 'distributionColors', {orangeCmap blueCmap});
ylabel('extension peak acceleration (cm/sec^2)')
set(gcf,'color','w');
makepretty;

pWilcox = ranksum(accelAvoid,accelReact); %Wilcoxon rank sum test (Mann-Whitney U Test)
p = pWilcox;
yLims = get(gca, 'YLim');
xLims = get(gca, 'XLim');
y = [yLims(1) yLims(2)];
if p < 0.001
   text(0, yLims(2)*0.9, '***', 'FontSize', 20);
elseif p < 0.01
   text(0, yLims(2)*0.9, '**', 'FontSize', 20);
elseif p < 0.05
   text(0, yLims(2)*0.9, '*', 'FontSize', 20);
else %p > 0.05
   text(0, yLims(2)*0.9, 'n.s.', 'FontSize', 20);
end
text(xLims(1)+0.1, yLims(2)*0.99, ['p = ' num2str(p,3)], 'FontSize', 20);


% mean accel across days %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; hold on;

for i = 1:length(mouseList.ID)

   if isempty(mouseList.expList{i})
       currentExpList = 1:length(find(numJumpsAcrossExp(i,:))); 
   else
       currentExpList = mouseList.expList{i};
   end
   for j = 1:4  %currentExpList %plot over each experiment separate to mark session boundaries
       if j == 1
           startJump = 1;
       end
       currentExpJumps = numJumpsAcrossExp(i,j);
       
       if ~isnan(currentExpJumps) 
           currentOptoIndeces = optoIndVsJumps{i}; 
           currentAccelAll = accelVsJumps{i};
           %currentAccelAll(currentOptoIndeces) = NaN(1,length(find(currentOptoIndeces))); %remove opto extensions for general accel calculation if applicable
           meanAccelByDay(i,j) = mean(currentAccelAll(startJump:startJump+currentExpJumps-1),'omitnan');
           startJump = startJump+currentExpJumps;
       else
           meanAccelByDay(i,j) = NaN;
       end
        
   end
   
   plot(meanAccelByDay(i,:), 'Color', lightGreyCmap);
   
end

[meanAccelDays,varAccelDays] = grpstats(meanAccelByDay,[],{'mean' 'sem'});
x = 1:1:length(meanAccelDays);
[hlA,hpA] = boundedline(x, meanAccelDays, varAccelDays, 'cmap', [0 0 0], 'alpha','y', 'LineWidth', 1); %#ok<*ASGLU>

ylabel({'mean peak extension','acceleration (cm/sec^2)'}, 'FontSize', fontSz);
xLims = get(gca, 'XLim');
x = [xLims(1) xLims(2)];
yLims = get(gca, 'YLim');
y = [yLims(1) yLims(2)];
plot(x, [0,0],'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 1); %dotted line for latency=0 (cold air on)
axis([1 4 0 160]);
set(gcf,'color','w');
xlabel('session #', 'FontSize', fontSz); 
hold off;



end %end if behaviorData

if kinematicData
%% KINEMATICS: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use pooled data across mice/experiments for summary plots: 
    
    kptNames = {'tail','r_toe','r_ankle','r_knee','l_toe','l_ankle','l_knee','r_paw','l_paw','nose','mouth','chest'}; %grouped by HL, FL, other
    %kptNamesSubset = {'l_knee','r_knee','l_ankle','r_ankle','l_toe','r_toe'};
    kptNamesSubset = {'l_knee'};

    % get median keypoints during ITI for stable centering & comparing across different conditions:
    [medianKpts] = getMedianItiKeypoints(kinVariableVsTrialsITI_allMice);

    for kptIdx = kptNames
        %plotKinHeatmapAcrossTrials(kinVariableVsTrialsAtVhex_allMice, kinVarTimescale, kptIdx{1}, notOptoTrialsPerMouse, medianKpts)
        kptTrialSimilarity.(kptIdx{1}) = plotKinCorrMatrixAcrossTrials(kinVariableVsTrialsAtVhex_allMice, kptIdx{1}, notOptoTrialsPerMouse);
    end
    
    % plot all diagonals as a measure of self-similarity of kinematics across trials
    plotKinCorrSimilarities(kptTrialSimilarity,kptNames);
    
    % joint angle calculations and plots:    
    plotKinAngleHeatmapAcrossTrials(kinAngleVsTrialsAtVhex_allMice, kinVarTimescale, kinParamName, jointNames, notOptoTrialsPerMouse) %Suppl Fig 2
    [syncKnee, syncAnkle, syncHLL, syncHLR] = plotKinAnglesCorrs(kinAngleVsTrialsAtVhex_allMice, kinVarTimescale, kinParamName, notOptoTrialsPerMouse, accelNotOptoTrials_allMice);
    
    % compare overall keypoint mean PSTHs across all mice and different events/conditions
     % @ extension
     %compareKinMeanHeatmaps(kinVariableVsTrialsAvoid_allMice, kinVariableVsTrialsReact_allMice, kinVarTimescale, avoidTrialsPerMouse, reactTrialsPerMouse, medianKpts)
     compareKinMeanHeatmaps(kinVariableVsTrialsReact_allMice, kinVariableVsTrialsOpto_allMice, kinVarTimescale, reactTrialsPerMouse, optoTrialsPerMouse, medianKpts)
     % @ cue
     compareKinMeanHeatmaps(kinVariableVsTrialsAvoidAtCue_allMice, kinVariableVsTrialsReactAtCue_allMice, kinVarTimescale, avoidTrialsPerMouse, reactTrialsPerMouse, medianKpts)
     compareKinMeanPositionsBeforeEvent(kinVariableVsTrialsAvoidAtCue_allMice, kinVariableVsTrialsReactAtCue_allMice, kinVarTimescale, avoidTrialsPerMouse, reactTrialsPerMouse, medianKpts) %Suppl Fig 3
    
%% plot PCA of kinematics:
    plotWindow = [-1 1]; %seconds around extension to plot trajectory
%     plotWindow = [-2.0 1.0]; %seconds around CUE to plot trajectory
    ptsPerMouse = ones(1,length(mouseNumsLoaded)).*36; %12kpts by X/Y/Z = 36 dimensions per mouse
    % @ VHEx:
    %plotKinAvoidVsReactPcTrajectoriesAllMice(allKptsMeanZ_allMice, squeeze(mean(avoidMeanZ_allMice,1)), squeeze(mean(reactMeanZ_allMice)), ptsPerMouse, timescaleSeg, pcaWindowVhex, plotWindow)
    plotKinAvoidVsReactPcDistancesAllMiceEqualTrials(allKptsMeanZ_allMice, avoidMeanZ_allMice, reactMeanZ_allMice, allKptsMeanZShuf1_allMice, allKptsMeanZShuf2_allMice, ptsPerMouse, timescaleSeg, pcaWindowVhex, plotWindow)
    % @ control vs ITI @ VHEx and cue
    plotKinAvoidVsReactPcDistancesAllMiceEqualTrials(allKptsMeanZ_allMice, avoidMeanZ_allMice, itiControl_allMice, allKptsMeanZShuf1_allMice, allKptsMeanZShuf2_allMice, ptsPerMouse, timescaleSeg, pcaWindowVhex, plotWindow)
    %plotKinAvoidVsReactPcDistancesAllMiceEqualTrials(allKptsAtCueMeanZ_allMice, avoidAtCueMeanZ_allMice, itiControl_allMice, allKptsMeanZShufAtCue1_allMice, allKptsMeanZShufAtCue2_allMice, ptsPerMouse, timescaleSeg, pcaWindowVhex, plotWindow)
    % @ cue:
    %plotKinAvoidVsReactPcDistancesAllMiceEqualTrials(allKptsAtCueMeanZ_allMice, avoidAtCueMeanZ_allMice, reactAtCueMeanZ_allMice, allKptsMeanZShufAtCue1_allMice, allKptsMeanZShufAtCue2_allMice, ptsPerMouse, timescaleSeg, pcaWindowCue, plotWindow)
    % shuffle control:
    %plotKinAvoidVsReactPcTrajectoriesAllMice(allKptsMeanZ_allMice, allKptsMeanZShuf1_allMice, allKptsMeanZShuf2_allMice, ptsPerMouse, timescaleSeg, pcaWindowVhex, plotWindow)


end % end kinematicData

if ephysData
%% EPHYS: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    % ephys constants 
    anatNames = {'CTX', 'BS', 'HB', 'MB', 'TH', 'HY'}; %, 'CNU', 'HPF' }; %based on Allen Reference taxonomy (see 'gui_data.st' in atlas script, https://atlas.brain-map.org/)
    %anatNames = {'Isocortex','CTX'};
    anatCmaps = linesCmap(1:length(anatNames),:); %corresponding colormaps to index for plotting
    % abbr.     name                    atlas taxonomy hierarchy depth
    % -----------------------------------------------------------------
    % BS        brainstem               2   (includes IB w/ TH & HY, plus MB & HB)
    % CTX       cerebral cortex         3   (also make "nonCTX" structure with compliment, as a special case)
    % CNU       cerebral nuclei         3
    % HB        hindbrain               3   (include pontine & medulla RF) 
    % MB        midbrain                3   (includes SC, PAG, SNr, VTA, RN, CUN, PPN, midbrain RF)
    % TH        thalamus                4
    % HY        hypothalamus            4
    % HPF       hippocampal formation   5

    % subset of anatomical areas on different probes to compare:
    probesSub = {'Probe0','Probe1','Probe0','Probe2'}; %matched pairs with below
    anatsSub =  {'CTX',   'CTX',   'BS',    'HB'};
    labelsSub = {'M1',    'PFC',   'TH/HY', 'HB'}; 
    lineCmaps = linspecer(4); %better colormap for only a few options

    % make structure that pools across all probes
    allProbes.allUnitsMeanZFR = [Probe0.allUnitsMeanZFR; Probe1.allUnitsMeanZFR; Probe2.allUnitsMeanZFR];
    allProbes.allUnitsAvoidMeanZFR = [Probe0.avoidMeanZFR; Probe1.avoidMeanZFR; Probe2.avoidMeanZFR];
    allProbes.allUnitsReactMeanZFR = [Probe0.reactMeanZFR; Probe1.reactMeanZFR; Probe2.reactMeanZFR];
    allProbes.allUnitsOptoMeanZFR = [Probe0.optoMeanZFR; Probe1.optoMeanZFR; Probe2.optoMeanZFR];
    allProbes.allUnitsItiMeanZFR = [Probe0.itiMeanZFR; Probe1.itiMeanZFR; Probe2.itiMeanZFR];
    allProbes.allUnitsAvoidAtCueMeanZFR = [Probe0.avoidAtCueMeanZFR; Probe1.avoidAtCueMeanZFR; Probe2.avoidAtCueMeanZFR];
    allProbes.allUnitsReactAtCueMeanZFR = [Probe0.reactAtCueMeanZFR; Probe1.reactAtCueMeanZFR; Probe2.reactAtCueMeanZFR];
    allProbes.allUnitsAvoidAtItiMeanZFR = [Probe0.avoidAtItiMeanZFR; Probe1.avoidAtItiMeanZFR; Probe2.avoidAtItiMeanZFR];
    allProbes.allUnitsReactAtItiMeanZFR = [Probe0.reactAtItiMeanZFR; Probe1.reactAtItiMeanZFR; Probe2.reactAtItiMeanZFR];
    allProbes.allUnitsOptoAtCueMeanZFR = [Probe0.optoAtCueMeanZFR; Probe1.optoAtCueMeanZFR; Probe2.optoAtCueMeanZFR];
    allProbes.allUnitsDiffAtCueMeanZFR = [Probe0.diffMeanZFR; Probe1.diffMeanZFR; Probe2.diffMeanZFR]; 
    allProbes.allUnitsDiffAtVhexMeanZFR = allProbes.allUnitsAvoidMeanZFR - allProbes.allUnitsReactMeanZFR;
    allProbes.allUnitsDiffAtItiMeanZFR = allProbes.allUnitsAvoidAtItiMeanZFR - allProbes.allUnitsReactAtItiMeanZFR;
    allProbes.allUnitsPlatformExtendMeanZFR = [Probe0.platformReleaseMeanZFR; Probe1.platformReleaseMeanZFR; Probe2.platformReleaseMeanZFR];
    allProbes.allUnitsPlatformCommandMeanZFR = [Probe0.platformCommandMeanZFR; Probe1.platformCommandMeanZFR; Probe2.platformCommandMeanZFR];
    allProbes.allUnitsColdAirOnMeanZFR = [Probe0.coldAirOnMeanZFR; Probe1.coldAirOnMeanZFR; Probe2.coldAirOnMeanZFR];
    allProbes.allUnitsMeanZFR_shuf1 = [Probe0.allUnitsMeanZFR_shuf1; Probe1.allUnitsMeanZFR_shuf1; Probe2.allUnitsMeanZFR_shuf1];
    allProbes.allUnitsMeanZFR_shuf2 = [Probe0.allUnitsMeanZFR_shuf2; Probe1.allUnitsMeanZFR_shuf2; Probe2.allUnitsMeanZFR_shuf2];
    allProbes.allUnitsAtCueMeanZFR_shuf1 = [Probe0.allUnitsAtCueMeanZFR_shuf1; Probe1.allUnitsAtCueMeanZFR_shuf1; Probe2.allUnitsAtCueMeanZFR_shuf1];
    allProbes.allUnitsAtCueMeanZFR_shuf2 = [Probe0.allUnitsAtCueMeanZFR_shuf2; Probe1.allUnitsAtCueMeanZFR_shuf2; Probe2.allUnitsAtCueMeanZFR_shuf2];
    allProbes.numProbeTypeUnits = [sum(Probe0.numGoodUnits); sum(Probe1.numGoodUnits); sum(Probe2.numGoodUnits)]; % for separating units across probe types
    allProbes.anatHiers = [Probe0.anatHiers Probe1.anatHiers Probe2.anatHiers]; % all anatomical information for all units
    allProbes.stDepths = [Probe0.stDepths; Probe1.stDepths; Probe2.stDepths]; 

%% plot raster map sort of all probes (***ASSUMES only 1 mouse/session loaded)
%     plotEphysClassicRasterByAnatomy(tsData_test, idxData_test, allStCellsAllProbes, allProbes.anatHiers, movementIndex, frameTimes) 
%     plotEphysTrialRaster(tsData_test, idxData_test, numData_test, shankData_test.stCell, shankData_test.anatHiers) %Fig 5 trial rasters

%% tuning plots
%     sortIdx = rastermapSort(allProbes.allUnitsMeanZFR'); % rastermap in large window around jump to capture entire trial structure
%     sortIdx = rastermapSort(allProbes.allUnitsDiffMeanZFR'); %rastermap to cluster by cov similarity of diff activity
    sortIdxAtCue = sortEphysZMeansAroundEvent(allProbes.allUnitsDiffAtCueMeanZFR, timescaleSeg, [-6 -5 0 1]);
    sortIdxAtVhex = sortEphysZMeansAroundEvent(allProbes.allUnitsDiffAtVhexMeanZFR, timescaleSeg, [-6 -5 -0.5 0]);
    sortIdxAtIti = sortEphysZMeansAroundEvent(allProbes.allUnitsDiffAtItiMeanZFR, timescaleSeg, [-1 0 0 1]);

    pcaWindowVhex = [-1.5 0.0]; %seconds around extension to calculate PCs
    pcaWindowCue = [0.0 1.5]; %seconds around CUE to calculate PCs; match window length for VHEx
    plotWindow = [-2 2]; %seconds around extension to plot trajectory 
    plotSamps = [find(plotWindow(1)>=timescaleSeg,1,'last') find(plotWindow(2)>=timescaleSeg,1,'last')];
    timescalePlot = timescaleSeg(plotSamps(1):plotSamps(2));
    itiWindow = [-6 -5]; %for sort compare to ITI

    plotBasicNeuronTuningGridPerProbe(allProbes, timescaleSeg, plotWindow, itiWindow) % Fig 5c
    plotAvoidAtCueNeuronTuningGrid(allProbes.allUnitsDiffAtCueMeanZFR, allProbes.allUnitsDiffAtVhexMeanZFR, allProbes.allUnitsPlatformExtendMeanZFR, allProbes.allUnitsColdAirOnMeanZFR, allProbes.allUnitsAvoidMeanZFR, allProbes.allUnitsReactMeanZFR, allProbes.anatHiers, timescaleSeg, plotWindow, sortIdxAtCue) 
    plotAvoidAtVhexNeuronTuningGrid(allProbes.allUnitsDiffAtVhexMeanZFR, allProbes.allUnitsDiffAtCueMeanZFR, allProbes.allUnitsPlatformExtendMeanZFR, allProbes.allUnitsColdAirOnMeanZFR, allProbes.allUnitsAvoidMeanZFR, allProbes.allUnitsReactMeanZFR, allProbes.anatHiers, timescaleSeg, plotWindow, sortIdxAtVhex) 

    % find neurons across areas that are up/down modulated in avoid vs react trials:
    plotAvoidVersusReactNeuronModulations(Probe0, Probe1, Probe2, timescaleSeg)

    % compare all VHEx to avoid & react
    plotEphysMeanZscoreFR3Events(allProbes.allUnitsDiffAtCueMeanZFR(sortIdxAtCue,plotSamps(1):plotSamps(2)), allProbes.allUnitsDiffAtVhexMeanZFR(sortIdxAtVhex,plotSamps(1):plotSamps(2)), allProbes.allUnitsDiffAtItiMeanZFR(sortIdxAtIti,plotSamps(1):plotSamps(2)), timescalePlot)

%% mean activity plots and "avoid/persistent neurons" and tuning pie charts
    % plot mean activity from struct
    plotMeanOverUnitsFromStruct(dataToDecodeWithNull, false)

    % anatomy plots for each probe
    plotUnitAnatomyPerProbe(Probe2, 'Probe2')

    % find "avoid neurons" by comparing avoid-react response to react only + plot tuning pie charts
    probesSub = {'Probe0','Probe1','Probe0','Probe2'}; %matched pairs with below
    anatsSub =  {'CTX',   'CTX',   'BS',    'HB'};
    labelsSub = {'M1',    'PFC',   'TH/HY', 'HB'}; 
    lineCmaps = linspecer(4);

    for probeAnatIdx = 1:length(probesSub)
        currData = dataToDecodeWithNull.(probesSub{probeAnatIdx}).(anatsSub{probeAnatIdx});
        currLabel = labelsSub{probeAnatIdx};
        currCmap = lineCmaps(probeAnatIdx,:);
        [specificAvoidIndeces.(probesSub{probeAnatIdx}).(anatsSub{probeAnatIdx}), nonSpecificAvoidIndeces.(probesSub{probeAnatIdx}).(anatsSub{probeAnatIdx})]...
            = findAvoidAtCueNeurons(currData, dataToDecodeWithNull.timescaleSeg, timescaleSeg, currLabel, currCmap);
%         [specificAvoidIndeces.(probesSub{probeAnatIdx}).(anatsSub{probeAnatIdx}), nonSpecificAvoidIndeces.(probesSub{probeAnatIdx}).(anatsSub{probeAnatIdx})]...
%             = findAvoidAtVhexNeurons(currData, timescaleSeg, currLabel);
    end

%% PCA plots from subspace
plotAvoidVsReactPcTrajectoriesFromSubspace(allDataSubspace,false); %example trajectories

plotAvoidVsReactPcDistancesFromSubspace(allDataSubspace,false); %control without null removed
plotAvoidVsReactPcDistancesFromSubspace(allDataSubspace,true); %with null removal

% plotAvoidVsReactPcDistancesWithOptoFromSubspace(allDataSubspace,false); %control without null removed
% plotAvoidVsReactPcDistancesWithOptoFromSubspace(allDataSubspace,true); %with null removal

%% PCA plots from trial data structure
plotAvoidVsReactPcDistancesFromStruct(dataToDecodeWithNull,false); %control without null removed
plotAvoidVsReactPcDistancesFromStruct(dataToDecodeWithNull,true); %with null removal

% plotAvoidVsReactPcDistancesWithOptoFromStruct(dataToDecodeWithNull,false); %control without null removed
% plotAvoidVsReactPcDistancesWithOptoFromStruct(dataToDecodeWithNull,true); %with null removal

end %end if ephysData

%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% loading and basic plotting functions
function cmap = blueOrange
    % colormap that linearly varies between triplet1 to white (=0) and triplet2
    % vmin & vmax are value limits for map
    triplet1 = [0 0.5 0.75]; %blue
    triplet2 = [0.9 0.7 0.1]; %orange
    numPoints = 256;
    R = [linspace(triplet1(1), 1, numPoints/2) linspace(1, triplet2(1), numPoints/2)]; 
    G = [linspace(triplet1(2), 1, numPoints/2) linspace(1, triplet2(2), numPoints/2)]; 
    B = [linspace(triplet1(3), 1, numPoints/2) linspace(1, triplet2(3), numPoints/2)]; 
    cmap = [R; G; B]';
end

function cmap = greenMag
    % colormap that linearly varies between triplet1 to white (=0) and triplet2
    % vmin & vmax are value limits for map
    triplet1 = [0.4660 0.6740 0.1880]; %green
    triplet2 = [0.7 0 0.7]; %magenta
    numPoints = 256;
    R = [linspace(triplet1(1), 1, numPoints/2) linspace(1, triplet2(1), numPoints/2)]; 
    G = [linspace(triplet1(2), 1, numPoints/2) linspace(1, triplet2(2), numPoints/2)]; 
    B = [linspace(triplet1(3), 1, numPoints/2) linspace(1, triplet2(3), numPoints/2)]; 
    cmap = [R; G; B]';
end

function plotBasicVhexDataNorm(tsData, idxData) %#ok<*DEFNU>
    %plot basic behavior data   
    jumpTimes = tsData.timescale(idxData.jumpIndeces);
    figure; 
    hold on;
    plot(tsData.timescale, tsData.jumpDistance ./ max(abs(tsData.jumpDistance)), 'b-', 'LineWidth', 2)
%     plot(tsData.timescale(1:end-1), tsData.jumpVel ./ max(abs(tsData.jumpVel)), 'b--', 'LineWidth', 2)
%     plot(tsData.timescale(1:end-2), tsData.jumpAccel ./ max(abs(tsData.jumpAccel)), 'g-', 'LineWidth', 2)
    plot(tsData.timescale, tsData.coldAirSolenoid, 'k-', 'LineWidth', 2)
    plot(tsData.timescale, tsData.beeps, 'r-', 'LineWidth', 2)
    plot(tsData.timescale, tsData.optoVoltage./5, 'c-', 'LineWidth', 2) %normalize by max 5V
    plot(jumpTimes, (ones(length(jumpTimes),1).*0.8), 'b*');
%     legend('distance', 'velocity', 'acceleration', 'cold air', 'beep', 'opto')
    legend('distance', 'cold air', 'beep', 'opto')
    set(gcf,'color','w');
    hold off;
    xlabel('time (sec)')
    ylabel('all signals a.u.')
    axis tight;
    makepretty;
end

function plotBasicVhexData(tsData) %#ok<*DEFNU>

    %plot basic behavior data    
    figure; 
    hold on;
    plot(tsData.timescale, tsData.jumpDistance, 'b-', 'LineWidth', 2)
    plot(tsData.timescale(1:end-1), tsData.jumpVel, 'b--', 'LineWidth', 2)
%     plot(tsData.timescale(1:end-2), tsData.jumpAccel, 'g-', 'LineWidth', 2)
    %plot(tsData.timescale(1:end-3), tsData.jumpJerk, 'g--', 'LineWidth', 2)
    plot(tsData.timescale, tsData.coldAirSolenoid, 'k-', 'LineWidth', 2)
    plot(tsData.timescale, tsData.beeps, 'r-', 'LineWidth', 2)
    plot(tsData.timescale, tsData.optoVoltage./5, 'c-', 'LineWidth', 2) %normalize by max 5V
    legend('distance', 'velocity', 'cold air', 'beep', 'opto')
    set(gcf,'color','w');
    hold off;
    xlabel('time (sec)')
    ylabel('cm, cm/sec, a.u.')
    axis tight;
    makepretty;
end

function zeroXDottedLine
% for current axes, draw a dotted line at X=0 extending to yLimits
    x = [0 0]; %dotted line at event = time 0
    yLims = get(gca, 'YLim');
    y = [yLims(1) yLims(2)];
    plot(x, y, 'k--', 'LineWidth', 0.1);
end

function zeroYDottedLine
% for current axes, draw a dotted line at Y=0 extending to xLimits
    y = [0 0]; %dotted line at event = time 0
    xLims = get(gca, 'XLim');
    x = [xLims(1) xLims(2)];
    plot(x, y, 'k--', 'LineWidth', 0.1);
end

function drawMouseSeparatorsCorr(trialsPerMouse)
%if multiple mice on correlation plot, draw gridlines to separate
    trialsPerMouseSeparators = cumsum(trialsPerMouse);
    if length(trialsPerMouse)>1 
        xLims = get(gca, 'XLim');
        x = [xLims(1) xLims(2)];
        yLims = get(gca, 'YLim');
        y = [yLims(1) yLims(2)];
        for m = 1:length(trialsPerMouse)-1
            mouseSeparator = [trialsPerMouseSeparators(m) trialsPerMouseSeparators(m)];
            plot(mouseSeparator,y,'k-', 'LineWidth', 2)
            plot(x,mouseSeparator,'k-', 'LineWidth', 2)
        end
    end
end

function drawMouseSeparatorsHeatmap(trialsPerMouse)
%if multiple mice on correlation plot, draw gridlines to separate
    trialsPerMouseSeparators = cumsum(trialsPerMouse);
    if length(trialsPerMouse)>1 
        xLims = get(gca, 'XLim');
        x = [xLims(1) xLims(2)];
        yLims = get(gca, 'YLim');
        y = [yLims(1) yLims(2)];
        for m = 1:length(trialsPerMouse)-1
            mouseSeparator = [trialsPerMouseSeparators(m) trialsPerMouseSeparators(m)] + 0.5; %add 0.5 so between heatmap rows
            plot(x,mouseSeparator,'k-', 'LineWidth', 2)
        end
    end
end

function drawMouseSeparatorsPlot(trialsPerMouse)
%if multiple mice on correlation plot, draw gridlines to separate
    trialsPerMouseSeparators = cumsum(trialsPerMouse);
    if length(trialsPerMouse)>1 
        xLims = get(gca, 'XLim');
        x = [xLims(1) xLims(2)];
        yLims = get(gca, 'YLim');
        y = [yLims(1) yLims(2)];
        for m = 1:length(trialsPerMouse)-1
            mouseSeparator = [trialsPerMouseSeparators(m) trialsPerMouseSeparators(m)];
            plot(mouseSeparator, y,'r-', 'LineWidth', 2)
        end
    end
end

function drawProbeSeparatorsPlot(unitsPerProbeType)
%if multiple probe types, draw gridlines to separate
    probeSeparators = cumsum(unitsPerProbeType);
    if length(unitsPerProbeType)>1 
        yLims = get(gca, 'YLim');
        for m = 1:length(unitsPerProbeType)-1
            probeSeparator = [probeSeparators(m) probeSeparators(m)];
            plot(probeSeparator, yLims,'k--', 'LineWidth', 2)
        end
    end
end

%% split timeseries into segments functions
function [avoidTrials, reactTrials, timescaleSeg] = getAvoidVsReactJumpSegments(dataToSplit, dataTimescale, behaviorTimescale, idxData, secBefore, secAfter)
% use index data to split time series data into segments by trial type, for a single session
% pass in appropriate dataToSplit and its corresponding timescale (since dataToSplit can have different clocks)
% also always pass in behavior timescale since idxData corresponds to it

    sampTime = dataTimescale(2) - dataTimescale(1);
    sampRate = 1/sampTime;
    numSampAfter = round(sampRate*secAfter);
    numSampBefore = round(sampRate*secBefore);
    timescaleSeg = (-numSampBefore*sampTime):sampTime:(numSampAfter*sampTime - sampTime);

    %break up data into avoid vs. react trials (non-opto):
    avoidTrials = zeros(length(find(idxData.jumpAntIndeces & ~idxData.jumpOptoIndeces)), numSampBefore + numSampAfter + 1, size(dataToSplit,2));
    reactTrials = zeros(length(find(~idxData.jumpAntIndeces & ~idxData.jumpOptoIndeces)), numSampBefore + numSampAfter + 1, size(dataToSplit,2));
    avoidCurrentTrial = 1;
    reactCurrentTrial = 1;
    for i = 1:length(idxData.jumpIndeces) %loop over all extensions for the session
        if ~idxData.jumpOptoIndeces(i) %only non-opto trials for now
            jumpTime = behaviorTimescale(idxData.jumpIndeces(i));
            jumpDataSamp = find(dataTimescale > jumpTime, 1, 'first'); %now switch to data timescale, but matches calibrated time in seconds
            if idxData.jumpAntIndeces(i)
                avoidTrials(avoidCurrentTrial,:,:) = dataToSplit(jumpDataSamp-numSampBefore:jumpDataSamp+numSampAfter,:);
                avoidCurrentTrial = avoidCurrentTrial + 1;
            else
                reactTrials(reactCurrentTrial,:,:) = dataToSplit(jumpDataSamp-numSampBefore:jumpDataSamp+numSampAfter,:);
                reactCurrentTrial = reactCurrentTrial + 1;
            end
        end
    end

end

function [avoidTrials, reactTrials, timescaleSeg] = getAvoidVsReactJumpConcat(dataToSplit, dataTimescale, behaviorTimescale, idxData, secBefore, secAfter)
% use index data to split time series data into segments by trial type, for a single session
% like segment function above, but concatenate segments across time to get matrix for PCA
% pass in appropriate dataToSplit and its corresponding timescale (since dataToSplit can have different clocks)
% also always pass in behavior timescale since idxData corresponds to it

    sampTime = dataTimescale(2) - dataTimescale(1);
    sampRate = 1/sampTime;
    numSampAfter = round(sampRate*secAfter);
    numSampBefore = round(sampRate*secBefore);
    segLen = numSampBefore + numSampAfter + 1;
    timescaleSeg = (-numSampBefore*sampTime):sampTime:(numSampAfter*sampTime - sampTime);

    %break up data into avoid vs. react trials (non-opto):
    avoidTrials = zeros(size(dataToSplit,2), segLen*length(find(idxData.jumpAntIndeces & ~idxData.jumpOptoIndeces)));
    reactTrials = zeros(size(dataToSplit,2), segLen*length(find(~idxData.jumpAntIndeces & ~idxData.jumpOptoIndeces)));
    avoidCurrentTrial = 1;
    reactCurrentTrial = 1;
    for i = 1:length(idxData.jumpIndeces) %loop over all extensions for the session
        if ~idxData.jumpOptoIndeces(i) %only non-opto trials for now
            jumpTime = behaviorTimescale(idxData.jumpIndeces(i));
            jumpDataSamp = find(dataTimescale > jumpTime, 1, 'first'); %now switch to data timescale, but matches calibrated time in seconds
            currSegData = dataToSplit(jumpDataSamp-numSampBefore:jumpDataSamp+numSampAfter,:);
            if idxData.jumpAntIndeces(i)
                startSamp = (avoidCurrentTrial-1)*segLen + 1;
                endSamp = avoidCurrentTrial*segLen;
                avoidTrials(:,startSamp:endSamp) = currSegData';%concatenate across time
                avoidCurrentTrial = avoidCurrentTrial + 1;
            else
                startSamp = (reactCurrentTrial-1)*segLen + 1;
                endSamp = reactCurrentTrial*segLen;
                reactTrials(:,startSamp:endSamp) = currSegData';%concatenate across time
                reactCurrentTrial = reactCurrentTrial + 1;
            end
        end
    end

end

function [optoTrials, notOptoTrials, timescaleSeg, numNotOptoTrials] = getOptoVsNotJumpConcat(dataToSplit, dataTimescale, behaviorTimescale, idxData, secBefore, secAfter)
% use index data to split time series data into segments by trial type, for a single session
% like segment function above, but concatenate segments across time to get matrix for PCA
% pass in appropriate dataToSplit and its corresponding timescale (since dataToSplit can have different clocks)
% also always pass in behavior timescale since idxData corresponds to it

    diffDataTimescale = diff(dataTimescale);
    sampTime = mean(diffDataTimescale);
    sampRate = 1/sampTime;
    numSampAfter = round(sampRate*secAfter);
    numSampBefore = round(sampRate*secBefore);
    segLen = numSampBefore + numSampAfter;
    timescaleSeg = (-numSampBefore*sampTime):sampTime:(numSampAfter*sampTime - sampTime);

    %break up data into opto vs. non-opto trials:
    optoTrials = zeros(size(dataToSplit,2), segLen*length(find(idxData.jumpOptoIndeces)));
    notOptoTrials = zeros(size(dataToSplit,2), segLen*length(find(~idxData.jumpOptoIndeces)));
    numNotOptoTrials = length(find(~idxData.jumpOptoIndeces));
    optoCurrentTrial = 1;
    notOptoCurrentTrial = 1;
    for i = 1:length(idxData.jumpIndeces) %loop over all extensions for the session
        if ~idxData.jumpOptoIndeces(i) %only non-opto trials for now
            jumpTime = behaviorTimescale(idxData.jumpIndeces(i));
            jumpDataSamp = find(dataTimescale > jumpTime, 1, 'first'); %now switch to data timescale, but matches calibrated time in seconds
            currSegData = dataToSplit(jumpDataSamp-numSampBefore:jumpDataSamp+numSampAfter-1,:);
            if idxData.jumpOptoIndeces(i)
                startSamp = (optoCurrentTrial-1)*segLen + 1;
                endSamp = optoCurrentTrial*segLen;
                optoTrials(:,startSamp:endSamp) = currSegData';%concatenate across time
                optoCurrentTrial = optoCurrentTrial + 1;
            else
                startSamp = (notOptoCurrentTrial-1)*segLen + 1;
                endSamp = notOptoCurrentTrial*segLen;
                notOptoTrials(:,startSamp:endSamp) = currSegData';%concatenate across time
                notOptoCurrentTrial = notOptoCurrentTrial + 1;
            end
        end
    end

end

function [optoTrials, notOptoTrials, timescaleSeg] = getOptoVsNotJumpSegments(dataToSplit, dataTimescale, behaviorTimescale, idxData, secBefore, secAfter)
% use index data to split time series data into segments by trial type, for a single session
% pass in appropriate dataToSplit and its corresponding timescale (since dataToSplit can have different clocks)
% also always pass in behavior timescale since idxData corresponds to it

    sampTime = dataTimescale(2) - dataTimescale(1);
    sampRate = 1/sampTime;
    numSampAfter = round(sampRate*secAfter);
    numSampBefore = round(sampRate*secBefore);
    timescaleSeg = (-numSampBefore*sampTime):sampTime:(numSampAfter*sampTime - sampTime);
    optoIndeces = idxData.jumpOptoIndeces;
    nonOptoIndeces = ~idxData.jumpOptoIndeces & ~idxData.jumpAntIndeces; %only compare apples/react

    %break up data into opto vs. non-opto trials:
    optoTrials = zeros(length(find(optoIndeces)), numSampBefore + numSampAfter + 1, size(dataToSplit,2));
    notOptoTrials = zeros(length(find(nonOptoIndeces)), numSampBefore + numSampAfter + 1, size(dataToSplit,2));
    optoCurrentTrial = 1;
    nonOptoCurrentTrial = 1;
    for i = 1:length(idxData.jumpIndeces) %loop over all extensions for the session
        jumpTime = behaviorTimescale(idxData.jumpIndeces(i));
        jumpDataSamp = find(dataTimescale > jumpTime, 1, 'first'); %now switch to data timescale, but matches calibrated time in seconds
        if optoIndeces(i)
            optoTrials(optoCurrentTrial,:,:) = dataToSplit(jumpDataSamp-numSampBefore:jumpDataSamp+numSampAfter,:);
            optoCurrentTrial = optoCurrentTrial + 1;
        elseif nonOptoIndeces(i)
            notOptoTrials(nonOptoCurrentTrial,:,:) = dataToSplit(jumpDataSamp-numSampBefore:jumpDataSamp+numSampAfter,:);
            nonOptoCurrentTrial = nonOptoCurrentTrial + 1;
        end
    end
end

function [trialCategories, shuffleTrialCategories, allTrialData, timescaleSeg] = getAllClassifierSegments(dataToSplit, dataTimescale, eventTimes1, eventTimes2, secBefore, secAfter, minTrials)
% function to get trial categories (e.g. avoid vs react) and data that can predict it, for feeding to trial type classifier
% choose a random number of largerEventTimes trials to match the number of smallerEventTimes trials, and create categorical outcome + shuffle control

    sampTime = dataTimescale(2) - dataTimescale(1);
    sampRate = 1/sampTime;
    numSampAfter = round(sampRate*secAfter);
    numSampBefore = round(sampRate*secBefore);
    timescaleSeg = (-numSampBefore*sampTime):sampTime:(numSampAfter*sampTime - sampTime);

    trialIdx = [];
    allTrialData = {};

    numShufEvents = minTrials; %use min of [#avoids /react /opto trials] to pick from larger distribution..this matches ephys analysis
    if length(eventTimes1) < length(eventTimes2)
        largerEventTimes = eventTimes2;
        smallerEventTimes = eventTimes1;
        eventTimes1Smaller = true;
    else
        largerEventTimes = eventTimes1;
        smallerEventTimes = eventTimes2;
        eventTimes1Smaller = false;
    end

    s = RandStream('dsfmt19937','Seed','shuffle');
    randIdx = randperm(s, length(largerEventTimes), numShufEvents); %draw randomly from largerEventTimes
    largerEventTimesMatched = largerEventTimes(randIdx);

    if smallerEventTimes > minTrials %need to make sure #opto trials isn't less - in this case we subsample both types of trials
        s = RandStream('dsfmt19937','Seed','shuffle');
        randIdx = randperm(s, length(smallerEventTimes), numShufEvents); %draw randomly from largerEventTimes
        smallerEventTimes = smallerEventTimes(randIdx);
    end

    for i = 1:numShufEvents
        eventSampS = find(dataTimescale > smallerEventTimes(i), 1, 'first');
        eventSampL = find(dataTimescale > largerEventTimesMatched(i), 1, 'first');
        currentDataS = dataToSplit(eventSampS-numSampBefore:eventSampS+numSampAfter-1,:);
        currentDataL = dataToSplit(eventSampL-numSampBefore:eventSampL+numSampAfter-1,:);
        allTrialData{end+1} = permute(currentDataS, [2 1]); %put units/keypoints in 1st dimension for MATLAB LSTM decoder; interleave eventTimes1 & eventTimes2
        allTrialData{end+1} = permute(currentDataL, [2 1]);
        trialIdx = [trialIdx eventTimes1Smaller]; %construct avoid/react indeces for only non-opto extensions
        trialIdx = [trialIdx ~eventTimes1Smaller];
    end

    trialCategories = categorical(trialIdx);
    s = RandStream('dsfmt19937','Seed','shuffle');
    randIdx = randperm(s, length(trialIdx), length(trialIdx)); %draw randomly from mix of eventTimes1 & eventTimes2
    shuffleTrialCategories = categorical(trialIdx(randIdx));
end

function [concatDataToClassify] = getAllClassifierSegmentsSubspaceConcat(subspaceData, probeNames, anatNames)
% function to (1) unpack per-mouse classifier cell array and combine across probes, thus preserving existing trial sampling procedure
% (2) re-pack data into probe structure, but using subspace-projected ephys predictor data 
% output is one cell containing all data (classifierTrialCategories, classifierTrialCategoriesShuf, classifierTrialData), with interleaved avoid & react trials for each sampling repeat
% trials are aligned at platformCommand cue

    numPCsToUse = 10; %number of PCs for decoder predictor data
    numSamplingRepeats = 10;

    for probeAnatIdx = 1:length(probeNames) %main loop over probe structures; ex. {'Probe1'} 

        thisAnatAvoidData = subspaceData.(probeNames{probeAnatIdx}).(anatNames{probeAnatIdx}).allUnitsTrialsAvoidAtCue; %this is cell with numMice, each having numSamplingRepeats
        thisAnatReactData = subspaceData.(probeNames{probeAnatIdx}).(anatNames{probeAnatIdx}).allUnitsTrialsReactAtCue;
        
        for rpt = 1:numSamplingRepeats

            thisRptAllTrialCategories = []; %empty arrays to concatenate trials across mice, independently per repeat
            thisRptAllTrialCategoriesShuf = [];
            thisRptAllTrialData = [];
            trialIdx = [];

            for i = 1:length(thisAnatAvoidData) %loop over numMice on inside loop to concatenate trials across mice

                thisMouseAvoidData = thisAnatAvoidData{i};
                thisMouseReactData = thisAnatReactData{i};
                thisRptAvoidData = thisMouseAvoidData{rpt}; %now we are unpacked to #trials cells with actual data matrices [numPCs x time]
                thisRptReactData = thisMouseReactData{rpt};
                numTrials = length(thisRptAvoidData);

                for j = 1:numTrials
                    thisTrialSubspaceAvoidData = thisRptAvoidData{j}; %thisRpt is [numPCs x time], units/PCs/keypoints in 1st dimension for MATLAB LSTM decoder
                    thisTrialSubspaceReactData = thisRptReactData{j};

                    thisRptAllTrialData = [thisRptAllTrialData {thisTrialSubspaceAvoidData} {thisTrialSubspaceReactData}]; % interleaved cells for avoid and react trials
                    thisRptAllTrialCategories = [thisRptAllTrialCategories {categorical(1)} {categorical(0)}];
                    trialIdx = [trialIdx 1 0];% keep an array as well to make shuffling easier below

                end %end trials
            end %end numMice

            % now make shuffled trial categories across all mice:
            s = RandStream('dsfmt19937','Seed','shuffle');
            randIdx = randperm(s, length(trialIdx), length(trialIdx)); %draw randomly from mix of avoid & react
            thisRptAllTrialCategoriesShuf = num2cell(categorical(trialIdx(randIdx)));

            % now re-pack all into single cell since dowstream functions expect cell input:
            concatDataToClassify.(probeNames{probeAnatIdx}).(anatNames{probeAnatIdx}).classifierTrialCategories{rpt} = thisRptAllTrialCategories;
            concatDataToClassify.(probeNames{probeAnatIdx}).(anatNames{probeAnatIdx}).classifierTrialCategoriesShuf{rpt} = thisRptAllTrialCategoriesShuf;
            concatDataToClassify.(probeNames{probeAnatIdx}).(anatNames{probeAnatIdx}).classifierTrialData{rpt} = thisRptAllTrialData;

        end %end rpt
    end %end main current probe loop

end

function [dataToDecodePlatform] = getPlatformDecoderConcat(dataToDecode, probeNames, anatNames, timescaleSeg, pcaWindow)
% function to (1) unpack per-mouse decoding cell array and combine across probes, thus preserving existing trial sampling procedure
% (2) compute shared subspace PCs and orthogonalized projections for each mouse/trial
% (3) re-pack data into probe structure, but concatenated for support vector regression

    numDecoderSampRepeats = 10; % num sampling repeats, should match definition above
    unitsThreshold = 20; %threshold below which decoding is skipped for particular anatomical region on particular probe; perhaps make equal to kinematics dimensionality, but at least as big as numPcsToUse
    usePcContrast = 0;
    alignPCs = 1;
    numPcsToUse = 10;
    smoothingSamples = 5; 
    pcaSamps = [find(pcaWindow(1)>=timescaleSeg,1,'last') find(pcaWindow(2)>=timescaleSeg,1,'last')];

    for probeIdx = probeNames %main loop over probe structures; ex. {'Probe1'} 
        if ~isempty(anatNames)  % select anatomical subsets of probe (ex. only cortex electrodes)
            for anatIdx = anatNames
                thisAnatNumGoodUnits = dataToDecode.(probeIdx{1}).(anatIdx{1}).numGoodUnits;
                if thisAnatNumGoodUnits<=unitsThreshold
                    continue %skip this anat if not enough units to give ~10 PCs to use after nullspace projection
                end
                numMice = length(thisAnatNumGoodUnits); %assume 1 session/probe per mouse
                thisAnatAllUnitTrialMeans = dataToDecode.(probeIdx{1}).(anatIdx{1}).allUnitsMeanZFR;
                thisAnatAllUnitAvoidMeans = squeeze(mean(dataToDecode.(probeIdx{1}).(anatIdx{1}).allUnitsAvoidMeanZFR,1)); %take mean across sampling repeats before PCA for avoid/react trials
                thisAnatAllUnitReactMeans = squeeze(mean(dataToDecode.(probeIdx{1}).(anatIdx{1}).allUnitsReactMeanZFR,1));

                if ~isempty(thisAnatNumGoodUnits) %check for any units found; data will only be saved by getProbeSubset() if above thresh there
                    % first find cold air nullspace from time period around cold air on when latency in long enough that movement doesn't interfere
                    % then use this nullspace to minimize cold air contribution by removing it from the PCs that will later be used to decode platform/extension distance
                    thisAnatAllUnitColdOnMeans = dataToDecode.(probeIdx{1}).(anatIdx{1}).allUnitsColdAirOnMeanZFR;
%                     thisAnatAllUnitPlatformReleaseMeans = dataToDecode.(probeIdx{1}).(anatIdx{1}).allUnitsPlatformReleaseMeanZFR;
                    nullspaceWindow = [-0.5 0.2]; %for cold air
%                     nullspaceWindow = [-1.5 0.1]; %for VHEx
                    nullSamps = [find(nullspaceWindow(1)>=timescaleSeg,1,'last') find(nullspaceWindow(2)>=timescaleSeg,1,'last')];
                    % get aligned, re-orthogonalized PCs for each mouse/probe, after removing covariance patterns from nullspace (e.g. cold air on)
                    [pcSeg, varExplained] = getEphysPCsWithNullspace(thisAnatAllUnitTrialMeans, thisAnatAllUnitTrialMeans, thisAnatAllUnitTrialMeans, thisAnatAllUnitColdOnMeans, thisAnatNumGoodUnits, numMice, numPcsToUse, pcaSamps, nullSamps, usePcContrast); % remove cold air info
%                     [pcSeg, varExplained] = getEphysPCsWithNullspace(thisAnatAllUnitTrialMeans, thisAnatAllUnitTrialMeans, thisAnatAllUnitTrialMeans, thisAnatAllUnitPlatformReleaseMeans, thisAnatNumGoodUnits, numMice, numPcsToUse, pcaSamps, nullSamps); %remove cue info
%                     [pcSeg, varExplained] = getEphysPCsWithNullspace(thisAnatAllUnitTrialMeans, thisAnatAllUnitTrialMeans, thisAnatAllUnitTrialMeans, thisAnatAllUnitAvoidMeans, thisAnatNumGoodUnits, numMice, numPcsToUse, pcaSamps, nullSamps); %control remove avoid VHEx info, should leave only cold

                    thisAnatPlatDistZAvoid = dataToDecode.(probeIdx{1}).(anatIdx{1}).platDistZAvoid; % at the top level these arrays have a cell for each mouse, then below that a cell for each sampling repeat
                    thisAnatPlatDistZReact = dataToDecode.(probeIdx{1}).(anatIdx{1}).platDistZReact;
                    thisAnatPlatDistZOpto = dataToDecode.(probeIdx{1}).(anatIdx{1}).platDistZOpto;
                    thisAnatPlatDistZIti = dataToDecode.(probeIdx{1}).(anatIdx{1}).platDistZIti;
                    thisAnatAllUnitsTrialsAvoid = dataToDecode.(probeIdx{1}).(anatIdx{1}).allUnitsTrialsAvoid;
                    thisAnatAllUnitsTrialsReact = dataToDecode.(probeIdx{1}).(anatIdx{1}).allUnitsTrialsReact;
                    thisAnatAllUnitsTrialsOpto = dataToDecode.(probeIdx{1}).(anatIdx{1}).allUnitsTrialsOpto;
                    thisAnatAllUnitsTrialsIti = dataToDecode.(probeIdx{1}).(anatIdx{1}).allUnitsTrialsIti;

                    %now remove empty cells where numUnits threshold not reached in getProbeSubset():
                    thisAnatPlatDistZAvoid = thisAnatPlatDistZAvoid(~cellfun('isempty',thisAnatPlatDistZAvoid));
                    thisAnatPlatDistZReact = thisAnatPlatDistZReact(~cellfun('isempty',thisAnatPlatDistZReact));
                    thisAnatPlatDistZOpto = thisAnatPlatDistZOpto(~cellfun('isempty',thisAnatPlatDistZOpto));
                    thisAnatPlatDistZIti = thisAnatPlatDistZIti(~cellfun('isempty',thisAnatPlatDistZIti));
                    thisAnatAllUnitsTrialsAvoid = thisAnatAllUnitsTrialsAvoid(~cellfun('isempty',thisAnatAllUnitsTrialsAvoid));
                    thisAnatAllUnitsTrialsReact = thisAnatAllUnitsTrialsReact(~cellfun('isempty',thisAnatAllUnitsTrialsReact));
                    thisAnatAllUnitsTrialsOpto = thisAnatAllUnitsTrialsOpto(~cellfun('isempty',thisAnatAllUnitsTrialsOpto));
                    thisAnatAllUnitsTrialsIti = thisAnatAllUnitsTrialsIti(~cellfun('isempty',thisAnatAllUnitsTrialsIti));

                    numTimepoints = size(thisAnatPlatDistZAvoid{1}{1},2);
                    totalTrials = 0;
                    for i = 1:numMice %quick loop to figure out total # trials across mice to preallocate arrays
                        currentMouse = thisAnatPlatDistZAvoid{i}{1};
                        currentTrials = size(currentMouse,1);
                        totalTrials = totalTrials + currentTrials;
                    end

                    % initialize arrays to collect concatenated behavior & subspace projections for predictors, for each sampling repeat
                    platDistZAvoidConcat = NaN(numDecoderSampRepeats, totalTrials*numTimepoints); 
                    platDistZReactConcat = NaN(numDecoderSampRepeats, totalTrials*numTimepoints); 
                    platDistZOptoConcat = NaN(numDecoderSampRepeats, totalTrials*numTimepoints); 
                    platDistZItiConcat = NaN(numDecoderSampRepeats, totalTrials*numTimepoints); 
                    allUnitsTrialsAvoidConcat = NaN(numDecoderSampRepeats, totalTrials*numTimepoints, numPcsToUse); 
                    allUnitsTrialsReactConcat = NaN(numDecoderSampRepeats, totalTrials*numTimepoints, numPcsToUse); 
                    allUnitsTrialsOptoConcat = NaN(numDecoderSampRepeats, totalTrials*numTimepoints, numPcsToUse); 
                    allUnitsTrialsItiConcat = NaN(numDecoderSampRepeats, totalTrials*numTimepoints, numPcsToUse);    

                    endMouseIdx = 0;
                    for i = 1:numMice %loop over mice/probes now (those with enough units) to apply PC projection and concatenate output

                        thisMouseAvoidEphysData = thisAnatAllUnitsTrialsAvoid{i}; %unpack to get per-mouse cells that include all sample repeats
                        thisMouseReactEphysData = thisAnatAllUnitsTrialsReact{i};
                        thisMouseOptoEphysData = thisAnatAllUnitsTrialsOpto{i};
                        thisMouseItiEphysData = thisAnatAllUnitsTrialsIti{i};
                        thisMousePlatDistZAvoid = thisAnatPlatDistZAvoid{i};
                        thisMousePlatDistZReact = thisAnatPlatDistZReact{i};
                        thisMousePlatDistZOpto = thisAnatPlatDistZOpto{i};
                        thisMousePlatDistZIti = thisAnatPlatDistZIti{i};
                        
                        for decRpt = 1:numDecoderSampRepeats
                            startIdx = endMouseIdx + 1; %for each mouse loop, keep concatenating from last

                            thisRptAvoidEphysData = thisMouseAvoidEphysData{decRpt}; %unpack one more time to get each repeat cell - now this cell for each trial, each having units x timepoints matrix
                            thisRptReactEphysData = thisMouseReactEphysData{decRpt};
                            thisRptOptoEphysData = thisMouseOptoEphysData{decRpt};
                            thisRptItiEphysData = thisMouseItiEphysData{decRpt};
                            thisRptPlatDistZAvoid = thisMousePlatDistZAvoid{decRpt}; %unpack once more - this is now # trials x timepoints
                            thisRptPlatDistZReact = thisMousePlatDistZReact{decRpt};
                            thisRptPlatDistZOpto = thisMousePlatDistZOpto{decRpt};
                            thisRptPlatDistZIti = thisMousePlatDistZIti{decRpt};
                            
                            numTrials = size(thisRptPlatDistZAvoid,1);

                            for j = 1:numTrials
                                trialProjAvoid = pcSeg{i}'*thisRptAvoidEphysData{j}; %apply PCs to individual trials; output is #PCsToUse x #timepoints
                                trialProjAvoid = smoothdata(trialProjAvoid, 2, 'gaussian', smoothingSamples); %ex. smooth window 25 is 25*20ms PSTH bins = 500ms
                                trialProjReact = pcSeg{i}'*thisRptReactEphysData{j}; %
                                trialProjReact = smoothdata(trialProjReact, 2, 'gaussian', smoothingSamples); %
                                trialProjOpto = pcSeg{i}'*thisRptOptoEphysData{j}; %
                                trialProjOpto = smoothdata(trialProjOpto, 2, 'gaussian', smoothingSamples); %
                                trialProjIti = pcSeg{i}'*thisRptItiEphysData{j}; %
                                trialProjIti = smoothdata(trialProjIti, 2, 'gaussian', smoothingSamples); %
                                platAvoidTrial = thisRptPlatDistZAvoid(j,:);
                                platReactTrial = thisRptPlatDistZReact(j,:);
                                platOptoTrial = thisRptPlatDistZOpto(j,:);
                                platItiTrial = thisRptPlatDistZIti(j,:);

                                endIdx = startIdx + numTimepoints - 1;
                    
                                % concatenate all trials over time
                                platDistZAvoidConcat(decRpt,startIdx:endIdx) = platAvoidTrial;
                                platDistZReactConcat(decRpt,startIdx:endIdx) = platReactTrial;
                                platDistZOptoConcat(decRpt,startIdx:endIdx) = platOptoTrial; 
                                platDistZItiConcat(decRpt,startIdx:endIdx) = platItiTrial;
                                allUnitsTrialsAvoidConcat(decRpt,startIdx:endIdx,:) = trialProjAvoid';
                                allUnitsTrialsReactConcat(decRpt,startIdx:endIdx,:) = trialProjReact';
                                allUnitsTrialsOptoConcat(decRpt,startIdx:endIdx,:) = trialProjOpto';
                                allUnitsTrialsItiConcat(decRpt,startIdx:endIdx,:) = trialProjIti';

                                startIdx = endIdx + 1;
                            end %end trials

                        end %end sampling repeats
                        endMouseIdx = endIdx;
                    end %end numMice

                    % now re-pack all into single cell matrices:
                    dataToDecodePlatform.(probeIdx{1}).(anatIdx{1}).platDistZAvoid = platDistZAvoidConcat;
                    dataToDecodePlatform.(probeIdx{1}).(anatIdx{1}).platDistZReact = platDistZReactConcat;
                    dataToDecodePlatform.(probeIdx{1}).(anatIdx{1}).platDistZOpto = platDistZOptoConcat;
                    dataToDecodePlatform.(probeIdx{1}).(anatIdx{1}).platDistZIti = platDistZItiConcat;
                    dataToDecodePlatform.(probeIdx{1}).(anatIdx{1}).allUnitsTrialsAvoid = allUnitsTrialsAvoidConcat;
                    dataToDecodePlatform.(probeIdx{1}).(anatIdx{1}).allUnitsTrialsReact = allUnitsTrialsReactConcat;
                    dataToDecodePlatform.(probeIdx{1}).(anatIdx{1}).allUnitsTrialsOpto = allUnitsTrialsOptoConcat;
                    dataToDecodePlatform.(probeIdx{1}).(anatIdx{1}).allUnitsTrialsIti = allUnitsTrialsItiConcat;   
                    
                end % end thisAnatNumUnits check
            end % anatIdx
        end %~isempty(anatName)
    end %end main current probe loop
end

function [dataToDecodePlatform] = getPlatformDecoderSubspaceConcat(dataToDecode, subspaceData, probeNames, anatNames, removeNull)
% function to (1) unpack per-mouse decoding cell array and combine across probes, thus preserving existing trial sampling procedure
% (2) re-pack data into probe structure, but concatenated for support vector regression

    numPcsToUse = 10;
    numDecoderSampRepeats = 10; % num sampling repeats, should match definition above
    decodeWindow = [-2 2];
    timescaleFull = subspaceData.timescaleSeg;
    decodeSamps = [find(decodeWindow(1)>=timescaleFull,1,'last') find(decodeWindow(2)>=timescaleFull,1,'last')-1];
    timescaleSeg = timescaleFull(decodeSamps(1):decodeSamps(2)); %subsample in case decoder segments are longer than saved platform segments
    numTimepoints = size(timescaleSeg,2);
    numMice = size(subspaceData.(probeNames{1}).(anatNames{1}).allUnitsTrialsAvoid,2);

    for probeAnatIdx = 1:length(probeNames) %main loop over probe structures; ex. {'Probe1'} 
 
        thisAnatPlatDistZAvoid = dataToDecode.(probeNames{probeAnatIdx}).(anatNames{probeAnatIdx}).platDistZAvoid; % at the top level these arrays have a cell for each mouse, then below that a cell for each sampling repeat
        thisAnatPlatDistZReact = dataToDecode.(probeNames{probeAnatIdx}).(anatNames{probeAnatIdx}).platDistZReact;
        thisAnatPlatDistZOpto = dataToDecode.(probeNames{probeAnatIdx}).(anatNames{probeAnatIdx}).platDistZOpto;
        if removeNull
            thisAnatAllUnitsTrialsAvoid = subspaceData.(probeNames{probeAnatIdx}).(anatNames{probeAnatIdx}).allUnitsTrialsAvoidNull;
            thisAnatAllUnitsTrialsReact = subspaceData.(probeNames{probeAnatIdx}).(anatNames{probeAnatIdx}).allUnitsTrialsReactNull;
            thisAnatAllUnitsTrialsOpto = subspaceData.(probeNames{probeAnatIdx}).(anatNames{probeAnatIdx}).allUnitsTrialsOptoNull;
        else
            thisAnatAllUnitsTrialsAvoid = subspaceData.(probeNames{probeAnatIdx}).(anatNames{probeAnatIdx}).allUnitsTrialsAvoid;
            thisAnatAllUnitsTrialsReact = subspaceData.(probeNames{probeAnatIdx}).(anatNames{probeAnatIdx}).allUnitsTrialsReact;
            thisAnatAllUnitsTrialsOpto = subspaceData.(probeNames{probeAnatIdx}).(anatNames{probeAnatIdx}).allUnitsTrialsOpto;
        end

        totalTrials = 0;
        for i = 1:numMice %quick loop to figure out total # trials across mice to preallocate arrays
            currentMouse = thisAnatPlatDistZAvoid{i}{1};
            currentTrials = size(currentMouse,1);
            totalTrials = totalTrials + currentTrials;
        end

        % initialize arrays to collect concatenated behavior & subspace projections for predictors, for each sampling repeat
        platDistZAvoidConcat = NaN(numDecoderSampRepeats, totalTrials*numTimepoints); 
        platDistZReactConcat = NaN(numDecoderSampRepeats, totalTrials*numTimepoints); 
        platDistZOptoConcat = NaN(numDecoderSampRepeats, totalTrials*numTimepoints); 
        allUnitsTrialsAvoidConcat = NaN(numDecoderSampRepeats, totalTrials*numTimepoints, numPcsToUse); 
        allUnitsTrialsReactConcat = NaN(numDecoderSampRepeats, totalTrials*numTimepoints, numPcsToUse); 
        allUnitsTrialsOptoConcat = NaN(numDecoderSampRepeats, totalTrials*numTimepoints, numPcsToUse); 

        endMouseIdx = 0;
        for i = 1:numMice %loop over mice/probes now (those with enough units) to apply PC projection and concatenate output

            thisMouseAvoidEphysData = thisAnatAllUnitsTrialsAvoid{i}; %unpack to get per-mouse cells that include all sample repeats
            thisMouseReactEphysData = thisAnatAllUnitsTrialsReact{i};
            thisMouseOptoEphysData = thisAnatAllUnitsTrialsOpto{i};
            thisMousePlatDistZAvoid = thisAnatPlatDistZAvoid{i};
            thisMousePlatDistZReact = thisAnatPlatDistZReact{i};
            thisMousePlatDistZOpto = thisAnatPlatDistZOpto{i};
            
            for decRpt = 1:numDecoderSampRepeats
                startIdx = endMouseIdx + 1; %for each mouse loop, keep concatenating from last

                thisRptAvoidEphysData = thisMouseAvoidEphysData{decRpt}; %unpack one more time to get each repeat cell - now this cell for each trial, each having units x timepoints matrix
                thisRptReactEphysData = thisMouseReactEphysData{decRpt};
                thisRptOptoEphysData = thisMouseOptoEphysData{decRpt};
                thisRptPlatDistZAvoid = thisMousePlatDistZAvoid{decRpt}; %unpack once more - this is now # trials x timepoints
                thisRptPlatDistZReact = thisMousePlatDistZReact{decRpt};
                thisRptPlatDistZOpto = thisMousePlatDistZOpto{decRpt};
                
                numTrials = size(thisRptPlatDistZAvoid,1);

                for j = 1:numTrials

                    platAvoidTrial = thisRptPlatDistZAvoid(j,decodeSamps(1):decodeSamps(2));
                    platReactTrial = thisRptPlatDistZReact(j,decodeSamps(1):decodeSamps(2));
                    platOptoTrial = thisRptPlatDistZOpto(j,decodeSamps(1):decodeSamps(2));
                    ephysAvoidTrial = thisRptAvoidEphysData{j}(:,decodeSamps(1):decodeSamps(2))'; %put subspace dimensions in columns for decoder
                    ephysReactTrial = thisRptReactEphysData{j}(:,decodeSamps(1):decodeSamps(2))';
                    ephysOptoTrial = thisRptOptoEphysData{j}(:,decodeSamps(1):decodeSamps(2))';

                    endIdx = startIdx + numTimepoints - 1;
        
                    % concatenate all trials over time
                    platDistZAvoidConcat(decRpt,startIdx:endIdx) = platAvoidTrial;
                    platDistZReactConcat(decRpt,startIdx:endIdx) = platReactTrial;
                    platDistZOptoConcat(decRpt,startIdx:endIdx) = platOptoTrial; 
                    allUnitsTrialsAvoidConcat(decRpt,startIdx:endIdx,:) = ephysAvoidTrial;
                    allUnitsTrialsReactConcat(decRpt,startIdx:endIdx,:) = ephysReactTrial;
                    allUnitsTrialsOptoConcat(decRpt,startIdx:endIdx,:) = ephysOptoTrial;

                    startIdx = endIdx + 1;
                end %end trials

            end %end sampling repeats
            endMouseIdx = endIdx;
        end %end numMice

        % now re-pack all into single cell matrices:
        dataToDecodePlatform.(probeNames{probeAnatIdx}).(anatNames{probeAnatIdx}).platDistZAvoid = platDistZAvoidConcat;
        dataToDecodePlatform.(probeNames{probeAnatIdx}).(anatNames{probeAnatIdx}).platDistZReact = platDistZReactConcat;
        dataToDecodePlatform.(probeNames{probeAnatIdx}).(anatNames{probeAnatIdx}).platDistZOpto = platDistZOptoConcat;
        dataToDecodePlatform.(probeNames{probeAnatIdx}).(anatNames{probeAnatIdx}).allUnitsTrialsAvoid = allUnitsTrialsAvoidConcat;
        dataToDecodePlatform.(probeNames{probeAnatIdx}).(anatNames{probeAnatIdx}).allUnitsTrialsReact = allUnitsTrialsReactConcat;
        dataToDecodePlatform.(probeNames{probeAnatIdx}).(anatNames{probeAnatIdx}).allUnitsTrialsOpto = allUnitsTrialsOptoConcat;
            
    end %end main current probe/anat loop

end

function [dataToDecodePlatform] = getPlatformDecoderConcatPerMouse(dataToDecode, probeNames, anatNames, removeNull)
% function to (1) unpack per-mouse decoding cell array and combine across probes, thus preserving existing trial sampling procedure
% (2) re-pack data into probe structure, but concatenated for support vector regression (per mouse, so that using raw neural data is possible)

    numDecoderSampRepeats = 10; % num sampling repeats, should match definition above
    decodeWindow = [-2 2];
    timescaleFull = dataToDecode.timescaleSeg;
    decodeSamps = [find(decodeWindow(1)>=timescaleFull,1,'last') find(decodeWindow(2)>=timescaleFull,1,'last')];
    timescaleSeg = timescaleFull(decodeSamps(1):decodeSamps(2)); %subsample in case decoder segments are longer than saved platform segments
    numTimepoints = size(timescaleSeg,2);
    numMice = size(dataToDecode.(probeNames{1}).(anatNames{1}).allUnitsTrialsAvoidMean,2);

    for probeAnatIdx = 1:length(probeNames) %main loop over probe structures; ex. {'Probe1'} 
 
        thisAnatPlatDistZAvoid = dataToDecode.(probeNames{probeAnatIdx}).(anatNames{probeAnatIdx}).platDistZAvoid; % at the top level these arrays have a cell for each mouse, then below that a cell for each sampling repeat
        thisAnatPlatDistZReact = dataToDecode.(probeNames{probeAnatIdx}).(anatNames{probeAnatIdx}).platDistZReact;
        thisAnatPlatDistZOpto = dataToDecode.(probeNames{probeAnatIdx}).(anatNames{probeAnatIdx}).platDistZOpto;
        if removeNull
            thisAnatAllUnitsTrialsAvoid = dataToDecode.(probeNames{probeAnatIdx}).(anatNames{probeAnatIdx}).allUnitsTrialsAvoidNull;
            thisAnatAllUnitsTrialsReact = dataToDecode.(probeNames{probeAnatIdx}).(anatNames{probeAnatIdx}).allUnitsTrialsReactNull;
            thisAnatAllUnitsTrialsOpto = dataToDecode.(probeNames{probeAnatIdx}).(anatNames{probeAnatIdx}).allUnitsTrialsOptoNull;
        else
            thisAnatAllUnitsTrialsAvoid = dataToDecode.(probeNames{probeAnatIdx}).(anatNames{probeAnatIdx}).allUnitsTrialsAvoid;
            thisAnatAllUnitsTrialsReact = dataToDecode.(probeNames{probeAnatIdx}).(anatNames{probeAnatIdx}).allUnitsTrialsReact;
            thisAnatAllUnitsTrialsOpto = dataToDecode.(probeNames{probeAnatIdx}).(anatNames{probeAnatIdx}).allUnitsTrialsOpto;
        end

        for i = 1:numMice %loop over mice/probes now (those with enough units) to apply PC projection and concatenate output

            thisMouseAvoidEphysData = thisAnatAllUnitsTrialsAvoid{i}; %unpack to get per-mouse cells that include all sample repeats
            thisMouseReactEphysData = thisAnatAllUnitsTrialsReact{i};
            thisMouseOptoEphysData = thisAnatAllUnitsTrialsOpto{i};
            thisMousePlatDistZAvoid = thisAnatPlatDistZAvoid{i};
            thisMousePlatDistZReact = thisAnatPlatDistZReact{i};
            thisMousePlatDistZOpto = thisAnatPlatDistZOpto{i};
            
            numTrials = size(thisMouseAvoidEphysData{1},2);
            numUnits = size(thisMouseAvoidEphysData{1}{1},1);

            % initialize arrays to collect concatenated behavior & neural data (predictors), for each mouse and sampling repeat
            platDistZAvoidConcat = NaN(numDecoderSampRepeats, numTrials*numTimepoints); 
            platDistZReactConcat = NaN(numDecoderSampRepeats, numTrials*numTimepoints); 
            platDistZOptoConcat = NaN(numDecoderSampRepeats, numTrials*numTimepoints); 
            allUnitsTrialsAvoidConcat = NaN(numDecoderSampRepeats, numTrials*numTimepoints, numUnits); 
            allUnitsTrialsReactConcat = NaN(numDecoderSampRepeats, numTrials*numTimepoints, numUnits); 
            allUnitsTrialsOptoConcat = NaN(numDecoderSampRepeats, numTrials*numTimepoints, numUnits); 

            for decRpt = 1:numDecoderSampRepeats
                startIdx = 1; %start over for each repeat

                thisRptAvoidEphysData = thisMouseAvoidEphysData{decRpt}; %unpack one more time to get each repeat cell - now this cell for each trial, each having units x timepoints matrix
                thisRptReactEphysData = thisMouseReactEphysData{decRpt};
                thisRptOptoEphysData = thisMouseOptoEphysData{decRpt};
                thisRptPlatDistZAvoid = thisMousePlatDistZAvoid{decRpt}; %unpack once more - this is now # trials x timepoints
                thisRptPlatDistZReact = thisMousePlatDistZReact{decRpt};
                thisRptPlatDistZOpto = thisMousePlatDistZOpto{decRpt};
                
                for j = 1:numTrials

                    platAvoidTrial = thisRptPlatDistZAvoid(j,:);
                    platReactTrial = thisRptPlatDistZReact(j,:);
                    platOptoTrial = thisRptPlatDistZOpto(j,:);
                    ephysAvoidTrial = thisRptAvoidEphysData{j}(:,decodeSamps(1):decodeSamps(2))'; %put unit dimensions in columns for decoder
                    ephysReactTrial = thisRptReactEphysData{j}(:,decodeSamps(1):decodeSamps(2))';
                    ephysOptoTrial = thisRptOptoEphysData{j}(:,decodeSamps(1):decodeSamps(2))';

                    endIdx = startIdx + numTimepoints - 1;
        
                    % concatenate all trials over time
                    platDistZAvoidConcat(decRpt,startIdx:endIdx) = platAvoidTrial;
                    platDistZReactConcat(decRpt,startIdx:endIdx) = platReactTrial;
                    platDistZOptoConcat(decRpt,startIdx:endIdx) = platOptoTrial; 
                    allUnitsTrialsAvoidConcat(decRpt,startIdx:endIdx,:) = ephysAvoidTrial;
                    allUnitsTrialsReactConcat(decRpt,startIdx:endIdx,:) = ephysReactTrial;
                    allUnitsTrialsOptoConcat(decRpt,startIdx:endIdx,:) = ephysOptoTrial;
                    
                    startIdx = endIdx + 1;
                end %end trials
            end %end sampling repeats

            % now re-pack all into single cell matrices per mouse:
            dataToDecodePlatform.(probeNames{probeAnatIdx}).(anatNames{probeAnatIdx}){i}.platDistZAvoid = platDistZAvoidConcat;
            dataToDecodePlatform.(probeNames{probeAnatIdx}).(anatNames{probeAnatIdx}){i}.platDistZReact = platDistZReactConcat;
            dataToDecodePlatform.(probeNames{probeAnatIdx}).(anatNames{probeAnatIdx}){i}.platDistZOpto = platDistZOptoConcat;
            dataToDecodePlatform.(probeNames{probeAnatIdx}).(anatNames{probeAnatIdx}){i}.allUnitsTrialsAvoid = allUnitsTrialsAvoidConcat;
            dataToDecodePlatform.(probeNames{probeAnatIdx}).(anatNames{probeAnatIdx}){i}.allUnitsTrialsReact = allUnitsTrialsReactConcat;
            dataToDecodePlatform.(probeNames{probeAnatIdx}).(anatNames{probeAnatIdx}){i}.allUnitsTrialsOpto = allUnitsTrialsOptoConcat;
        end %end numMice  
    end %end main current probe/anat loop

end

function [concatDataToDecode] = getAllClassifierSegmentsSubspaceConcatKin(dataToClassifyKinData, timescaleSeg, pcaWindow, kinMeans)
% function to (1) unpack per-mouse decoding cell array and combine across mice, thus preserving existing trial sampling procedure
% (2) compute shared subspace PCs and orthogonalized projections for each mouse/trial
% (3) re-pack data into cell structure, but only one cell containing all trials (classifierTrialCategories, classifierTrialCategoriesShuf, classifierTrialData)

    usePcContrast = 0;
    alignPCs = 1;
    smoothingSamples = 1; 
    numSampRepeats = 10;
    pcaSamps = [find(pcaWindow(1)>=timescaleSeg,1,'last') find(pcaWindow(2)>=timescaleSeg,1,'last')];
    numMice = length(dataToClassifyKinData.trialLabels{1});

    trialMeans = kinMeans.allKptsAtCueMeanZ_allMice; 
    avoidMeans = kinMeans.avoidAtCueMeanZ_allMice;
    reactMeans = kinMeans.reactAtCueMeanZ_allMice;
    numKpts = 36;

    % get aligned, re-orthogonalized PCs for each mouse
    [pcSeg, varExplained] = getKinPCs(trialMeans, avoidMeans, reactMeans, ones(numMice,1).*numKpts, numMice, pcaSamps, usePcContrast, alignPCs);

    trialCategories = dataToClassifyKinData.trialLabels; %this is cell array (numMice) with all categorical trial labels
    trialCategoriesShuf = dataToClassifyKinData.trialLabelsShuffled;
    trialData = dataToClassifyKinData.trialData;

    for rpt = 1:numSampRepeats %for each sampling repeat, get collects cells for each trial, each containing [PCs x time]

        allTrialCategories = []; 
        allTrialCategoriesShuf = [];
        allTrialData = [];

        for i = 1:numMice %main loop over mice 
            numTrials = size(trialCategories{i}{rpt},2); %for kinematic data rpts are nested within mouse cells
            thisMouseTrialData = trialData{i}{rpt};

            for j = 1:numTrials
                trialProj = pcSeg{i}'*thisMouseTrialData{j}; %apply PCs to individual trials
                trialProj = smoothdata(trialProj, 2, 'gaussian', smoothingSamples); %ex. smooth window 25 is 25*20ms PSTH bins = 500ms
                allTrialData = [allTrialData {trialProj}];
                allTrialCategories = [allTrialCategories {trialCategories{i}{rpt}(j)}];
                allTrialCategoriesShuf = [allTrialCategoriesShuf {trialCategoriesShuf{i}{rpt}(j)}];
            end
        end %end main current mouse loop
    

        % now re-pack all into single cell since downstream functions expect cell input:
        concatDataToDecode.trialLabels{rpt} = allTrialCategories;
        concatDataToDecode.trialLabelsShuffled{rpt} = allTrialCategoriesShuf;
        concatDataToDecode.trialData{rpt} = allTrialData;
    end
end

function [dataToDecodeKinPlatformData] = getPlatformDecoderConcatKin(platDistZAvoid_allMice, platDistZReact_allMice, platDistZOpto_allMice, platDistZIti_allMice, kinKptsTrialsAvoid_allMice, kinKptsTrialsReact_allMice, kinKptsTrialsOpto_allMice, kinKptsTrialsIti_allMice, timescaleSeg, kinMeans, pcaWindow)
% function unpack per-mouse decoding cell array and combine across mice, using full keypoints

    numMice = length(platDistZAvoid_allMice);

    platDistZAvoidConcat = []; 
    platDistZReactConcat = [];
    platDistZOptoConcat = []; 
    platDistZItiConcat = [];
    kinKptsTrialsAvoidConcat = [];
    kinKptsTrialsReactConcat = [];
    kinKptsTrialsOptoConcat = [];
    kinKptsTrialsItiConcat = [];

    for i = 1:numMice %main loop over mice 
        numTrials = size(platDistZAvoid_allMice{i},1);
        thisMouseAvoidPlatData = platDistZAvoid_allMice{i};
        thisMouseAvoidKptData = kinKptsTrialsAvoid_allMice{i};
        thisMouseReactPlatData = platDistZReact_allMice{i};
        thisMouseReactKptData = kinKptsTrialsReact_allMice{i};
        thisMouseOptoPlatData = platDistZOpto_allMice{i};
        thisMouseOptoKptData = kinKptsTrialsOpto_allMice{i};
        thisMouseItiPlatData = platDistZIti_allMice{i};
        thisMouseItiKptData = kinKptsTrialsIti_allMice{i};

        for j = 1:numTrials
            trialProjAvoid = thisMouseAvoidKptData{j}; 
            trialProjReact = thisMouseReactKptData{j}; 
            trialProjOpto = thisMouseOptoKptData{j};
            trialProjIti = thisMouseItiKptData{j}; 
            platAvoidTrial = thisMouseAvoidPlatData(j,:);
            platReactTrial = thisMouseReactPlatData(j,:);
            platOptoTrial = thisMouseOptoPlatData(j,:);
            platItiTrial = thisMouseItiPlatData(j,:);

            % concatenate all trials over time
            platDistZAvoidConcat = [platDistZAvoidConcat; platAvoidTrial']; %make column vector of concatenated platform distances
            platDistZReactConcat = [platDistZReactConcat; platReactTrial'];
            platDistZOptoConcat = [platDistZOptoConcat; platOptoTrial']; 
            platDistZItiConcat = [platDistZItiConcat; platItiTrial'];
            kinKptsTrialsAvoidConcat = [kinKptsTrialsAvoidConcat; trialProjAvoid'];
            kinKptsTrialsReactConcat = [kinKptsTrialsReactConcat; trialProjReact'];
            kinKptsTrialsOptoConcat = [kinKptsTrialsOptoConcat; trialProjOpto'];
            kinKptsTrialsItiConcat = [kinKptsTrialsItiConcat; trialProjIti'];
        end
    end %end main current mouse loop

    % now re-pack all into single cell matrices:
    dataToDecodeKinPlatformData.platDistZAvoid = platDistZAvoidConcat;
    dataToDecodeKinPlatformData.platDistZReact = platDistZReactConcat;
    dataToDecodeKinPlatformData.platDistZOpto = platDistZOptoConcat;
    dataToDecodeKinPlatformData.platDistZIti = platDistZItiConcat;
    dataToDecodeKinPlatformData.kinKptsTrialsAvoid = kinKptsTrialsAvoidConcat;
    dataToDecodeKinPlatformData.kinKptsTrialsReact = kinKptsTrialsReactConcat;
    dataToDecodeKinPlatformData.kinKptsTrialsOpto = kinKptsTrialsOptoConcat;
    dataToDecodeKinPlatformData.kinKptsTrialsIti = kinKptsTrialsItiConcat;   

end

function [dataToDecodeKinPlatformData] = getPlatformDecoderSubspaceConcatKin(platDistZAvoid_allMice, platDistZReact_allMice, platDistZOpto_allMice, platDistZIti_allMice, kinKptsTrialsAvoid_allMice, kinKptsTrialsReact_allMice, kinKptsTrialsOpto_allMice, kinKptsTrialsIti_allMice, timescaleSeg, kinMeans, pcaWindow)
% function unpack per-mouse decoding cell array and combine across mice
    numDecoderSampRepeats = 10; % num sampling repeats, should match definition above
    usePcContrast = 0;
    alignPCs = 1;
    smoothingSamples = 1; 
    pcaSamps = [find(pcaWindow(1)>=timescaleSeg,1,'last') find(pcaWindow(2)>=timescaleSeg,1,'last')];
    

    trialMeans = kinMeans.allKptsMeanZ_allMice;
    avoidMeans = squeeze(mean(kinMeans.avoidMeanZ_allMice,1)); %take mean across repeats before PCA
    reactMeans = squeeze(mean(kinMeans.reactMeanZ_allMice,1)); %take mean across repeats before PCA
    numMice = size(platDistZAvoid_allMice,2)/numDecoderSampRepeats;
    numKpts = 36;

    % get aligned, re-orthogonalized PCs for each mouse
    [pcSeg, varExplained] = getKinPCs(trialMeans, avoidMeans, reactMeans, ones(numMice,1).*numKpts, numMice, pcaSamps, usePcContrast, alignPCs);

    for decRpt = 1:numDecoderSampRepeats
        platDistZAvoidConcat = []; 
        platDistZReactConcat = [];
        platDistZOptoConcat = []; 
        platDistZItiConcat = [];
        kinKptsTrialsAvoidConcat = [];
        kinKptsTrialsReactConcat = [];
        kinKptsTrialsOptoConcat = [];
        kinKptsTrialsItiConcat = [];
        for i = 1:numMice %main loop over mice 
            thisMouseIdx = 10*(i-1)+decRpt;
            numTrials = size(platDistZAvoid_allMice{thisMouseIdx},2);
            thisMouseAvoidPlatData = platDistZAvoid_allMice{thisMouseIdx};
            thisMouseAvoidKptData = kinKptsTrialsAvoid_allMice{thisMouseIdx};
            thisMouseReactPlatData = platDistZReact_allMice{thisMouseIdx};
            thisMouseReactKptData = kinKptsTrialsReact_allMice{thisMouseIdx};
            thisMouseOptoPlatData = platDistZOpto_allMice{thisMouseIdx};
            thisMouseOptoKptData = kinKptsTrialsOpto_allMice{thisMouseIdx};
            thisMouseItiPlatData = platDistZIti_allMice{thisMouseIdx};
            thisMouseItiKptData = kinKptsTrialsIti_allMice{thisMouseIdx};
    
            for j = 1:numTrials
                trialProjAvoid = pcSeg{i}'*thisMouseAvoidKptData{j}; %apply PCs to individual trials
                trialProjAvoid = smoothdata(trialProjAvoid, 2, 'gaussian', smoothingSamples); %ex. smooth window 25 is 25*20ms PSTH bins = 500ms
                trialProjReact = pcSeg{i}'*thisMouseReactKptData{j}; %apply PCs to individual trials
                trialProjReact = smoothdata(trialProjReact, 2, 'gaussian', smoothingSamples); %ex. smooth window 25 is 25*20ms PSTH bins = 500ms
                trialProjOpto = pcSeg{i}'*thisMouseOptoKptData{j}; %apply PCs to individual trials
                trialProjOpto = smoothdata(trialProjOpto, 2, 'gaussian', smoothingSamples); %ex. smooth window 25 is 25*20ms PSTH bins = 500ms
                trialProjIti = pcSeg{i}'*thisMouseItiKptData{j}; %apply PCs to individual trials
                trialProjIti = smoothdata(trialProjIti, 2, 'gaussian', smoothingSamples); %ex. smooth window 25 is 25*20ms PSTH bins = 500ms
                platAvoidTrial = thisMouseAvoidPlatData{j};
                platReactTrial = thisMouseReactPlatData{j};
                platOptoTrial = thisMouseOptoPlatData{j};
                platItiTrial = thisMouseItiPlatData{j};
    
                % concatenate all trials over time
                platDistZAvoidConcat = [platDistZAvoidConcat; platAvoidTrial']; %make column vector of concatenated platform distances
                platDistZReactConcat = [platDistZReactConcat; platReactTrial'];
                platDistZOptoConcat = [platDistZOptoConcat; platOptoTrial']; 
                platDistZItiConcat = [platDistZItiConcat; platItiTrial'];
                kinKptsTrialsAvoidConcat = [kinKptsTrialsAvoidConcat; trialProjAvoid'];
                kinKptsTrialsReactConcat = [kinKptsTrialsReactConcat; trialProjReact'];
                kinKptsTrialsOptoConcat = [kinKptsTrialsOptoConcat; trialProjOpto'];
                kinKptsTrialsItiConcat = [kinKptsTrialsItiConcat; trialProjIti'];
            end
        end %end main current mouse loop
         
        % now re-pack all into single matrices:
        dataToDecodeKinPlatformData.platDistZAvoid(decRpt,:) = platDistZAvoidConcat;
        dataToDecodeKinPlatformData.platDistZReact(decRpt,:) = platDistZReactConcat;
        dataToDecodeKinPlatformData.platDistZOpto(decRpt,:) = platDistZOptoConcat;
        dataToDecodeKinPlatformData.platDistZIti(decRpt,:) = platDistZItiConcat;
        dataToDecodeKinPlatformData.kinKptsTrialsAvoid(decRpt,:,:) = kinKptsTrialsAvoidConcat;
        dataToDecodeKinPlatformData.kinKptsTrialsReact(decRpt,:,:) = kinKptsTrialsReactConcat;
        dataToDecodeKinPlatformData.kinKptsTrialsOpto(decRpt,:,:) = kinKptsTrialsOptoConcat;
        dataToDecodeKinPlatformData.kinKptsTrialsIti(decRpt,:,:) = kinKptsTrialsItiConcat;   
    end % end sampling repeats
end

%% kinematics functions
function [kinKptsTrials, timescaleSegKin] = getKinTrials(kinDataKpts, tsData, eventTimes, window)
% function to collect kinematics matrix of windowed trials
% input is struct of raw X/Y/Z positions for each keypoint, so convert back to matrix here but don't normalize
% output is matrix of [#keypoints*trials x timepointsInWindow] around event of interest
    if isempty(eventTimes)
        kinKptsTrials = [];
        timescaleSegKin = [];
        return
    end

    eventTimes = eventTimes(eventTimes<1200); %for datasets with optoTest at end but no video, restrict eventTimes to video recording
    dataTimescale = tsData.timescale(tsData.frameTimeIndeces(1:480000));
    diffDataTimescale = diff(dataTimescale);
    timeBin = mean(diffDataTimescale);
    kptSamplesBefore = round(window(1)/timeBin);
    kptSamplesAfter = round(window(2)/timeBin);
    timescaleSegKin = window(1):timeBin:window(2);
    
    kptNames = {'tail','nose','mouth','r_paw','r_toe','r_ankle','r_knee','l_paw','l_toe','l_ankle','l_knee','chest'};
    for k = kptNames
        kinTrials = zeros(length(eventTimes)*3,abs(kptSamplesBefore)+kptSamplesAfter);
        
        j = 1; 
        for eventIdx = 1:length(eventTimes)
            kptSample = find(dataTimescale>eventTimes(eventIdx), 1, 'first');
            kinTrials(j:j+2,:) = kinDataKpts.(k{1})(:,kptSample+kptSamplesBefore:kptSample+kptSamplesAfter-1);
            j = j + 3; %iterate for X/Y/Z positions
        end   
        kinKptsTrials.(k{1}) = kinTrials;
    end
end

function [behTrialsZ, kinKptsTrials, timescaleSeg] = getBehaviorKinDecEventsBinned(kinDataKptMat, behaviorData, behaviorTimescale, eventTimes, eventWindow, minTrials)
% function to collect z-scored behavior measure in matrix of windows trials, but binned to match decoder predictor timescale
% also collect cell array of [kpts/predictors x positions] for each trial (to match ephys case where # predictors can change for each mouse)
    if isempty(eventTimes)
        behTrialsZ = [];
        kinKptsTrials = [];
        timescaleSeg = [];
        return
    end
    
    % randomly choose minTrials to include
    numEvents = length(eventTimes);
    s = RandStream('dsfmt19937','Seed','shuffle');
    randIdx = randperm(s, numEvents, minTrials); %draw randomly # minimum of avoid/react/opto trials for this mouse
    eventTimesShuf = eventTimes(randIdx);
    numKpts = 36;

    % for kinematics, downsample in time to match ephys PSTH bin size:
    % frame rate = 400Hz, PSTH bin rate = 1/.02 sec = 50Hz, so take mean keypoint over 400/50 = 8 frames
    framesPerPsthBin = 8; %for 20ms ephys bin size (kinematics data s/b matched to this as well) & 10kHz sampRate
    sampTime = (behaviorTimescale(2) - behaviorTimescale(1))*framesPerPsthBin;
    timescaleSeg = eventWindow(1):sampTime:eventWindow(2);
    samplesBefore = round(eventWindow(1)/sampTime);
    samplesAfter = round(eventWindow(2)/sampTime);

    % first normalize data to z-scored version
    behZ = normalize(behaviorData);

    % now bin to same timescale as predictors (i.e. ephys bins) using mean 
    newLength = length(behZ) / framesPerPsthBin;
    ind = [];
    for indk = 1:newLength
        indAdd = ones(1,framesPerPsthBin)*indk;
        ind = [ind indAdd];
    end
    newBehZ = accumarray(ind',behZ(1:length(ind)),[int32(newLength) 1],@mean)'; %use mean of each bin for new value
    timescaleNewBeh = behaviorTimescale(1:framesPerPsthBin:end);
    newKptData = NaN(numKpts, int32(newLength)); %use NaN since it will throw decoder error if not written over correctly
    for kptIdx = 1:numKpts
        newKptData(kptIdx,:) = accumarray(ind(1:480000)',kinDataKptMat(kptIdx, 1:480000),[int32(newLength) 1],@mean)';
    end

    % now collect across events - both behavior & kinematics at same downsampled timescale now, so can use same samples
    kinKptsTrials = [];
    behTrialsZ = [];
    thisBehTrial = zeros(length(timescaleSeg));
    thisKinTrial = zeros(numKpts,length(timescaleSeg));
    for j = 1:length(eventTimesShuf)
        currentSample = find(timescaleNewBeh>eventTimesShuf(j), 1, 'first');
        thisBehTrial = newBehZ(currentSample+samplesBefore:currentSample+samplesAfter-1); 
        behTrialsZ = [behTrialsZ {thisBehTrial}];
        thisKinTrial = newKptData(:,currentSample+samplesBefore:currentSample+samplesAfter-1);
        kinKptsTrials = [kinKptsTrials {thisKinTrial}];
    end

end

function [kinAnglesTrials, timescaleSegKin] = getKinAnglesTrials(kinData, tsData, eventTimes, window, kinParamName, jointNames)
% function to collect kinematics struct of windowed trials for more abstract quantities like joint angle & velocity
% output is matrix of [#jointAngles/etc*trials x timepointsInWindow] around event of interest
    if isempty(eventTimes)
        kinAnglesTrials = [];
        timescaleSegKin = [];
        return
    end

    eventTimes = eventTimes(eventTimes<1200); %for datasets with optoTest at end but no video, restrict eventTimes to video recording
    dataTimescale = tsData.timescale(tsData.frameTimeIndeces(1:480000));
    diffDataTimescale = diff(dataTimescale);
    timeBin = mean(diffDataTimescale);
    kptSamplesBefore = round(window(1)/timeBin);
    kptSamplesAfter = round(window(2)/timeBin);
    timescaleSegKin = window(1):timeBin:window(2);
    
    %names for accessing kinematics struct hierarchy:
    kinNames1 = kinParamName; % jointAnglesNorm, jointAngles, jointVel
    kinNames2 = {'HL'}; 
    kinNames3 = {'L', 'R'}; 
    kinNames4 = jointNames;

    for k1 = kinNames1
        for k2 = kinNames2
            for k3 = kinNames3
                for k4 = kinNames4
                    kinTrials = zeros(length(eventTimes),abs(kptSamplesBefore)+kptSamplesAfter);
                    
                    for eventIdx = 1:length(eventTimes)
                        kptSample = find(dataTimescale>eventTimes(eventIdx), 1, 'first');
                        kinTrials(eventIdx,:) = kinData.(k1{1}).(k2{1}).(k3{1}).(k4{1})(kptSample+kptSamplesBefore:kptSample+kptSamplesAfter-1);
                    end   
                    kinAnglesTrials.(k1{1}).(k2{1}).(k3{1}).(k4{1}) = kinTrials;
                end
            end
        end
    end
end

function [kinKptsTrialsOut] = concatKinTrials(kinKptsTrialsIn, newkinKptsTrials)
% function to concatenate kinematics struct of windowed trials
    kptNames = {'tail','nose','mouth','r_paw','r_toe','r_ankle','r_knee','l_paw','l_toe','l_ankle','l_knee','chest'};
    if isempty(kinKptsTrialsIn)
        kinKptsTrialsOut = newkinKptsTrials; %initialize struct
    elseif isempty(newkinKptsTrials)
        kinKptsTrialsOut = kinKptsTrialsIn;
    else 
        for k = kptNames 
            kinKptsTrialsOut.(k{1}) = [kinKptsTrialsIn.(k{1}); newkinKptsTrials.(k{1})];
        end
    end
end

function [kinAnglesTrialsOut] = concatKinAngleTrials(kinAngleVsTrialsIn, newkinAngleVsTrials, kinParamName, jointNames)
% function to concatenate kinematics angles/etc struct of windowed trials
    
    %names for accessing kinematics struct hierarchy; must match 'getKinAnglesTrials' function definition
    kinNames1 = kinParamName; % jointAnglesNorm, jointAngles, jointVel
    kinNames2 = {'HL'}; 
    kinNames3 = {'L', 'R'}; 
    kinNames4 = jointNames;

    if isempty(kinAngleVsTrialsIn)
        kinAnglesTrialsOut = newkinAngleVsTrials; %initialize struct
    elseif isempty(newkinAngleVsTrials)
        kinAnglesTrialsOut = kinAngleVsTrialsIn;
    else 
        for k1 = kinNames1
            for k2 = kinNames2
                for k3 = kinNames3
                    for k4 = kinNames4
                        kinAnglesTrialsOut.(k1{1}).(k2{1}).(k3{1}).(k4{1}) = [kinAngleVsTrialsIn.(k1{1}).(k2{1}).(k3{1}).(k4{1}); newkinAngleVsTrials.(k1{1}).(k2{1}).(k3{1}).(k4{1})];
                    end
                end
            end
        end
    end

end

function [kinTrialsZ, timescaleSeg] = getKinZscoreEvents(kinData, frameTimes, eventTimes, eventWindow)
% function to collect z-scored kinematics matrix of windowed trials
% output is matrix of [keypoints x trials x timepointsInWindow]
    
    keypoints = kinData.kptMat; %already normalized
    numKeypoints = size(keypoints,1);
    sampTime = frameTimes(2) - frameTimes(1);
    timescaleSeg = eventWindow(1):sampTime:eventWindow(2);
    samplesBefore = round(eventWindow(1)/sampTime);
    samplesAfter = round(eventWindow(2)/sampTime);
    kinTrialsZ = zeros(numKeypoints, length(eventTimes),length(timescaleSeg)); %neurons x trials x time
    
    % now collect across events
    for i = 1:numKeypoints
        for j = 1:length(eventTimes)
            currentSample = find(frameTimes>eventTimes(j), 1, 'first');
            kinTrialsZ(i,j,:) = keypoints(i,currentSample+samplesBefore:currentSample+samplesAfter-1);
        end
    end

end

function plotKinAvoidVsReactPcTrajectoriesAllMice(allUnitsMean, allUnitsAvoidMean, allUnitsReactMean, numKptsIn, timescaleSeg, pcaWindow, plotWindow)
% very similar to function for ephys but with minor modifications for kinematics
% input is all N x T matrices that concatenate trial averages (over avoid / react trials) across mice/sessions
% also input numGoodUnits matrix that keeps track of # neuron splits between mice/sessions, and time windows for PCA calculation and plotting
    numMice = size(numKptsIn,2);
    pcaSamps = [find(pcaWindow(1)>=timescaleSeg,1,'last') find(pcaWindow(2)>=timescaleSeg,1,'last')];
    plotSamps = [find(plotWindow(1)>=timescaleSeg,1,'last') find(plotWindow(2)>=timescaleSeg,1,'last')];
    timescalePlot = timescaleSeg(plotSamps(1):plotSamps(2));
    plot3d = 1;
    pcaIdx1 = 1; %make this variable to easy to change which PCs to view if necessary
    pcaIdx2 = 2;
    pcaIdx3 = 3;

    % for re-orthogonalization, must limit PCs to maximum of numKpts
    numPcsToUse = 10;
    smoothingSamples = 20; %ex. smooth window samples = 10 is 10*2.5ms frames = 25ms; for causal half-kernel use [X 0]
    usePcContrast = 0; % whether we perform contrastive PCA on data1 vs data2, or vanilla PCA on allData 
    alignPCs = 0; % whether to align all mice/sessions into common subspace for computing PCs; otherwise PCs are computed independently for each mouse/session (default for simplicity)

    %plotting constants
    fontSz = 16;
    numColorPts = 256; %standard colormap pts
    cmapIdx = numColorPts - round(numColorPts/4); %for single colors for individual PC plots
    R = linspace(0.3,0.9,numColorPts)'; %go from half-black to full color
    G = linspace(0.1,0.5,numColorPts)';
    B = zeros(numColorPts,1);
    cmap1 = [R G B];
    R = zeros(numColorPts,1);
    G = linspace(0.0,0.4,numColorPts)';
    B = linspace(0.3,0.9,numColorPts)';
    cmap2 = [R G B];
    psub1 = [];
    psub2 = [];
    cmapMice = linspecer(numMice);
 
    % optional sanity plot, look at data before PCA:
%     climAll = [-0.5 0.5];
%     xData = avoidDataCent(1:3:end-2,:) - reactDataCent(1:3:end-2,:);
%     yData = avoidDataCent(1:3:end-1,:) - reactDataCent(1:3:end-1,:);
%     zData = avoidDataCent(3:3:end,:) - reactDataCent(3:3:end,:);
%     for k = 1:12
%         figure;
%         hold on;
%         kptIdx = k;
%         for i = 1:numMice
%             plot(zData(kptIdx,:))
%             kptIdx = kptIdx + 12;
%         end
%     end
% 
%     figure;
%     tiledlayout(3,1,'TileSpacing', 'compact')
%     nexttile
%     hold on;
%     imagesc(timescaleSeg, 1:1:size(xData,1), xData)
%     set(gca, 'YDir', 'reverse');
%     axis tight;
%     zeroXDottedLine;
%     drawMouseSeparatorsHeatmap(ones(1,numMice).*12)
%     colormap(greenMag)
%     clim(climAll);
%     nexttile
%     hold on;
%     imagesc(timescaleSeg, 1:1:size(yData,1), yData)
%     set(gca, 'YDir', 'reverse');
%     axis tight;
%     zeroXDottedLine;
%     drawMouseSeparatorsHeatmap(ones(1,numMice).*12)
%     colormap(greenMag)
%     clim(climAll);
%     nexttile
%     hold on;
%     imagesc(timescaleSeg, 1:1:size(zData,1), zData)
%     set(gca, 'YDir', 'reverse');
%     axis tight;
%     zeroXDottedLine;
%     drawMouseSeparatorsHeatmap(ones(1,numMice).*12)
%     colormap(greenMag)
%     clim(climAll);

    [pcSeg, varExplained] = getKinPCs(allUnitsMean, allUnitsAvoidMean, allUnitsReactMean, numKptsIn, numMice, pcaSamps, usePcContrast, alignPCs); 
    
    % optional plot to look at variance explained & distribution of PC weights
    hVarFig = figure; 
    plotCap = 5; %cap the num PCs to plot to see detail
    varMin = 90;
    hold on;
    if alignPCs
        plot(cumsum(varExplained), 'k-', 'LineWidth', 3);
    else
        for i = 1:numMice
            plot(cumsum(varExplained{i}), 'k--', 'LineWidth', 0.5);
            varExplainedAllMice(i,:) = varExplained{i}(1:plotCap); %collect for mean
        end
        varExplainedMeanAcrossMice = mean(varExplainedAllMice,1);
%         cumulativeMeanVarExplainedAllMice = mean(cumsum(varExplainedAllMice(:,1:3),2),1) %optionally look at mean variance explained
%         cumulativeStdDevVarExplainedAllMice = std(cumsum(varExplainedAllMice(:,1:3),2),1)
        plot(cumsum(varExplainedMeanAcrossMice), 'k-', 'LineWidth', 3);
    end
%     plot([numPcsToUse numPcsToUse], [varMin 100], 'b--');
    hold off;
    axis([0.9 plotCap varMin 100])
    xlabel('PCs')
    ylabel('% variance explained')
    set(gcf,'color','w'); %set figure background white
    makepretty



    % plot PC weights to help interpret PCs:------------------------------------
    numKpts = numKptsIn(1);
    hlIdxOffsets = [1:21]; % for each mouse, these are hindlimb points + base of tail which covaries with hindlimbs 7x3 = 21/36
    flIdxOffsets = [22:27]; % these are forelimb paw points 6/36
    otherIdxOffsets = [28:36]; % these are nose / mouth / chest points 9/36
    hlSeparator = [length(hlIdxOffsets) length(hlIdxOffsets)] + 0.5;
    flSeparator = [length(hlIdxOffsets)+length(flIdxOffsets) length(hlIdxOffsets)+length(flIdxOffsets)] + 0.5;

    %gather up all PC coefficients by keypoint type across mice:
    hlPCs = zeros(numMice,3,length(hlIdxOffsets));
    flPCs = zeros(numMice,3,length(flIdxOffsets));
    otherPCs = zeros(numMice,3,length(otherIdxOffsets));
    for i = 1:numMice
        hlPCs(i,1,:) = pcSeg{i}(hlIdxOffsets, pcaIdx1);
        hlPCs(i,2,:) = pcSeg{i}(hlIdxOffsets, pcaIdx2);
        hlPCs(i,3,:) = pcSeg{i}(hlIdxOffsets, pcaIdx3);
        flPCs(i,1,:) = pcSeg{i}(flIdxOffsets, pcaIdx1);
        flPCs(i,2,:) = pcSeg{i}(flIdxOffsets, pcaIdx2);
        flPCs(i,3,:) = pcSeg{i}(flIdxOffsets, pcaIdx3);
        otherPCs(i,1,:) = pcSeg{i}(otherIdxOffsets, pcaIdx1);
        otherPCs(i,2,:) = pcSeg{i}(otherIdxOffsets, pcaIdx2);
        otherPCs(i,3,:) = pcSeg{i}(otherIdxOffsets, pcaIdx3);
    end

    % plot spread of PC coefficients across first 3 PCs 
    figure; 
    tiledlayout(3,1,'TileSpacing', 'tight'); nexttile(1); hold on; nexttile(2); hold on; nexttile(3); hold on;
    for i = 1:numMice
%         if (i == 1) | (i == 5) | (i == 5) %plot particular mouse in another color for debug
%             nexttile(1)
%             plot([squeeze(hlPCs(i,1,:))' squeeze(flPCs(i,1,:))' squeeze(otherPCs(i,1,:))'], 'Color', [1 0 0]); % cmapMice(i,:)) %lines options to see diff mice
%             nexttile(2)
%             plot([squeeze(hlPCs(i,2,:))' squeeze(flPCs(i,2,:))' squeeze(otherPCs(i,2,:))'], 'Color', [1 0 0]); % cmapMice(i,:))
%             nexttile(3)
%             plot([squeeze(hlPCs(i,3,:))' squeeze(flPCs(i,3,:))' squeeze(otherPCs(i,3,:))'], 'Color', [1 0 0]); % cmapMice(i,:))     
%         else
            nexttile(1)
            plot([squeeze(hlPCs(i,1,:))' squeeze(flPCs(i,1,:))' squeeze(otherPCs(i,1,:))'], 'Color', cmapMice(i,:)); % cmapMice(i,:)) %lines options to see diff mice
            nexttile(2)
            plot([squeeze(hlPCs(i,2,:))' squeeze(flPCs(i,2,:))' squeeze(otherPCs(i,2,:))'], 'Color', cmapMice(i,:)); % cmapMice(i,:))
            nexttile(3)
            plot([squeeze(hlPCs(i,3,:))' squeeze(flPCs(i,3,:))' squeeze(otherPCs(i,3,:))'], 'Color', cmapMice(i,:)); % cmapMice(i,:))
%         end
    end
    nexttile(1)
    yLims = get(gca, 'YLim');
    y = [yLims(1) yLims(2)];
    plot(hlSeparator, y, 'r--', 'LineWidth', 2); %mark end of HL points
    plot(flSeparator, y, 'r--', 'LineWidth', 2); %mark end of FL points
    for j = 1:3:numKpts-1
        plot([j j], y, 'k--', 'LineWidth', 1); %mark all left-right (X in kptMat) coordinates, X front-back, Z side-side
        plot([j+1 j+1], y, 'b--', 'LineWidth', 1); %mark all front-back (Y in kptMat) coordinates, X front-back, Z side-side
        plot([j+2 j+2], y, 'r--', 'LineWidth', 1); %mark all vertical (Z in kptMat) coordinates, X front-back, Z side-side
    end
    xticks([1:36])
    xticklabels({})
    axis tight;
    zeroYDottedLine
    hold off;
    ylabel(['PC' num2str(pcaIdx1) ' coeff.'])
%     makepretty

    nexttile(2)
    yLims = get(gca, 'YLim');
    y = [yLims(1) yLims(2)];
    plot(hlSeparator, y, 'r--', 'LineWidth', 2); 
    plot(flSeparator, y, 'r--', 'LineWidth', 2); 
    for j = 1:3:numKpts-1
        plot([j j], y, 'k--', 'LineWidth', 1); %mark all left-right (X in kptMat) coordinates, X front-back, Z side-side
        plot([j+1 j+1], y, 'b--', 'LineWidth', 1); %mark all front-back (Y in kptMat) coordinates, X front-back, Z side-side
        plot([j+2 j+2], y, 'g--', 'LineWidth', 1); %mark all vertical (Z in kptMat) coordinates, X front-back, Z side-side
    end
    xticks([1:36])
    xticklabels({})
    axis tight;
    zeroYDottedLine
    hold off;
    ylabel(['PC' num2str(pcaIdx2) ' coeff.'])
%     makepretty

    nexttile(3)
    yLims = get(gca, 'YLim');
    y = [yLims(1) yLims(2)];
    plot(hlSeparator, y, 'r--', 'LineWidth', 2); 
    plot(flSeparator, y, 'r--', 'LineWidth', 2);
    for j = 1:3:numKpts-1
        plot([j j], y, 'k--', 'LineWidth', 1); %mark all left-right (X in kptMat) coordinates, X front-back, Z side-side
        plot([j+1 j+1], y, 'b--', 'LineWidth', 1); %mark all front-back (Y in kptMat) coordinates, X front-back, Z side-side
        plot([j+2 j+2], y, 'g--', 'LineWidth', 1); %mark all vertical (Z in kptMat) coordinates, X front-back, Z side-side
    end
    xticks([1:36])
    xticklabels({[],[],'base tail vert',[],[],'R toe vert',[],[],'R ankle vert',[],[],'R knee vert',[],[],'L toe vert',[],[],'L ankle vert',[],[],'L knee vert',[],[],'R finger vert',[],[],'L finger vert',[],[],'nose vert',[],[],'mouth vert',[],[],'chest vert'});
    xtickangle(60)
    axis tight;
    zeroYDottedLine
    hold off;
    ylabel(['PC' num2str(pcaIdx3) ' coeff.'])
    set(gcf,'color','w'); %set figure background white  
%     makepretty
    % --------------------------------------------------------------------------
    
    currentEndUnit = 0;
    if plot3d
        pc3dFigure = figure; 
        tiledlayout('flow');
    end
    pc1OnlyFigure = figure; 
    tiledlayout(1,3,'TileSpacing', 'compact');

    for i = 1:numMice
        currentStartUnit = currentEndUnit + 1;
        currentEndUnit = currentEndUnit + numKptsIn(i);

        % project mean avoid/ react activity for each mouse onto its respective PC matrix to get trajectories for each mouse & overlay
        avoidMeanProj = pcSeg{i}'*allUnitsAvoidMean(currentStartUnit:currentEndUnit,plotSamps(1):plotSamps(2)); %apply PCs to the mean responses across trials
        reactMeanProj = pcSeg{i}'*allUnitsReactMean(currentStartUnit:currentEndUnit,plotSamps(1):plotSamps(2));
        avoidMeanProj = smoothdata(avoidMeanProj, 2, 'gaussian', smoothingSamples); %ex. smooth window 10 is 10*50ms PSTH bins = 500ms
        reactMeanProj = smoothdata(reactMeanProj, 2, 'gaussian', smoothingSamples);

        colorsIdx = floor(linspace(1, numColorPts, size(avoidMeanProj,2))); %evenly spaced colors that show time dimension
    
        % loop for plotting PC trajectories in 3D
        if plot3d
            figure(pc3dFigure)
            nexttile; 
            hold on;
            zeroSample = find(timescalePlot>0,1,'first'); %for plotting midpoint
            for dataPoint = 1:size(avoidMeanProj,2)-1
    
%                 if dataPoint == size(avoidMeanProj,2)-1
%                     psub1 = plot3(avoidMeanProj(pcaIdx1,dataPoint:dataPoint+1),avoidMeanProj(pcaIdx2,dataPoint:dataPoint+1),avoidMeanProj(pcaIdx3,dataPoint:dataPoint+1),'Color',cmap1(colorsIdx(dataPoint),:), 'LineWidth', 3, 'DisplayName','avoid'); %subset for legend display name
%                     psub2 = plot3(reactMeanProj(pcaIdx1,dataPoint:dataPoint+1),reactMeanProj(pcaIdx2,dataPoint:dataPoint+1),reactMeanProj(pcaIdx3,dataPoint:dataPoint+1),'Color',cmap2(colorsIdx(dataPoint),:),  'LineWidth', 3, 'DisplayName','react'); %subset for legend display name
%                 else
%                     plot3(avoidMeanProj(pcaIdx1,dataPoint:dataPoint+1),avoidMeanProj(pcaIdx2,dataPoint:dataPoint+1),avoidMeanProj(pcaIdx3,dataPoint:dataPoint+1),'Color',cmap1(colorsIdx(dataPoint),:), 'LineWidth', 3); 
%                     plot3(reactMeanProj(pcaIdx1,dataPoint:dataPoint+1),reactMeanProj(pcaIdx2,dataPoint:dataPoint+1),reactMeanProj(pcaIdx3,dataPoint:dataPoint+1),'Color',cmap2(colorsIdx(dataPoint),:),  'LineWidth', 3); 
%                     if dataPoint == zeroSample %plot marker at t=0 extension threshold
%                         plot3(avoidMeanProj(1,dataPoint),avoidMeanProj(2,dataPoint),avoidMeanProj(3,dataPoint), 'Marker', 'o', 'MarkerFaceColor','m', 'MarkerSize', 9, 'MarkerEdgeColor','none'); 
%                         plot3(reactMeanProj(1,dataPoint),reactMeanProj(2,dataPoint),reactMeanProj(3,dataPoint),'Marker', 'o','MarkerFaceColor','m', 'MarkerSize', 9, 'MarkerEdgeColor','none'); 
%                     end
%                     if dataPoint == 1 %plot marker beginning
%                         plot3(avoidMeanProj(1,dataPoint),avoidMeanProj(2,dataPoint),avoidMeanProj(3,dataPoint), 'Marker', 'o', 'MarkerFaceColor','g', 'MarkerSize', 9, 'MarkerEdgeColor','none'); 
%                         plot3(reactMeanProj(1,dataPoint),reactMeanProj(2,dataPoint),reactMeanProj(3,dataPoint),'Marker', 'o','MarkerFaceColor','g', 'MarkerSize', 9, 'MarkerEdgeColor','none'); 
%                     end
%                 end
                
                % simpler 2d plot option
                if dataPoint == size(avoidMeanProj,2)-1
                    psub1 = plot(avoidMeanProj(pcaIdx1,dataPoint:dataPoint+1),avoidMeanProj(pcaIdx2,dataPoint:dataPoint+1),'Color',cmap1(colorsIdx(dataPoint),:), 'LineWidth', 3, 'DisplayName','avoid'); %subset for legend display name
                    psub2 = plot(reactMeanProj(pcaIdx1,dataPoint:dataPoint+1),reactMeanProj(pcaIdx2,dataPoint:dataPoint+1),'Color',cmap2(colorsIdx(dataPoint),:),  'LineWidth', 3, 'DisplayName','react'); %subset for legend display name
                else
                    plot(avoidMeanProj(pcaIdx1,dataPoint:dataPoint+1),avoidMeanProj(pcaIdx2,dataPoint:dataPoint+1),'Color',cmap1(colorsIdx(dataPoint),:), 'LineWidth', 3); 
                    plot(reactMeanProj(pcaIdx1,dataPoint:dataPoint+1),reactMeanProj(pcaIdx2,dataPoint:dataPoint+1),'Color',cmap2(colorsIdx(dataPoint),:),  'LineWidth', 3); 
                    if dataPoint == zeroSample %plot marker at t=0 extension threshold
                        plot(avoidMeanProj(1,dataPoint),avoidMeanProj(2,dataPoint), 'Marker', 'o', 'MarkerFaceColor','m', 'MarkerSize', 9, 'MarkerEdgeColor','none'); 
                        plot(reactMeanProj(1,dataPoint),reactMeanProj(2,dataPoint),'Marker', 'o','MarkerFaceColor','m', 'MarkerSize', 9, 'MarkerEdgeColor','none'); 
                    end
                    if dataPoint == 1 %plot marker beginning
                        plot(avoidMeanProj(1,dataPoint),avoidMeanProj(2,dataPoint), 'Marker', 'o', 'MarkerFaceColor','g', 'MarkerSize', 9, 'MarkerEdgeColor','none'); 
                        plot(reactMeanProj(1,dataPoint),reactMeanProj(2,dataPoint),'Marker', 'o','MarkerFaceColor','g', 'MarkerSize', 9, 'MarkerEdgeColor','none'); 
                    end
                end


            end % end dataPoint loop
            hold off;
%             view(45,30);
            axis tight;
        end % end if(plot3d)

        % separate single PC trajectories plots - check for temporal structure in higher PCs
        figure(pc1OnlyFigure); 
        nexttile(1); hold on;
        psub3 = plot(timescalePlot, avoidMeanProj(pcaIdx1,:),'Color',cmap1(cmapIdx,:), 'LineWidth', 0.1, 'DisplayName','avoid'); 
        psub4 = plot(timescalePlot, reactMeanProj(pcaIdx1,:),'Color',cmap2(cmapIdx,:), 'LineWidth', 0.1, 'DisplayName','react'); 
        %cla; % uncomment & set breakpoint here to look at one mouse at a time
        pc1ProjsAvoid(i,:) = avoidMeanProj(pcaIdx1,:); % save for bounded line / mean plot below
        pc1ProjsReact(i,:) = reactMeanProj(pcaIdx1,:);
        nexttile(2); hold on;
        psub5 = plot(timescalePlot, avoidMeanProj(pcaIdx2,:),'Color',cmap1(cmapIdx,:), 'LineWidth', 0.1, 'DisplayName','avoid'); 
        psub6 = plot(timescalePlot, reactMeanProj(pcaIdx2,:),'Color',cmap2(cmapIdx,:), 'LineWidth', 0.1, 'DisplayName','react'); 
        pc2ProjsAvoid(i,:) = avoidMeanProj(pcaIdx2,:); % save for bounded line / mean plot below
        pc2ProjsReact(i,:) = reactMeanProj(pcaIdx2,:);
        nexttile(3); hold on;
        psub7 = plot(timescalePlot, avoidMeanProj(pcaIdx3,:),'Color',cmap1(cmapIdx,:), 'LineWidth', 0.1, 'DisplayName','avoid'); 
        psub8 = plot(timescalePlot, reactMeanProj(pcaIdx3,:),'Color',cmap2(cmapIdx,:), 'LineWidth', 0.1, 'DisplayName','react'); 
        pc3ProjsAvoid(i,:) = avoidMeanProj(pcaIdx3,:); % save for bounded line / mean plot below
        pc3ProjsReact(i,:) = reactMeanProj(pcaIdx3,:);

    end %end loop over mice

    % finish 3d PCs figure:
    if plot3d
        figure(pc3dFigure)
        xlabel(['PC' num2str(pcaIdx1) ' (a.u.)'], 'FontSize', fontSz)
        ylabel(['PC' num2str(pcaIdx2) ' (a.u.)'], 'FontSize', fontSz)
        zlabel(['PC' num2str(pcaIdx3) ' (a.u.)'], 'FontSize', fontSz)
        if ~isempty(psub1)
            legend([psub1 psub2], 'Box', 'off');
        end
        set(gcf,'color','w'); %set figure background white
    end

    % finish single PC figures:
    figure(pc1OnlyFigure);
    nexttile(1); 
    [meanAvoidPC1, semAvoidPC1] = grpstats(pc1ProjsAvoid,[],{'mean' 'sem'});
    [meanReactPC1, semReactPC1] = grpstats(pc1ProjsReact,[],{'mean' 'sem'});
    plot(timescalePlot, meanAvoidPC1, 'Color',cmap1(cmapIdx,:), 'LineWidth', 3);
    plot(timescalePlot, meanReactPC1, 'Color',cmap2(cmapIdx,:), 'LineWidth', 3);
%     boundedline(timescalePlot, meanAvoidPC1, semAvoidPC1, 'cmap', cmap1(cmapIdx,:));
%     boundedline(timescalePlot, meanReactPC1, semReactPC1, 'cmap', cmap2(cmapIdx,:));
    zeroXDottedLine;
    hold off;
    xlabel('time (sec)', 'FontSize', fontSz)
    ylabel(['PC' num2str(pcaIdx1) ' (a.u.)'], 'FontSize', fontSz)
%     legend([psub3 psub4], 'Box', 'off');
    axis tight;

    nexttile(2); 
    [meanAvoidPC2, semAvoidPC2] = grpstats(pc2ProjsAvoid,[],{'mean' 'sem'});
    [meanReactPC2, semReactPC2] = grpstats(pc2ProjsReact,[],{'mean' 'sem'});
    plot(timescalePlot, meanAvoidPC2, 'Color',cmap1(cmapIdx,:), 'LineWidth', 3);
    plot(timescalePlot, meanReactPC2, 'Color',cmap2(cmapIdx,:), 'LineWidth', 3);
%     boundedline(timescalePlot, meanAvoidPC2, semAvoidPC2, 'cmap', cmap1(cmapIdx,:));
%     boundedline(timescalePlot, meanReactPC2, semReactPC2, 'cmap', cmap2(cmapIdx,:));
    zeroXDottedLine;
    hold off;
    xlabel('time (sec)', 'FontSize', fontSz)
    ylabel(['PC' num2str(pcaIdx2) ' (a.u.)'], 'FontSize', fontSz)
%     legend([psub5 psub6], 'Box', 'off');
    axis tight;

    nexttile(3); 
    [meanAvoidPC3, semAvoidPC3] = grpstats(pc3ProjsAvoid,[],{'mean' 'sem'});
    [meanReactPC3, semReactPC3] = grpstats(pc3ProjsReact,[],{'mean' 'sem'});
    plot(timescalePlot, meanAvoidPC3, 'Color',cmap1(cmapIdx,:), 'LineWidth', 3);
    plot(timescalePlot, meanReactPC3, 'Color',cmap2(cmapIdx,:), 'LineWidth', 3);
%     boundedline(timescalePlot, meanAvoidPC3, semAvoidPC3, 'cmap', cmap1(cmapIdx,:));
%     boundedline(timescalePlot, meanReactPC3, semReactPC3, 'cmap', cmap2(cmapIdx,:));
    zeroXDottedLine;
    hold off;
    xlabel('time (sec)', 'FontSize', fontSz)
    ylabel(['PC' num2str(pcaIdx3) ' (a.u.)'], 'FontSize', fontSz)
    legend([psub7 psub8], 'Box', 'off');
    axis tight;
    set(gcf,'color','w'); %set figure background white

end

function plotKinPcTrajectoriesAllMiceAcrossTrials(allKptsZConcat, timescaleSeg)
% very similar to function above but with minor modifications for kinematics
% input is all N x T matrices that concatenate all trials across mice/sessions
% windows for plotting and PCA are alreasy set by timescaleSeg & getTrials function
    trialSamps = length(timescaleSeg);
    numTrials = size(allKptsZConcat,2)/trialSamps;

    %plotting constants
    fontSz = 16;
    numColorPts = numTrials; %standard colormap pts
    R = zeros(numColorPts,1); %go from half-black to full color
    G = linspace(0,1,numColorPts)';
    B = linspace(0,1,numColorPts)';
    T = ones(numColorPts,1).*0.5; %4 element for transparency
    cmap1 = [R G B T];
    pcaIdx1 = 2; %make this variable to easy to change which PCs to view if necessary
    pcaIdx2 = 4;
    pcaIdx3 = 8;

    % first subtract grand mean over time, even though mean z-score data should be close to centered already
    allData = allKptsZConcat;
    allDataCent = allData - mean(allData,2);

    [PC,V,varExplained] = pcaCov(allDataCent);  %PCs in columns in decreasing VarExplained 
    kptProj = PC'*allKptsZConcat; %project keypoints for trials onto respective PCs matrix to get trajectories across trials
    kptProjSmooth = smoothdata(kptProj, 2, 'gaussian', 20); %ex. smooth window 10 is 10*2.5ms frame time = 25ms

    % now split projection across trials and plot
%     zeroSample = find(timescaleSeg>0,1,'first'); %for plotting midpoint
    figure;
    hold on;
    trialStartIdx = 1;
    for i = 1:numTrials
        trialEndIdx = trialStartIdx + trialSamps - 1;
%         if mod(i,10) == 0 %plot subset of trials
            plot3(kptProjSmooth(pcaIdx1,trialStartIdx:trialEndIdx),kptProjSmooth(pcaIdx2,trialStartIdx:trialEndIdx),kptProjSmooth(pcaIdx3,trialStartIdx:trialEndIdx),'Color',cmap1(i,:), 'LineWidth', 1); 
%         end
        trialStartIdx = trialStartIdx + trialSamps;
    end

    hold off;
    view(45,30);
    xlabel(['PC' num2str(pcaIdx1) ' (a.u.)'], 'FontSize', fontSz)
    ylabel(['PC' num2str(pcaIdx2) ' (a.u.)'], 'FontSize', fontSz)
    zlabel(['PC' num2str(pcaIdx3) ' (a.u.)'], 'FontSize', fontSz)
    set(gcf,'color','w'); %set figure background white

    % also plot trial overlay for first 8 PCs vs time to make sure we're not missing anything with first 3:
    figure; 
    pcsToPlot = [2 4 8];
    tiledlayout(length(pcsToPlot),1);
    for pcIdx = pcsToPlot
        nexttile
        hold on
        trialStartIdx = 1;
        for i = 1:numTrials
            trialEndIdx = trialStartIdx + trialSamps - 1;
            plot(timescaleSeg, kptProjSmooth(pcIdx,trialStartIdx:trialEndIdx),'Color',cmap1(i,:), 'LineWidth', 1); 
            trialStartIdx = trialStartIdx + trialSamps;
        end
        hold off
        xlabel(['time (sec)   var. expl. = ' num2str(varExplained(pcIdx))], 'FontSize', fontSz)
        ylabel(['PC' num2str(pcIdx) ' (a.u.)'], 'FontSize', fontSz)
    end
    set(gcf,'color','w'); %set figure background white

end

function plotKinAvoidVsReactPcDistancesAllMice(allUnitsMean, allUnitsAvoidMean, allUnitsReactMean, allUnitsMeanZFR_shuf1, allUnitsMeanZFR_shuf2, numKpts, timescaleSeg, pcaWindow, plotWindow)
% only plot Euclidean distances, with shuffle controls
% input is all N x T matrices that concatenate trial averages (over avoid / react trials) across mice/sessions
% also input numGoodUnits matrix that keeps track of # neuron splits between mice/sessions, and time windows for PCA calculation and plotting
    numMice = size(numKpts,2);
    pcaSamps = [find(pcaWindow(1)>=timescaleSeg,1,'last') find(pcaWindow(2)>=timescaleSeg,1,'last')];
    plotSamps = [find(plotWindow(1)>=timescaleSeg,1,'last') find(plotWindow(2)>=timescaleSeg,1,'last')];
    timescalePlot = timescaleSeg(plotSamps(1):plotSamps(2));
    
    usePcContrast = 1; % whether we perform contrastive PCA on data1 vs data2, or vanilla PCA on allData 
    alignPCs = 1; % whether to align all mice/sessions into common subspace for computing PCs; otherwise PCs are computed independently for each mouse/session (default for simplicity)
    smoothingSamples = 20; %ex. smooth window samples 20*2.5ms frames = 50ms; for causal half-kernel use [X 0]
    %plotting constants
    fontSz = 16;
    
    [pcSeg, varExplained] = getKinPCs(allUnitsMean, allUnitsAvoidMean, allUnitsReactMean, numKpts, numMice, pcaSamps, usePcContrast, alignPCs); 

    currentEndUnit = 0;
    savedEuclDistance = [];
    euclFig = figure; hold on;

    for i = 1:numMice
        currentStartUnit = currentEndUnit + 1;
        currentEndUnit = currentEndUnit + numKpts(i);

        % project mean avoid/ react activity for each mouse onto its respective PC matrix to get trajectories for each mouse & overlay
        avoidMeanProj = pcSeg{i}'*allUnitsAvoidMean(currentStartUnit:currentEndUnit,plotSamps(1):plotSamps(2)); %apply PCs to the mean responses across trials
        reactMeanProj = pcSeg{i}'*allUnitsReactMean(currentStartUnit:currentEndUnit,plotSamps(1):plotSamps(2));
        avoidMeanProj = smoothdata(avoidMeanProj, 2, 'gaussian', smoothingSamples); %ex. smooth window 10 is 10*2.5ms frames = 25ms
        reactMeanProj = smoothdata(reactMeanProj, 2, 'gaussian', smoothingSamples);
        
        euclDistance = zeros(size(avoidMeanProj,2)-1,1);
        for dataPoint = 1:size(avoidMeanProj,2)-1
            % keep track of Euclidean distances between trajectories for plotting below, using all dimensions available
            euclDistance(dataPoint) = norm(avoidMeanProj(:,dataPoint) - reactMeanProj(:,dataPoint));
        end
        
%         plot(timescalePlot(1:end-1), euclDistance./sqrt(numKpts(i)), 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 0.1)
        savedEuclDistance(i, :) = euclDistance./sqrt(numKpts(i)); %save for bounded line below
    end %end loop over mice

    if ~isempty(savedEuclDistance)
        [meanEuclDistance, semEuclDistance] = grpstats(savedEuclDistance,[],{'mean' 'sem'}); %mean across columns/timepoints
%         nexttile(1);
        [hl1,~] = boundedline(timescalePlot(1:end-1), meanEuclDistance, semEuclDistance, 'cmap', [0.4940 0.1840 0.5560]);

    end

    %%%%%%%%%%%%%%%% NOW PLOT CONTROL / SHUFFLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    savedEuclDistanceShuf = [];
    currentEndUnit = 0;
    for i = 1:numMice
        currentStartUnit = currentEndUnit + 1;
        currentEndUnit = currentEndUnit + numKpts(i);

        % project mean shuffle activity onto the same PCs used to differentiate avoid & react
        shuf1MeanProj = pcSeg{i}'*allUnitsMeanZFR_shuf1(currentStartUnit:currentEndUnit,plotSamps(1):plotSamps(2)); %apply PCs to the mean responses across trials
        shuf2MeanProj = pcSeg{i}'*allUnitsMeanZFR_shuf2(currentStartUnit:currentEndUnit,plotSamps(1):plotSamps(2));
        shuf1MeanProj = smoothdata(shuf1MeanProj, 2, 'gaussian', smoothingSamples); %ex. smooth window 10 is 10*50ms PSTH bins = 500ms
        shuf2MeanProj = smoothdata(shuf2MeanProj, 2, 'gaussian', smoothingSamples);
        
        euclDistance = zeros(size(shuf1MeanProj,2)-1,1);
        for dataPoint = 1:size(shuf1MeanProj,2)-1
            % keep track of Euclidean distances between trajectories for plotting below
            euclDistance(dataPoint) = norm(shuf1MeanProj(:,dataPoint) - shuf2MeanProj(:,dataPoint));
        end

        savedEuclDistanceShuf(i, :) = euclDistance./sqrt(numKpts(i)); %save for bounded line below
    end %end loop over mice

    if ~isempty(savedEuclDistanceShuf)
        [meanEuclDistance, semEuclDistance] = grpstats(savedEuclDistanceShuf,[],{'mean' 'sem'}); %mean across columns/timepoints
        [hl2,~] = boundedline(timescalePlot(1:end-1), meanEuclDistance, semEuclDistance, 'cmap', [0 0 0]);
        axis tight
        ylim([0 0.8])
        yticks([0 0.4 0.8])
        xticks([-1 -0.5 0 0.5])
        zeroXDottedLine;
        hold off;
        legend([hl1, hl2], {'avoid - react', 'shuffle'}, 'Box', 'off');
        xlabel('time (sec)', 'FontSize', fontSz)
        ylabel('distance between PC trajectories (a.u.)', 'FontSize', fontSz)
        set(gcf,'color','w'); %set figure background white
        makepretty;
    end

end

function plotKinAvoidVsReactPcDistancesAllMiceEqualTrials(allUnitsMean, allUnitsAvoidMean, allUnitsReactMean, allUnitsMean_shuf1, allUnitsMean_shuf2, numKpts, timescaleSeg, pcaWindow, plotWindow)
% only plot Euclidean distances, with shuffle controls
% input is all N x T matrices that concatenate trial averages (over avoid / react trials) across mice/sessions
% also input numGoodUnits matrix that keeps track of # neuron splits between mice/sessions, and time windows for PCA calculation and plotting

    numMice = size(numKpts,2);
    pcaSamps = [find(pcaWindow(1)>=timescaleSeg,1,'last') find(pcaWindow(2)>=timescaleSeg,1,'last')];
    plotSamps = [find(plotWindow(1)>=timescaleSeg,1,'last') find(plotWindow(2)>=timescaleSeg,1,'last')];
    timescalePlot = timescaleSeg(plotSamps(1):plotSamps(2));
    
    usePcContrast = 0; % whether we perform contrastive PCA on data1 vs data2, or vanilla PCA on allData 
    alignPCs = 0; % whether to align all mice/sessions into common subspace for computing PCs; otherwise PCs are computed independently for each mouse/session (default for simplicity)
    smoothingSamples = [40 0]; %ex. smooth window samples 40*2.5ms frames = 100ms; for causal half-kernel use [X 0]
    %plotting constants
    fontSz = 16;
     
    savedEuclDistanceRpts = [];
    savedEuclDistanceShufRpts = [];
    numRepeats = 10; %should match number in other EqualTrials functions
    for rpt = 1:numRepeats 
        % for each sampling repeat, take get PCs and distances for each mouse
        [pcSeg, varExplained] = getKinPCs(allUnitsMean, squeeze(allUnitsAvoidMean(rpt,:,:)), squeeze(allUnitsReactMean(rpt,:,:)), numKpts, numMice, pcaSamps, usePcContrast, alignPCs); 

        currentEndUnit = 0;
        for i = 1:numMice %plot all mice separately
            currentStartUnit = currentEndUnit + 1;
            currentEndUnit = currentEndUnit + numKpts(i);

            % REAL DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % project mean avoid/ react activity for each mouse onto its respective PC matrix to get trajectories for each mouse & overlay
            avoidMeanProj = pcSeg{i}'*squeeze(allUnitsAvoidMean(rpt,currentStartUnit:currentEndUnit,plotSamps(1):plotSamps(2))); %apply PCs to the mean responses across trials
            reactMeanProj = pcSeg{i}'*squeeze(allUnitsReactMean(rpt,currentStartUnit:currentEndUnit,plotSamps(1):plotSamps(2)));
            avoidMeanProj = smoothdata(avoidMeanProj, 2, 'gaussian', smoothingSamples); %ex. smooth window 25 is 25*20ms PSTH bins = 500ms
            reactMeanProj = smoothdata(reactMeanProj, 2, 'gaussian', smoothingSamples);
            euclDistance = zeros(size(avoidMeanProj,2)-1,1);
            for dataPoint = 1:size(avoidMeanProj,2)-1
                euclDistance(dataPoint) = norm(avoidMeanProj(:,dataPoint) - reactMeanProj(:,dataPoint));
            end
            % calculate Euclidean distance between trajectories over time, normalized by sqrt(# neurons recorded), since Euclidean distance scales by the sqrt of the number of dimensions
            savedEuclDistanceRpts(rpt, i, :) = euclDistance./sqrt(numKpts(i)); %save for bounded line below
       
            % SHUFFLED DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % project mean shuffled activity for each mouse onto its respective PC matrix to get trajectories for each mouse & overlay
            shuf1MeanProj = pcSeg{i}'*squeeze(allUnitsMean_shuf1(rpt,currentStartUnit:currentEndUnit,plotSamps(1):plotSamps(2))); %option to apply same PCs as above to shuffled means
            shuf2MeanProj = pcSeg{i}'*squeeze(allUnitsMean_shuf2(rpt,currentStartUnit:currentEndUnit,plotSamps(1):plotSamps(2)));
            shuf1MeanProj = smoothdata(shuf1MeanProj, 2, 'gaussian', smoothingSamples); %ex. smooth window 25 is 25*20ms PSTH bins = 500ms
            shuf2MeanProj = smoothdata(shuf2MeanProj, 2, 'gaussian', smoothingSamples);
            euclDistanceShuf = zeros(size(shuf1MeanProj,2)-1,1);
            for dataPoint = 1:size(shuf1MeanProj,2)-1
                euclDistanceShuf(dataPoint) = norm(shuf1MeanProj(:,dataPoint) - shuf2MeanProj(:,dataPoint));
            end
            savedEuclDistanceShufRpts(rpt, i, :) = euclDistanceShuf./sqrt(numKpts(i)); %save for bounded line below

        end %end loop over mice
    end %end sampling repeats

    % now take mean across repeats for each mouse:
    savedEuclDistance = squeeze(mean(savedEuclDistanceRpts,1));
    savedEuclDistanceShuf = squeeze(mean(savedEuclDistanceShufRpts,1));

    f = figure; hold on;
%     % optionally plot individual mice to see variability:
%     for i = 1:numMiceProbes
%         plot(timescalePlot(1:end-1), savedEuclDistance(i,:), 'Color', cmapAnat, 'LineWidth', 0.1)
%         plot(timescalePlot(1:end-1), savedEuclDistanceShuf(i,:), 'Color', [0 0 0], 'LineWidth', 0.1)
%     end
    [meanEuclDistance, semEuclDistance] = grpstats(savedEuclDistance,[],{'mean' 'sem'}); %mean across columns/timepoints
    [hl1,~] = boundedline(timescalePlot(1:end-1), meanEuclDistance, semEuclDistance, 'cmap', [0.4940 0.1840 0.5560]);
    [meanEuclDistance, semEuclDistance] = grpstats(savedEuclDistanceShuf,[],{'mean' 'sem'}); %mean across columns/timepoints
    [hl2,~] = boundedline(timescalePlot(1:end-1), meanEuclDistance, semEuclDistance, 'cmap', [0 0 0]);
    axis tight
    ylim([0 1])
    xlim([-1 1])
%     xlim([timescalePlot(1)+0.5 timescalePlot(end)]) %adjust for smoothing edge effects
    %yticks([0 0.4 0.8])
    %xticks([-1 -0.5 0 0.5])
    zeroXDottedLine;
    hold off;
    legend([hl1, hl2], {'avoid - react', 'shuffle'}, 'Box', 'off');
    xlabel('time (sec)', 'FontSize', fontSz)
    ylabel('distance between PC trajectories (a.u.)', 'FontSize', fontSz)
    makepretty;
  
end

function plotKinPcaAvoidVsReactJumps(PC, idxData, tsData, kinData, titleStr) 
%PCA across only jump segments
    data = kinData.kptMat'; %use z-scored data so that large mm values don't bias
    dataTimescale = tsData.timescale(tsData.frameTimeIndeces);
    behaviorTimescale = tsData.timescale;

    %now break up projected data into avoid vs. react trials (non-opto):
    secBefore = 2;
    secAfter = 2;
    [avoidTrials, reactTrials, ~] = getAvoidVsReactJumpSegments(data, dataTimescale, behaviorTimescale, idxData, secBefore, secAfter);
    avoidTrialMean = squeeze(mean(avoidTrials,1)); %mean across trials
    reactTrialMean = squeeze(mean(reactTrials,1));
    avoidMeanProj = PC'*avoidTrialMean'; %apply PCs to the mean responses across trials
    reactMeanProj = PC'*reactTrialMean';
    avoidMeanProj = smoothdata(avoidMeanProj, 2, 'gaussian', 10); %ex. smooth window 10 is 10*50ms PSTH bins = 500ms
    reactMeanProj = smoothdata(reactMeanProj, 2, 'gaussian', 10);

%     cmap1 = autumn;
%     cmap2 = winter;
    numColorPts = 256; %standard colormap pts
    R = linspace(0.5,0.9,numColorPts)'; %go from half-black to full color
    G = linspace(0.1,0.5,numColorPts)';
    B = zeros(numColorPts,1);
    cmap1 = [R G B];
    R = zeros(numColorPts,1);
    G = linspace(0.0,0.4,numColorPts)';
    B = linspace(0.5,0.9,numColorPts)';
    cmap2 = [R G B];
    colorsIdx = floor(linspace(1, numColorPts, size(avoidMeanProj,2))); %evenly spaced colors that show time dimension

    figure;
    hold on;
    for dataPoint = 1:size(avoidMeanProj,2)-1
        if dataPoint == size(avoidMeanProj,2)-1
            psub1 = plot3(avoidMeanProj(1,dataPoint:dataPoint+1),avoidMeanProj(2,dataPoint:dataPoint+1),avoidMeanProj(3,dataPoint:dataPoint+1),'Color',cmap1(colorsIdx(dataPoint),:), 'LineWidth', 2, 'DisplayName','avoid'); %subset for legend display name
            psub2 = plot3(reactMeanProj(1,dataPoint:dataPoint+1),reactMeanProj(2,dataPoint:dataPoint+1),reactMeanProj(3,dataPoint:dataPoint+1),'Color',cmap2(colorsIdx(dataPoint),:),  'LineWidth', 2, 'DisplayName','react'); %subset for legend display name
        else
            plot3(avoidMeanProj(1,dataPoint:dataPoint+1),avoidMeanProj(2,dataPoint:dataPoint+1),avoidMeanProj(3,dataPoint:dataPoint+1),'Color',cmap1(colorsIdx(dataPoint),:), 'LineWidth', 2); 
            plot3(reactMeanProj(1,dataPoint:dataPoint+1),reactMeanProj(2,dataPoint:dataPoint+1),reactMeanProj(3,dataPoint:dataPoint+1),'Color',cmap2(colorsIdx(dataPoint),:),  'LineWidth', 2); 
        end
    end
    xlabel('PC1 (a.u.)')
    ylabel('PC2 (a.u.)')
    zlabel('PC3 (a.u.)')
    legend([psub1 psub2], 'Box', 'off');
    hold off;
    view(45,30);
    title(titleStr)
    axis tight;
%     makepretty;
    set(gcf,'color','w'); %set figure background white

end

function plotKinPcaOptoVsNotJumps(PC, idxData, tsData, kinData, titleStr) 
%PCA across only jump segments
    data = kinData.kptMat'; %use z-scored data so that large mm values don't bias
    dataTimescale = tsData.timescale(tsData.frameTimeIndeces);
    behaviorTimescale = tsData.timescale;

    %now break up projected data into opto vs. non-opto trials:
    secBefore = 2;
    secAfter = 2;
    [optoTrials, notOptoTrials, ~] = getOptoVsNotJumpSegments(data, dataTimescale, behaviorTimescale, idxData, secBefore, secAfter);
    optoTrialMean = squeeze(mean(optoTrials,1)); %mean across trials
    notOptoTrialMean = squeeze(mean(notOptoTrials,1));
    optoMeanProj = PC'*optoTrialMean'; %apply PCs to the mean responses across trials
    notOptoMeanProj = PC'*notOptoTrialMean';
    optoMeanProj = smoothdata(optoMeanProj, 2, 'gaussian', 10); %ex. smooth window 10 is 10*50ms PSTH bins = 500ms
    notOptoMeanProj = smoothdata(notOptoMeanProj, 2, 'gaussian', 10);

    numColorPts = 256; %standard colormap pts
    R = zeros(numColorPts,1); %go from half-black to full color
    G = zeros(numColorPts,1);
    B = linspace(0.5,0.9,numColorPts)';
    cmap1 = [R G B];
    R = linspace(0.0,0.4,numColorPts)';
    G = linspace(0.0,0.4,numColorPts)';
    B = linspace(0.0,0.4,numColorPts)';
    cmap2 = [R G B];
    colorsIdx = floor(linspace(1, numColorPts, size(optoMeanProj,2))); %evenly spaced colors that show time dimension

    figure;
    hold on;
    for dataPoint = 1:size(optoMeanProj,2)-1
        if dataPoint == size(optoMeanProj,2)-1
            psub1 = plot3(optoMeanProj(1,dataPoint:dataPoint+1),optoMeanProj(2,dataPoint:dataPoint+1),optoMeanProj(3,dataPoint:dataPoint+1),'Color',cmap1(colorsIdx(dataPoint),:), 'LineWidth', 2, 'DisplayName','opto'); %subset for legend display name
            psub2 = plot3(notOptoMeanProj(1,dataPoint:dataPoint+1),notOptoMeanProj(2,dataPoint:dataPoint+1),notOptoMeanProj(3,dataPoint:dataPoint+1),'Color',cmap2(colorsIdx(dataPoint),:),  'LineWidth', 2, 'DisplayName','non-opto'); %subset for legend display name
        else
            plot3(optoMeanProj(1,dataPoint:dataPoint+1),optoMeanProj(2,dataPoint:dataPoint+1),optoMeanProj(3,dataPoint:dataPoint+1),'Color',cmap1(colorsIdx(dataPoint),:), 'LineWidth', 2); 
            plot3(notOptoMeanProj(1,dataPoint:dataPoint+1),notOptoMeanProj(2,dataPoint:dataPoint+1),notOptoMeanProj(3,dataPoint:dataPoint+1),'Color',cmap2(colorsIdx(dataPoint),:),  'LineWidth', 2); 
        end
    end
    xlabel('PC1 (a.u.)')
    ylabel('PC2 (a.u.)')
    zlabel('PC3 (a.u.)')
    legend([psub1 psub2], 'Box', 'off');
    hold off;
    view(45,30);
    title(titleStr)
    axis tight;
%     makepretty;
    set(gcf,'color','w'); %set figure background white

end

function [allKptsMean, allKptsStd, timescaleSeg] = getKinMeanStdEvent(kinData, tsData, eventTimes, window)
% plot meanFR around some events either zscored or not; zscore (can lead to problems with overrepresenting small noise fluctuations in stable keypoints) is across entire session % assume data already in desired sort order
    numKpts = size(kinData.kptMat,1);
    dataTimescale = tsData.timescale(tsData.frameTimeIndeces(1:480000));
    diffDataTimescale = diff(dataTimescale);
    timeBin = mean(diffDataTimescale);
    timescaleSeg = window(1):timeBin:window(2);
    kptSamplesBefore = round(window(1)/timeBin);
    kptSamplesAfter = round(window(2)/timeBin);
    eventTimes = eventTimes(eventTimes<1200); %for datasets with optoTest at end but no video, restrict eventTimes to video recording
    thisKptXAllEvents = zeros(length(eventTimes),abs(kptSamplesBefore)+kptSamplesAfter);
    thisKptYAllEvents = zeros(length(eventTimes),abs(kptSamplesBefore)+kptSamplesAfter);
    thisKptZAllEvents = zeros(length(eventTimes),abs(kptSamplesBefore)+kptSamplesAfter);
    allKptsMean = [];
    allKptsStd = [];

    kptNames = {'tail','r_toe','r_ankle','r_knee','l_toe','l_ankle','l_knee','r_paw','l_paw','nose','mouth','chest'}; %grouped by HL, FL, other
    for kptIdx = kptNames
        for j = 1:length(eventTimes)
            kptSample = find(dataTimescale>eventTimes(j), 1, 'first');
            thisKptXAllEvents(j,:) = kinData.kptRig.(kptIdx{1})(1,kptSample+kptSamplesBefore:kptSample+kptSamplesAfter-1);
            thisKptYAllEvents(j,:) = kinData.kptRig.(kptIdx{1})(2,kptSample+kptSamplesBefore:kptSample+kptSamplesAfter-1);
            thisKptZAllEvents(j,:) = kinData.kptRig.(kptIdx{1})(3,kptSample+kptSamplesBefore:kptSample+kptSamplesAfter-1);
        end
        % center each kpt dimension for each mouse to the grand mean across trials
        thisKptsMean = [mean(thisKptXAllEvents,1) - mean(thisKptXAllEvents,'all'); mean(thisKptYAllEvents,1) - mean(thisKptYAllEvents,'all'); mean(thisKptZAllEvents,1) - mean(thisKptZAllEvents,'all')]; %take mean across trials, and subtract grand mean to center data for each mouse/kptDimension
        thisKptsStd = [std(thisKptXAllEvents,0,1); std(thisKptYAllEvents,0,1); std(thisKptZAllEvents,0,1)];
        allKptsMean = [allKptsMean; thisKptsMean];
        allKptsStd = [allKptsStd; thisKptsStd];
    end
end

function [allKptsMean, allKptsStd, timescaleSeg] = getKinMeanStdEventEqualTrials(kinData, tsData, eventTimes1, eventTimes2, window)
% get mean kpts around some events, but match # trials between event types
% NOTE that this assumes window(1) is signed appropriately
    numRepeats = 10; %since we do random sampling, repeat the process to allow averaging out sampling bias later
    for rpt = 1:numRepeats
        %use min of the two # of events to pick from each distribution
        numEvents1 = length(eventTimes1);
        numEvents2 = length(eventTimes2);
        if numEvents1 < numEvents2
            eventTimes = eventTimes1; %if less of current event types, just use all of them
        else
            s = RandStream('dsfmt19937','Seed','shuffle');
            randIdx = randperm(s, numEvents1, numEvents2); %else take random subset euqal to smaller # of trials
            eventTimes = eventTimes1(randIdx);
        end
    
        % now call original function to compute trial means:
        [allKptsMean(rpt,:,:), allKptsStd(rpt,:,:), timescaleSeg] = getKinMeanStdEvent(kinData, tsData, eventTimes, window);
    end
end

function [allKptsMeanShuf, allKptsStdShuf, timescaleSeg] = getKinMeanStdShuffle(kinData, tsData, eventTimes1, eventTimes2, useEvents1, window)
% get mean from shuffled trials (random shuffle between eventTimes1 & eventTimes2) for control
% this version draws randomly so is likely to best represent differences that occur from different # of trials of one type or the other (ex. always more react trials)
% useEvents1 decides whether to sample the # of events from eventTime1 vs eventTimes2; this allows making 2 distributions that match respective trial numbers

    %use min of the two # of events to pick from each distribution
    concatEventTimes = [eventTimes1 eventTimes2];
    numEvents = length(concatEventTimes);
    if useEvents1
        numShufEvents = length(eventTimes1); 
    else
        numShufEvents = length(eventTimes2); 
    end
    % numShufEvents = min(length(eventTimes1), length(eventTimes2)); %use min of the two #'s of events to pick from each distribution

    s = RandStream('dsfmt19937','Seed','shuffle');
    randIdx = randperm(s, numEvents, numShufEvents); %draw randomly from mix of avoid and react times
    eventTimesShuf = concatEventTimes(randIdx);

    numKpts = size(kinData.kptMat,1);
    dataTimescale = tsData.timescale(tsData.frameTimeIndeces(1:480000));
    diffDataTimescale = diff(dataTimescale);
    timeBin = mean(diffDataTimescale);
    timescaleSeg = window(1):timeBin:window(2);
    kptSamplesBefore = round(window(1)/timeBin);
    kptSamplesAfter = round(window(2)/timeBin);

    thisKptXAllEvents = zeros(length(eventTimesShuf),abs(kptSamplesBefore)+kptSamplesAfter);
    thisKptYAllEvents = zeros(length(eventTimesShuf),abs(kptSamplesBefore)+kptSamplesAfter);
    thisKptZAllEvents = zeros(length(eventTimesShuf),abs(kptSamplesBefore)+kptSamplesAfter);
    allKptsMeanShuf = [];
    allKptsStdShuf = [];

    kptNames = {'tail','r_toe','r_ankle','r_knee','l_toe','l_ankle','l_knee','r_paw','l_paw','nose','mouth','chest'}; %grouped by HL, FL, other
    for kptIdx = kptNames
        for j = 1:length(eventTimesShuf)
            kptSample = find(dataTimescale>eventTimesShuf(j), 1, 'first');
            thisKptXAllEvents(j,:) = kinData.kptRig.(kptIdx{1})(1,kptSample+kptSamplesBefore:kptSample+kptSamplesAfter-1);
            thisKptYAllEvents(j,:) = kinData.kptRig.(kptIdx{1})(2,kptSample+kptSamplesBefore:kptSample+kptSamplesAfter-1);
            thisKptZAllEvents(j,:) = kinData.kptRig.(kptIdx{1})(3,kptSample+kptSamplesBefore:kptSample+kptSamplesAfter-1);
        end

        % center each kpt dimension for each mouse to the grand mean across trials:
        thisKptsMean = [mean(thisKptXAllEvents,1) - mean(thisKptXAllEvents,'all'); mean(thisKptYAllEvents,1) - mean(thisKptYAllEvents,'all'); mean(thisKptZAllEvents,1) - mean(thisKptZAllEvents,'all')]; %take mean across trials, and subtract grand mean to center data for each mouse/kptDimension
        thisKptsStd = [std(thisKptXAllEvents,0,1); std(thisKptYAllEvents,0,1); std(thisKptZAllEvents,0,1)];
        allKptsMeanShuf = [allKptsMeanShuf; thisKptsMean];
        allKptsStdShuf = [allKptsStdShuf; thisKptsStd];
    end
end 

function [allKptsMeanShuf, allKptsStdShuf, timescaleSeg] = getKinMeanStdShuffleEqualTrials(kinData, tsData, eventTimes1, eventTimes2, window)
% get mean from shuffled trials for control
% this 'EqualTrials' version uses equal numbers of events from Times 1 & 2 to make a shuffle distribution with # trials equal to the lesser # trials 

    numShuffles = 10; 
    for sh = 1:numShuffles
        %use min of the two # of events to pick from each distribution, to match normal EqualTrials sampling procedure
        numEvents1 = length(eventTimes1);
        numEvents2 = length(eventTimes2);
        numShufEvents = min(numEvents1, numEvents2);
    
        if mod(numShufEvents,2)==0 %if even
            numEventsToChoose1 = numShufEvents/2;
            numEventsToChoose2 = numShufEvents/2;
        else %if odd
            numEventsToChoose1 = (numShufEvents+1)/2;
            numEventsToChoose2 = (numShufEvents-1)/2;
        end
    
        s = RandStream('dsfmt19937','Seed','shuffle');
        randIdx1 = randperm(s, numEvents1, numEventsToChoose1); %take equal # random trials from each type of event
        randIdx2 = randperm(s, numEvents2, numEventsToChoose2);
        eventTimesShuf = [eventTimes1(randIdx1) eventTimes2(randIdx2)]; %1/2 times from event1, 1/2 times from event2
    
%         % option to choose randomly from all trials rather than 1/2 of each type, which may unfairly mask trial varaibility in one type but doesn't seem to make a difference:
%         concatEventTimes = [eventTimes1 eventTimes2];
%         randIdx = randperm(s, length(concatEventTimes), numShufEvents);
%         eventTimesShuf = concatEventTimes(randIdx);

        % now call original function to compute trial means:
        [allKptsMeanShuf(sh,:,:), allKptsStdShuf(sh,:,:), timescaleSeg] = getKinMeanStdEvent(kinData, tsData, eventTimesShuf, window);
    end
end 

function plotKinMeanZscoreEvent(allKptsMeanZ, timescaleSeg, titleStr) 
%plot output from above
    figure;
    hold on;
    imagesc(timescaleSeg, 1:1:size(allKptsMeanZ,1), allKptsMeanZ)
    set(gca, 'YDir', 'normal');
    axis tight;
    x = [0 0]; %dotted line at event = time 0
    yLims = get(gca, 'YLim');
    y = [yLims(1) yLims(2)];
    plot(x, y, 'k--', 'LineWidth', 0.1);

    colormap(greenMag)
%     clim([-4 4]);
%     clim([0 2]);
    cb = colorbar;
    cb.Label.String = 'mean z-score position';
    ylabel('keypoint #');
    xlabel('time (sec)');
    title(['PSTH ' titleStr]);
    hold off;
    set(gcf,'color','w'); %set figure background white
end

function plotKinMeanAngleZscoreEvent(allKptsMeanZ, timescaleSeg, titleStr) 
%plot output from above
    fontSz = 16;
    figure;
    hold on;
    imagesc(timescaleSeg, 1:1:size(allKptsMeanZ,1), allKptsMeanZ)
    set(gca, 'YDir', 'normal');
    axis tight;
    zeroXDottedLine;
    set(gca, 'YTick', []);
    colormap(greenMag)
    clim([-3 3]);
    cb = colorbar;
    cb.Label.String = 'mean z-score angle';
    cb.Label.FontSize = fontSz;
    ylabel('hindlimb joints (ankle & knee, all mice)', 'FontSize', fontSz);
    xlabel('time (sec)', 'FontSize', fontSz);
    title(titleStr, 'FontSize', fontSz+2);
    hold off;
    set(gcf,'color','w'); %set figure background white
end

function plotKinMeanAngleZscoreEvent2SideBySide(kptsMeanZ1, kptsMeanZ2, timescaleSeg) 
%plot output from above
    fontSz = 16;
    figure;
    tiledlayout(1,2,'TileSpacing', 'tight')
    
    nexttile
    hold on;
    imagesc(timescaleSeg, 1:1:size(kptsMeanZ1,1), kptsMeanZ1)
    set(gca, 'YDir', 'normal');
    axis tight;
    zeroXDottedLine;
    colormap(greenMag)
    clim([-3 3]);
    ylabel('hindlimb joints (ankle & knee, all mice)', 'FontSize', fontSz);
    xlabel('time (sec)', 'FontSize', fontSz);
    title('non opto VHEx', 'FontSize', fontSz+2);
    hold off;

    nexttile
    hold on;
    imagesc(timescaleSeg, 1:1:size(kptsMeanZ2,1), kptsMeanZ2)
    set(gca, 'YDir', 'normal');
    axis tight;
    zeroXDottedLine;
    colormap(greenMag)
    clim([-3 3]);
    xlabel('time (sec)', 'FontSize', fontSz);
    title('opto - non opto', 'FontSize', fontSz+2);
    hold off;
    cb = colorbar;
    cb.Label.String = 'mean z-score angle';
    cb.Label.FontSize = fontSz;
    set(gcf,'color','w'); %set figure background white
end

function plotKinMeanAngleZscoreEvent3SideBySide(antKptsMeanZ, reactKptsMeanZ, diffKptsMeanZ, timescaleSeg) 
%plot output from above
    fontSz = 16;
    figure;
    tiledlayout(1,3,'TileSpacing', 'compact')
%     climAll = [-2 2]; %for jointAnglesNorm
    climAll = [-30 30]; %for raw angles
    numMice = size(antKptsMeanZ,1)/4; %for use with all HL angles
    
    nexttile
    hold on;
    imagesc(timescaleSeg, 1:1:size(antKptsMeanZ,1), antKptsMeanZ)
    set(gca, 'YDir', 'normal');
    axis tight;
    zeroXDottedLine;
    drawMouseSeparatorsHeatmap(ones(1,numMice).*4)
%     set(gca, 'YTick', []);
    colormap(greenMag)
    clim(climAll);
    ylabel('hindlimb joints (ankle & knee, all mice)', 'FontSize', fontSz);
%     xlabel('time (sec)', 'FontSize', fontSz);
    title('avoid only', 'FontSize', fontSz+2);
    hold off;

    nexttile
    hold on;
    imagesc(timescaleSeg, 1:1:size(reactKptsMeanZ,1), reactKptsMeanZ)
    set(gca, 'YDir', 'normal');
    axis tight;
    zeroXDottedLine;
    drawMouseSeparatorsHeatmap(ones(1,numMice).*4)
%     set(gca, 'YTick', []);
    colormap(greenMag)
    clim(climAll);
    xlabel('time (sec)', 'FontSize', fontSz);
    title('react only', 'FontSize', fontSz+2);
    hold off;

    nexttile
    hold on;
    imagesc(timescaleSeg, 1:1:size(diffKptsMeanZ,1), diffKptsMeanZ)
    set(gca, 'YDir', 'normal');
    axis tight;
    zeroXDottedLine;
    drawMouseSeparatorsHeatmap(ones(1,numMice).*4)
%     set(gca, 'YTick', []);
    colormap(greenMag)
    clim(climAll);
%     xlabel('time (sec)', 'FontSize', fontSz);
    title('avoid - react', 'FontSize', fontSz+2);
    hold off;
    cb = colorbar;
    cb.Label.String = 'mean angle';
    cb.Label.FontSize = fontSz;
    set(gcf,'color','w'); %set figure background white or black
end

function [kinParamMean, kinParamSem, kinParamTimescale] = getKinMeanSemEvent(kinDataSubset, tsData, eventTimes, window) 
% get bounded line data at events for individual kinematic parameters
% kinDataSubset should be 1D array of length(tsData.frameTimeIndeces)
    dataTimescale = tsData.timescale(tsData.frameTimeIndeces(1:480000));
    eventTimes = eventTimes(eventTimes<1200); %for datasets with optoTest at end but no video, restrict eventTimes to video recording
    diffDataTimescale = diff(dataTimescale);
    timeBin = mean(diffDataTimescale);
    kinParamTimescale = window(1):timeBin:window(2);
    kptSamplesBefore = round(window(1)/timeBin);
    kptSamplesAfter = round(window(2)/timeBin);
    thisKptAllEvents = zeros(length(eventTimes),abs(kptSamplesBefore)+kptSamplesAfter);
    for j = 1:length(eventTimes)
        kptSample = find(dataTimescale>eventTimes(j), 1, 'first');
        thisKptAllEvents(j,:) = kinDataSubset(kptSample+kptSamplesBefore:kptSample+kptSamplesAfter-1);
        thisKptAllEvents(j,:) = thisKptAllEvents(j,:) - mean(thisKptAllEvents(j,:),2); %optional centering around mean for diverging colormaps
    end   
    if length(eventTimes) > 1
        [kinParamMean, kinParamSem] = grpstats(thisKptAllEvents,[],{'mean' 'sem'}); %mean across columns/timepoints
    else
        kinParamMean = thisKptAllEvents;
        kinParamSem = zeros(1,length(kinParamTimescale));
    end
end

function [kinAnglesMean, kinAnglesSem, kinParamTimescale] = getKinMeanAnglesEvent(kinData, tsData, eventTimes, window,  kinParamName, jointNames) 
% get mean data at events for all angles
% kinDataSubset should be 1D array of length(tsData.frameTimeIndeces)
    dataTimescale = tsData.timescale(tsData.frameTimeIndeces(1:480000));
    diffDataTimescale = diff(dataTimescale);
    timeBin = mean(diffDataTimescale);
    kinParamTimescale = window(1):timeBin:window(2);
    kptSamplesBefore = round(window(1)/timeBin);
    kptSamplesAfter = round(window(2)/timeBin);
    eventTimes = eventTimes(eventTimes<1200); %for datasets with optoTest at end but no video, restrict eventTimes to video recording

    %names for accessing kinematics struct hierarchy:
    kinNames1 = kinParamName; % jointAnglesNorm, jointAngles, jointVel
    kinNames2 = {'HL'}; 
    kinNames3 = {'L', 'R'}; 
    kinNames4 = jointNames;
    
    kinAnglesMean = [];
    kinAnglesSem = [];

    for k1 = kinNames1
        for k2 = kinNames2
            for k3 = kinNames3
                for k4 = kinNames4
                    kinTrials = zeros(length(eventTimes),abs(kptSamplesBefore)+kptSamplesAfter);
                    
                    for eventIdx = 1:length(eventTimes)
                        kptSample = find(dataTimescale>eventTimes(eventIdx), 1, 'first');
                        kinTrials(eventIdx,:) = kinData.(k1{1}).(k2{1}).(k3{1}).(k4{1})(kptSample+kptSamplesBefore:kptSample+kptSamplesAfter-1);
                        kinTrials(eventIdx,:) = kinTrials(eventIdx,:) - mean(kinTrials(eventIdx,:),2); %optional centering around mean for diverging colormaps
                    end 
                    if size(kinTrials,1)==1 %if only one trial
                        kinAnglesMean = kinTrials;
                        kinAnglesSem = [];
                    else
                        [thisKinAngleMean, thisKinAngleSem] = grpstats(kinTrials,[],{'mean' 'sem'}); %mean across columns/timepoints
                        kinAnglesMean = [kinAnglesMean; thisKinAngleMean];
                        kinAnglesSem = [kinAnglesSem; thisKinAngleSem];
                    end
                end
            end
        end
    end

end

function [pts] = kptStruct2Matrix(kptStruct)
% convert 3D keypoint structure to matrix of catersian points
    %kptNames = {'tail','nose','mouth','r_paw','r_toe','r_ankle','r_knee','l_paw','l_toe','l_ankle','l_knee','chest'};
    pts(:,:,1) = kptStruct.tail;
    pts(:,:,2) = kptStruct.nose;
    pts(:,:,3) = kptStruct.mouth;
    pts(:,:,4) = kptStruct.r_paw;
    pts(:,:,5) = kptStruct.r_toe;
    pts(:,:,6) = kptStruct.r_ankle;
    pts(:,:,7) = kptStruct.r_knee;
    pts(:,:,8) = kptStruct.l_paw;
    pts(:,:,9) = kptStruct.l_toe;
    pts(:,:,10) = kptStruct.l_ankle;
    pts(:,:,11) = kptStruct.l_knee;
    pts(:,:,12) = kptStruct.chest;
end

function [kptStruct] = kptMatrix2Struct(X)
% convert 3D keypoints to interpretable structure; just hard code names instead of using eval()
    %kptNames = {'tail','nose','mouth','r_paw','r_toe','r_ankle','r_knee','l_paw','l_toe','l_ankle','l_knee','chest'};
    kptStruct.tail    = X(:,:,1);
    kptStruct.nose    = X(:,:,2);
    kptStruct.mouth   = X(:,:,3);
    kptStruct.r_paw   = X(:,:,4);
    kptStruct.r_toe   = X(:,:,5);
    kptStruct.r_ankle = X(:,:,6);
    kptStruct.r_knee  = X(:,:,7);
    kptStruct.l_paw   = X(:,:,8);
    kptStruct.l_toe   = X(:,:,9);
    kptStruct.l_ankle = X(:,:,10);
    kptStruct.l_knee  = X(:,:,11);
    kptStruct.chest   = X(:,:,12);
end

function [movementIndex] = getMovementVelocity(kinVariablesVelocityMatrix)
    % function to calculate movement across all 4 paws & mouth as an interpretable measure of gross movement
    % output is matrix of summedVelocity (normalized) x time
    kptNamesSubset = {'l_paw','r_paw','l_toe','r_toe','mouth'}; %use paws & mouth as surrogate for exploratory movements
    numKpts = length(kptNamesSubset);
    numTimepoints = size(kinVariablesVelocityMatrix.(kptNamesSubset{1}),1);
    smoothWindow = 200; %for 400Hz, 200 pts is 500ms window - this gives approximate bimodal velocity with good SNR that cooresponds to platform movement
    
    movementIndex = NaN(numTimepoints,1);
    movementIndexAllKpts = NaN(numKpts, numTimepoints);
    kpt = 1;
    for kptIdx = kptNamesSubset
        movementIndexAllKpts(kpt,:) = smoothdata(abs(kinVariablesVelocityMatrix.(kptIdx{1})), 'gaussian', smoothWindow); %take absolute value since we don't care about direction of movement here
        kpt = kpt + 1;
    end
    
    movementIndexSum = sum(movementIndexAllKpts,1); 
    movementIndex = movementIndexSum ./ (max(movementIndexSum));

end

function meanMovementIdx = getMeanMovementIndex(movementIndex, frameTimes, eventTimes , window)

    allMovements = [];
    for i = 1:length(eventTimes)
        startSample = find(frameTimes > (eventTimes(i)+window(1)), 1, 'first');
        endSample = find(frameTimes > (eventTimes(i)+window(2)), 1, 'first');
        thisMovementMean = mean(movementIndex(startSample:endSample)); %mean over window for this trial
        allMovements = [allMovements thisMovementMean];
    end
    meanMovementIdx = mean(allMovements); %mean for this mouse/session across all trials
end

function plotKinPathLengthAcrossTrials(kinVariablesTrialMatrix, trialsPerMouse)
    % function to calculate cumulative path length of a keypoint during a particular window, across trials
    % input is trials x time struct
    kptNamesSubset = {'l_paw','r_paw','l_toe','r_toe'}; %use paws as surrogate for exploratory movements
    maxBinsPlotted = 16; %standard learning curve ~180 trials
    trialsBinned = 10; % bin several trials to reduce noise
    lightGreyCmap = [0.7 0.7 0.7];

    distVsTrialsAllKpts = [];
    for kptIdx = kptNamesSubset
        % calculate Euclidean distances from X/Y/Z positions of keypoint of interest
        kptMat = kinVariablesTrialMatrix.(kptIdx{1});
        numRows = size(kptMat,1);
        numTrials = numRows/3; %always have X/Y/Z positions
        xIdx = 1:3:numRows-2;
        yIdx = 2:3:numRows-1;
        zIdx = 3:3:numRows;
        kptMatX = kptMat(xIdx,:);
        kptMatY = kptMat(yIdx,:);
        kptMatZ = kptMat(zIdx,:);
        distVsTrials = zeros(numTrials,1);
        smoothWindow = 1; %Gaussian width in frames (ex. 20 frames ~50msec @400Hz)
        for i = 1:numTrials
            currTrialDiff = [diff(kptMatX(i,:)); diff(kptMatY(i,:)); diff(kptMatZ(i,:))]; %raw X/X/Z frame-to-frame distances
            currTrialDiffSmooth = smoothdata(currTrialDiff, 2, 'gaussian', smoothWindow); %small noise in keypoints amplified by derivative, so smooth as original keypoint smoothing
            currTrialNorm = squeeze(vecnorm(currTrialDiffSmooth,2,1));  % 2-norm is the Euclidean distance frame-to-frame
            distVsTrials(i) = sum(currTrialNorm);  % integrate distance traveled over trial segment 
        end
        distVsTrialsAllKpts = [distVsTrialsAllKpts distVsTrials];
    end

    distVsTrialsAllKptsMean = mean(distVsTrialsAllKpts,2); %mean across keypoint subset
  
    % separate mice
    trialsPerMouseSeparators = cumsum(trialsPerMouse);
    numMice = length(trialsPerMouse);
    maxTrials = max(trialsPerMouse);
    distVsTrialsPerMouse = NaN(numMice,maxTrials);
    currStartTrial = 1;
    figure; hold on;
    for j = 1:numMice
        currEndTrial = trialsPerMouseSeparators(j);
        currMouseDistVsTrials = distVsTrialsAllKptsMean(currStartTrial:currEndTrial);
        distVsTrialsPerMouse(j,1:trialsPerMouse(j)) = currMouseDistVsTrials;
        currStartTrial = currEndTrial + 1;
%         plot(currMouseDistVsTrials(1:maxTrialsToPlot))
    end

    % bin across trials:
    for i = 1:numMice
        binIdx = 1;
        for j = 1:trialsBinned:size(distVsTrialsPerMouse,2)-1
            distVsTrialsPerMouseBinned(i,binIdx) = mean(distVsTrialsPerMouse(i,j:j+trialsBinned-1));
            binIdx = binIdx + 1;
        end
        plot(distVsTrialsPerMouseBinned(i,1:maxBinsPlotted),'Color', lightGreyCmap)
    end

%     [meanDist,semDist] = grpstats(distVsTrialsPerMouse,[],{'mean' 'sem'}); %stats are across columns/timepoints
%     boundedline(1:length(meanDist(1:maxTrialsToPlot)), meanDist(1:maxTrialsToPlot), semDist(1:maxTrialsToPlot), 'cmap', [0 0 0])
    [meanDist,semDist] = grpstats(distVsTrialsPerMouseBinned,[],{'mean' 'sem'}); %stats are across columns/timepoints
    boundedline(1:length(meanDist(1:maxBinsPlotted)), meanDist(1:maxBinsPlotted), semDist(1:maxBinsPlotted), 'cmap', [0 0 0], 'alpha','y', 'transparency', 0.5, 'LineWidth', 2)
    fontSz = 16;
    xlabel('trials', 'FontSize', fontSz);
    ylabel({'mean paw displacement (cm)', 'before extension'}, 'FontSize', fontSz);
    xTickSpacer = 4;
    xticks(0:xTickSpacer:maxBinsPlotted)
    xticklabels({0:trialsBinned*xTickSpacer:maxBinsPlotted*trialsBinned})
    hold off;
    axis tight;
    set(gcf,'color','w'); %set figure background white

    % all mice on 1 continuous plot option:
%     fontSz = 16;
%     figure
%     hold on;
%     plot(distVsTrials, 'b-', 'LineWidth', 1)
%     axis tight;
%     drawMouseSeparatorsPlot(trialsPerMouse)
%     xlabel('trial #', 'FontSize', fontSz);
%     ylabel('total keypoint path length @ segment', 'FontSize', fontSz);
%     hold off;
%     set(gcf,'color','w'); %set figure background white

end

function plotKinHeatmapAcrossTrials(kinVariablesTrialMatrix, timescaleSeg, kptName, trialsPerMouse, medianKpts)
    % function to plot heatmap of kinematic variable across trials/learning
    % input is trials x time struct, timescaleSeg associated

    useCenteredOnly = 1; %use grand median/mean-centered coordinates, or raw / z-score / derivative options below
%     derivativeNum = 1; % for looking at velocity, acceleration of keypoints
%     smoothWindow = 20; % smoothing for derivatives; Gaussian width in frames (ex. 20 frames ~50msec @400Hz)

    % first extract only X/Y/Z positions of keypoint of interest
    kptMat = kinVariablesTrialMatrix.(kptName);
    numRows = size(kptMat,1);
    numTrials = numRows/3; %always have X/Y/Z positions
    xIdx = 1:3:numRows-2;
    yIdx = 2:3:numRows-1;
    zIdx = 3:3:numRows;
    cLimAll = [-2 2]; %use same limits for all plots for fair visual comparison (+/- 3cm is about max head-fixed body movement)
%     cLimAll = [-0.03 0.03]; %for velocity

    if useCenteredOnly %just center all data for each kpt around grand mean for all mice
        % use median since positions are real values and resting position should be median, more resistant to outliers / position biases
%         kptMatX = kptMat(xIdx,:) - medianKpts.(kptName)(1); %raw position wrt (0,0) mean nose point, in cm; subtract grand median during ITI for better visualization with diverging colormap, but preserving trial differences
%         kptMatY = kptMat(yIdx,:) - medianKpts.(kptName)(2);
%         kptMatZ = kptMat(zIdx,:) - medianKpts.(kptName)(3);
        % can also subtract just the grand median within this time period which works better for diverging colormap but not for comparing across conditions with different medians:
        kptMatX = kptMat(xIdx,:) - median(kptMat(xIdx,:),'all'); 
        kptMatY = kptMat(yIdx,:) - median(kptMat(yIdx,:),'all');
        kptMatZ = kptMat(zIdx,:) - median(kptMat(zIdx,:),'all');
        kptMatY = -kptMatY;  %reverse / flip Y & Z so that direction matches positive down/forward for extension angle plots
        kptMatZ = -kptMatZ;
    else % 
        kptMatXRaw = kptMat(xIdx,:); %raw position wrt (0,0) mean nose point, in cm
        kptMatYRaw = kptMat(yIdx,:);
        kptMatZRaw = kptMat(zIdx,:);
        kptMatX = kptMatXRaw;
        kptMatY = -kptMatYRaw; %reverse / flip Y & Z so that direction matches positive down/forward for extension angle plots
        kptMatZ = -kptMatZRaw;    
%         kptMatX = (kptMatXRaw - mean(kptMatXRaw,2)) ./ std(kptMatXRaw,0,2); %z-score: take mean & stddev across time for each trial, but this wipes absolute position varaiance across trials
%         kptMatY = (kptMatYRaw - mean(kptMatYRaw,2)) ./ std(kptMatYRaw,0,2);
%         kptMatZ = (kptMatZRaw - mean(kptMatZRaw,2)) ./ std(kptMatZRaw,0,2);
%         kptMatX = smoothdata([diff(kptMatXRaw,derivativeNum,2) NaN(numTrials,derivativeNum)],2, 'gaussian', smoothWindow); %also try velocity or acceleration of points to find time period where movement is maximum
%         kptMatY = smoothdata([diff(kptMatYRaw,derivativeNum,2) NaN(numTrials,derivativeNum)],2, 'gaussian', smoothWindow);
%         kptMatZ = smoothdata([diff(kptMatZRaw,derivativeNum,2) NaN(numTrials,derivativeNum)],2, 'gaussian', smoothWindow);
    end

    fontSz = 16;
    fPos = figure;
    tiledlayout(1,3,'TileSpacing', 'compact')

    nexttile
    hold on;
    imagesc(timescaleSeg, 1:1:length(xIdx), kptMatX)
    set(gca, 'YDir', 'reverse');
    axis tight;
    zeroXDottedLine;
    drawMouseSeparatorsHeatmap(trialsPerMouse)
    colormap(greenMag)
    clim(cLimAll)
    cb = colorbar;
    xlim([-2 2])
    ylabel('trial #', 'FontSize', fontSz);
    xlabel('time (sec)', 'FontSize', fontSz);
    title([kptName, 'left-right'], 'FontSize', fontSz+2, 'Interpreter', 'none'); %fix interpreter for underscore
    hold off;

    nexttile
    hold on;
    imagesc(timescaleSeg, 1:1:length(yIdx), kptMatY)
    set(gca, 'YDir', 'reverse');
    axis tight;
    zeroXDottedLine;
    drawMouseSeparatorsHeatmap(trialsPerMouse)
    colormap(greenMag)
    clim(cLimAll)
    cb = colorbar;
    xlim([-2 2])
    xlabel('time (sec)', 'FontSize', fontSz);
    title(['front-back'], 'FontSize', fontSz+2);
    hold off;

    nexttile
    hold on;
    imagesc(timescaleSeg, 1:1:length(zIdx), kptMatZ)
    set(gca, 'YDir', 'reverse');
    axis tight;
    zeroXDottedLine;
    drawMouseSeparatorsHeatmap(trialsPerMouse)
    colormap(greenMag)
    clim(cLimAll)
    cb = colorbar;
    xlim([-2 2])
    cb.Label.String = 'position (cm)';
    xlabel('time (sec)', 'FontSize', fontSz);
    title(['up-down'], 'FontSize', fontSz+2);
    hold off;
    fPos.Position(3:4) = [1200 400]; %set figure size appropriate for plots
    set(gcf,'color','w'); %set figure background white

end

function [medianKpts] = getMedianItiKeypoints(kinVariablesTrialMatrixITI)
    % function to compute grand medians of keypoints across mice, during mid-ITI 'at rest', to provide a stable centering metric
    % out put is structure with just X/Y/Z medians for each keypoint
    kptNames = {'tail','r_toe','r_ankle','r_knee','l_toe','l_ankle','l_knee','r_paw','l_paw','nose','mouth','chest'}; %grouped by HL, FL, other
    numRows = size(kinVariablesTrialMatrixITI.r_toe,1);

    for kptIdx = kptNames
        kptMat = kinVariablesTrialMatrixITI.(kptIdx{1}); % get one kpt matrix at a time and take median across trials
        xIdx = 1:3:numRows-2;
        yIdx = 2:3:numRows-1;
        zIdx = 3:3:numRows;

        medianKpts.(kptIdx{1})(1) = median(kptMat(xIdx,:),'all');
        medianKpts.(kptIdx{1})(2) = median(kptMat(yIdx,:),'all');
        medianKpts.(kptIdx{1})(3) = median(kptMat(zIdx,:),'all');
    end %kpt loop 
end

function compareKinMeanHeatmaps(kinVariablesTrialMatrix1, kinVariablesTrialMatrix2, timescaleSeg, trialsPerMouse1, trialsPerMouse2, medianKpts)
    % function to plot mean heatmap comparison of kinematic variable at diferrent events
    kptNames = {'r_toe','r_ankle','r_knee','l_toe','l_ankle','l_knee'}; %,'r_paw','l_paw','nose','mouth','chest'}; %grouped by HL, FL, other
%     kptNames = {'tail'}; 
%     kptNames = {'nose','mouth'}; % use  climAll = [-0.01 0.01];
    useCenteredOnly = 1; %use grand median/mean-centered coordinates, or raw / z-score / derivative options below
    derivativeNum = 1; % for looking at velocity, acceleration of keypoints
    smoothWindow = 20; % smoothing for derivatives; Gaussian width in frames (ex. 20 frames ~50msec @400Hz)
    climAll = [-2 2]; %use same limits for all plots for fair visual comparison (+/- 2cm is about max head-fixed body movement, allows some saturation)
    kptDim = length(kptNames);
    numMice = length(trialsPerMouse1);
    kptMeansX1 = []; 
    kptMeansX2 = []; 
    kptMeansY1 = [];
    kptMeansY2 = [];
    kptMeansZ1 = [];
    kptMeansZ2 = [];
    forwardPosMean1 = NaN(numMice,length(kptNames));
    forwardPosMean2 = NaN(numMice,length(kptNames));
    blueCmap = [0 0.4470 0.7410];
    medGreyCmap = [0.4 0.4 0.4];

    currStartRow1 = 1;
    currStartRow2 = 1;
    for mouseIdx = 1:numMice %collect per mouse to see compare individual differences
        for eventType = 1:2
            kptNum=1;
            for kptIdx = kptNames
                
                if eventType == 1
                    kptMat = kinVariablesTrialMatrix1.(kptIdx{1}); % get one kpt matrix at a time and take mean across trials
                    trialsPerMouse = trialsPerMouse1;
                    currStartRow = currStartRow1;
                else
                    kptMat = kinVariablesTrialMatrix2.(kptIdx{1}); 
                    trialsPerMouse = trialsPerMouse2;
                    currStartRow = currStartRow2;
                end
                numTrialsThisMouse = trialsPerMouse(mouseIdx); 
                numRows = numTrialsThisMouse*3; %always have X/Y/Z positions
                kptMat = kptMat(currStartRow:currStartRow+numRows-1,:); %take current mouse subset
                xIdx = 1:3:numRows-2;
                yIdx = 2:3:numRows-1;
                zIdx = 3:3:numRows;
        
                if useCenteredOnly %just center all data for each kpt around grand mean for all mice
                    % use median since positions are real values and resting position should be median, more resistant to outliers / position biases
                    kptMatX = kptMat(xIdx,:) - medianKpts.(kptIdx{1})(1); %raw position wrt (0,0) mean nose point, in cm; subtract grand median during ITI for better visualization with diverging colormap, but preserving trial differences
                    kptMatY = kptMat(yIdx,:) - medianKpts.(kptIdx{1})(2);
                    kptMatZ = kptMat(zIdx,:) - medianKpts.(kptIdx{1})(3);
                    kptMatY = -kptMatY;  %reverse / flip Y & Z so that direction matches positive down/forward for extension angle plots
                    kptMatZ = -kptMatZ;
                else % 
                    kptMatXRaw = kptMat(xIdx,:); %raw position wrt (0,0) mean nose point, in cm
                    kptMatYRaw = kptMat(yIdx,:);
                    kptMatZRaw = kptMat(zIdx,:);
                    kptMatX = kptMatXRaw;
                    kptMatY = -kptMatYRaw; %reverse / flip Y & Z so that direction matches positive down/forward for extension angle plots
                    kptMatZ = -kptMatZRaw;                    
%                     kptMatX = (kptMatXRaw - mean(kptMatXRaw,2)) ./ std(kptMatXRaw,0,2); %z-score: take mean & stddev across time for each trial, but this wipes absolute position varaiance across trials
%                     kptMatY = (kptMatYRaw - mean(kptMatYRaw,2)) ./ std(kptMatYRaw,0,2);
%                     kptMatZ = (kptMatZRaw - mean(kptMatZRaw,2)) ./ std(kptMatZRaw,0,2);
%                     kptMatX = smoothdata([diff(kptMatXRaw,derivativeNum,2) NaN(numTrials,derivativeNum)],2, 'gaussian', smoothWindow); %also try velocity or acceleration of points to find time period where movement is maximum
%                     kptMatY = smoothdata([diff(kptMatYRaw,derivativeNum,2) NaN(numTrials,derivativeNum)],2, 'gaussian', smoothWindow);
%                     kptMatZ = smoothdata([diff(kptMatZRaw,derivativeNum,2) NaN(numTrials,derivativeNum)],2, 'gaussian', smoothWindow);
                end
                
                kptMeanX = mean(kptMatX,1); %mean across trials
                kptMeanY = mean(kptMatY,1);
                kptMeanZ = mean(kptMatZ,1);
                if eventType == 1
                    kptMeansX1 = [kptMeansX1; kptMeanX];
                    kptMeansY1 = [kptMeansY1; kptMeanY];
                    kptMeansZ1 = [kptMeansZ1; kptMeanZ];
                    forwardPosMean1(mouseIdx, kptNum) = mean(kptMeanY);
                else
                    kptMeansX2 = [kptMeansX2; kptMeanX];
                    kptMeansY2 = [kptMeansY2; kptMeanY];
                    kptMeansZ2 = [kptMeansZ2; kptMeanZ];
                    forwardPosMean2(mouseIdx, kptNum) = mean(kptMeanY);
                end
                kptNum=kptNum+1;
            end %kpt loop 

            if eventType == 1
                currStartRow1 = currStartRow1 + numRows;
            else
                currStartRow2 = currStartRow2 + numRows;
            end
        end %event type loop 1:2
        
    end %mouse loop

    diffKptsMeanX = kptMeansX1 - kptMeansX2;
    diffKptsMeanY = kptMeansY1 - kptMeansY2;
    diffKptsMeanZ = kptMeansZ1 - kptMeansZ2;

    fontSz = 16;
    f = figure;
    tiledlayout(3,3,'TileSpacing', 'compact')

    % X left right ------------------------------------------------------------------------
    nexttile
    hold on;
    imagesc(timescaleSeg, 1:1:size(kptMeansX1,1), kptMeansX1)
    set(gca, 'YDir', 'reverse');
    axis tight;
    zeroXDottedLine;
    drawMouseSeparatorsHeatmap(ones(1,numMice).*kptDim)
%     set(gca, 'YTick', []);
    colormap(greenMag)
    clim(climAll);
    xlim([-2 2])
    ylabel('keypoints left-right', 'FontSize', fontSz);
%     xlabel('time (sec)', 'FontSize', fontSz);
%     title('avoid only', 'FontSize', fontSz+2);
    title('not opto', 'FontSize', fontSz+2);
    hold off;

    nexttile
    hold on;
    imagesc(timescaleSeg, 1:1:size(kptMeansX2,1), kptMeansX2)
    set(gca, 'YDir', 'reverse');
    axis tight;
    zeroXDottedLine;
    drawMouseSeparatorsHeatmap(ones(1,numMice).*kptDim)
%     set(gca, 'YTick', []);
    colormap(greenMag)
    xlim([-2 2])
    clim(climAll);
    xlabel('time (sec)', 'FontSize', fontSz);
%     title('react only', 'FontSize', fontSz+2);
    title('opto', 'FontSize', fontSz+2);
    hold off;

    nexttile
    hold on;
    imagesc(timescaleSeg, 1:1:size(diffKptsMeanX,1), diffKptsMeanX)
    set(gca, 'YDir', 'reverse');
    axis tight;
    zeroXDottedLine;
    drawMouseSeparatorsHeatmap(ones(1,numMice).*kptDim)
%     set(gca, 'YTick', []);
    colormap(greenMag)
    clim(climAll);
%     xlabel('time (sec)', 'FontSize', fontSz);
%     title('avoid - react', 'FontSize', fontSz+2);
    title('not opto - opto', 'FontSize', fontSz+2);
    hold off;
    cb = colorbar;
    xlim([-2 2])
    cb.Label.String = 'mean position (cm)';
    cb.Label.FontSize = fontSz;

    % Y front back ------------------------------------------------------------------------
    nexttile
    hold on;
    imagesc(timescaleSeg, 1:1:size(kptMeansY1,1), kptMeansY1)
    set(gca, 'YDir', 'reverse');
    axis tight;
    zeroXDottedLine;
    drawMouseSeparatorsHeatmap(ones(1,numMice).*kptDim)
%     set(gca, 'YTick', []);
    colormap(greenMag)
    xlim([-2 2])
    clim(climAll);
    ylabel('keypoints front-back', 'FontSize', fontSz);
%     xlabel('time (sec)', 'FontSize', fontSz);
    hold off;

    nexttile
    hold on;
    imagesc(timescaleSeg, 1:1:size(kptMeansY2,1), kptMeansY2)
    set(gca, 'YDir', 'reverse');
    axis tight;
    zeroXDottedLine;
    drawMouseSeparatorsHeatmap(ones(1,numMice).*kptDim)
%     set(gca, 'YTick', []);
    colormap(greenMag)
    xlim([-2 2])
    clim(climAll);
    xlabel('time (sec)', 'FontSize', fontSz);
    hold off;

    nexttile
    hold on;
    imagesc(timescaleSeg, 1:1:size(diffKptsMeanY,1), diffKptsMeanY)
    set(gca, 'YDir', 'reverse');
    axis tight;
    zeroXDottedLine;
    drawMouseSeparatorsHeatmap(ones(1,numMice).*kptDim)
%     set(gca, 'YTick', []);
    colormap(greenMag)
    xlim([-2 2])
    clim(climAll);
%     xlabel('time (sec)', 'FontSize', fontSz);
    hold off;
    cb = colorbar;
    cb.Label.String = 'mean position (cm)';
    cb.Label.FontSize = fontSz;
    
    % Z up down ------------------------------------------------------------------------
    nexttile
    hold on;
    imagesc(timescaleSeg, 1:1:size(kptMeansZ1,1), kptMeansZ1)
    set(gca, 'YDir', 'reverse');
    axis tight;
    zeroXDottedLine;
    drawMouseSeparatorsHeatmap(ones(1,numMice).*kptDim)
%     set(gca, 'YTick', []);
    colormap(greenMag)
    xlim([-2 2])
    clim(climAll);
    ylabel('keypoints up-down', 'FontSize', fontSz);
%     xlabel('time (sec)', 'FontSize', fontSz);
    hold off;

    nexttile
    hold on;
    imagesc(timescaleSeg, 1:1:size(kptMeansZ2,1), kptMeansZ2)
    set(gca, 'YDir', 'reverse');
    axis tight;
    zeroXDottedLine;
    drawMouseSeparatorsHeatmap(ones(1,numMice).*kptDim)
%     set(gca, 'YTick', []);
    colormap(greenMag)
    xlim([-2 2])
    clim(climAll);
    xlabel('time (sec)', 'FontSize', fontSz);
    hold off;

    nexttile
    hold on;
    imagesc(timescaleSeg, 1:1:size(diffKptsMeanZ,1), diffKptsMeanZ)
    set(gca, 'YDir', 'reverse');
    axis tight;
    zeroXDottedLine;
    drawMouseSeparatorsHeatmap(ones(1,numMice).*kptDim)
%     set(gca, 'YTick', []);
    colormap(greenMag)
    xlim([-2 2])
    clim(climAll);
%     xlabel('time (sec)', 'FontSize', fontSz);
    hold off;
    cb = colorbar;
    cb.Label.String = 'mean position (cm)';
    cb.Label.FontSize = fontSz;

    f.Position(3:4) = [1200 1200]; %set figure size appropriate for plots
    movegui(f,'center')
    set(gcf,'color','w'); %set figure background white

    % separate plot to quantify forward position opto vs. not
    figure; 
    hold on;
    forwardPosOpto = squeeze(mean(forwardPosMean2,2)); %mean across keypoints
    forwardPosNoOpto = squeeze(mean(forwardPosMean1,2));
    
    x1 = ones(1,length(forwardPosNoOpto))*-0.5;
    x2 = ones(1,length(forwardPosOpto))*0.5;
    scatter(x1, forwardPosNoOpto, 200, medGreyCmap, 'filled');
    scatter(x2, forwardPosOpto, 200, blueCmap, 'filled');
    plot([x1; x2], [forwardPosNoOpto'; forwardPosOpto'], 'LineStyle', '-', 'Color', medGreyCmap)
    xlim([-1 1]);
    xticks([-0.5 0.5])
    xticklabels({'no opto','opto'})
    ylabel('mean forward position before extension (cm)')
    makepretty;
    
    pWilcoxSign = signrank(forwardPosOpto,forwardPosNoOpto); %Wilcoxon signed rank test for paired measurements
    p = pWilcoxSign;
    yLims = get(gca, 'YLim');
    y = [yLims(1) yLims(2)];
    if p < 0.001
       text(0, yLims(2)*0.9, '***', 'FontSize', 20);
    elseif p < 0.01
       text(0, yLims(2)*0.9, '**', 'FontSize', 20);
    elseif p < 0.05
       text(0, yLims(2)*0.9, '*', 'FontSize', 20);
    else %p > 0.05
       text(0, yLims(2)*0.9, 'n.s.', 'FontSize', 20);
    end




end

function compareKinMeanPositionsBeforeEvent(kinVariablesTrialMatrix1, kinVariablesTrialMatrix2, timescaleSeg, trialsPerMouse1, trialsPerMouse2, medianKpts)
    % function to plot mean heatmap comparison of kinematic variable at diferrent events
    kptNames = {'tail','r_toe','r_ankle','r_knee','l_toe','l_ankle','l_knee','r_paw','l_paw','nose','mouth','chest'}; %grouped by HL, FL, other
    zeroIdx = find(timescaleSeg>0,1,'first');
    useCenteredOnly = 1; %use grand median/mean-centered coordinates, or raw / z-score / derivative options below
    kptDim = length(kptNames);
    numTimesteps = length(timescaleSeg);
    numMice = length(trialsPerMouse1);
    kptMeansX1 = NaN(numMice, kptDim, numTimesteps); 
    kptMeansX2 = NaN(numMice, kptDim, numTimesteps);
    kptMeansY1 = NaN(numMice, kptDim, numTimesteps);
    kptMeansY2 = NaN(numMice, kptDim, numTimesteps);
    kptMeansZ1 = NaN(numMice, kptDim, numTimesteps);
    kptMeansZ2 = NaN(numMice, kptDim, numTimesteps);

    currStartRow1 = 1;
    currStartRow2 = 1;
    for mouseIdx = 1:numMice %collect per mouse to see compare individual differences
        for eventType = 1:2
            k = 1;
            for kptIdx = kptNames
                if eventType == 1
                    kptMat = kinVariablesTrialMatrix1.(kptIdx{1}); % get one kpt matrix at a time and take mean across trials
                    trialsPerMouse = trialsPerMouse1;
                    currStartRow = currStartRow1;
                else
                    kptMat = kinVariablesTrialMatrix2.(kptIdx{1}); 
                    trialsPerMouse = trialsPerMouse2;
                    currStartRow = currStartRow2;
                end
                numTrialsThisMouse = trialsPerMouse(mouseIdx); 
                numRows = numTrialsThisMouse*3; %always have X/Y/Z positions
                kptMat = kptMat(currStartRow:currStartRow+numRows-1,:); %take current mouse subset
                xIdx = 1:3:numRows-2;
                yIdx = 2:3:numRows-1;
                zIdx = 3:3:numRows;
        
                if useCenteredOnly %just center all data for each kpt around grand mean for all mice
                    % use median since positions are real values and resting position should be median, more resistant to outliers / position biases
                    kptMatX = kptMat(xIdx,:) - medianKpts.(kptIdx{1})(1); %raw position wrt (0,0) mean nose point, in cm; subtract grand median during ITI for better visualization with diverging colormap, but preserving trial differences
                    kptMatY = kptMat(yIdx,:) - medianKpts.(kptIdx{1})(2);
                    kptMatZ = kptMat(zIdx,:) - medianKpts.(kptIdx{1})(3);
                    kptMatY = -kptMatY;  %reverse / flip Y & Z so that direction matches positive down/forward for extension angle plots
                    kptMatZ = -kptMatZ;
                else % 
                    kptMatXRaw = kptMat(xIdx,:); %raw position wrt (0,0) mean nose point, in cm
                    kptMatYRaw = kptMat(yIdx,:);
                    kptMatZRaw = kptMat(zIdx,:);
                    kptMatX = kptMatXRaw;
                    kptMatY = -kptMatYRaw; %reverse / flip Y & Z so that direction matches positive down/forward for extension angle plots
                    kptMatZ = -kptMatZRaw;                    
%                     kptMatX = (kptMatXRaw - mean(kptMatXRaw,2)) ./ std(kptMatXRaw,0,2); %z-score: take mean & stddev across time for each trial, but this wipes absolute position varaiance across trials
%                     kptMatY = (kptMatYRaw - mean(kptMatYRaw,2)) ./ std(kptMatYRaw,0,2);
%                     kptMatZ = (kptMatZRaw - mean(kptMatZRaw,2)) ./ std(kptMatZRaw,0,2);
%                     kptMatX = smoothdata([diff(kptMatXRaw,derivativeNum,2) NaN(numTrials,derivativeNum)],2, 'gaussian', smoothWindow); %also try velocity or acceleration of points to find time period where movement is maximum
%                     kptMatY = smoothdata([diff(kptMatYRaw,derivativeNum,2) NaN(numTrials,derivativeNum)],2, 'gaussian', smoothWindow);
%                     kptMatZ = smoothdata([diff(kptMatZRaw,derivativeNum,2) NaN(numTrials,derivativeNum)],2, 'gaussian', smoothWindow);
                end
                
                kptMeanX = mean(kptMatX,1); %mean across trials
                kptMeanY = mean(kptMatY,1);
                kptMeanZ = mean(kptMatZ,1);
                if eventType == 1
                    kptMeansX1(mouseIdx,k,:) = kptMeanX;
                    kptMeansY1(mouseIdx,k,:) = kptMeanY;
                    kptMeansZ1(mouseIdx,k,:) = kptMeanZ;
                else
                    kptMeansX2(mouseIdx,k,:) = kptMeanX;
                    kptMeansY2(mouseIdx,k,:) = kptMeanY;
                    kptMeansZ2(mouseIdx,k,:) = kptMeanZ;
                end
                k = k+1;
            end %kpt loop 

            if eventType == 1
                currStartRow1 = currStartRow1 + numRows;
            else
                currStartRow2 = currStartRow2 + numRows;
            end
     
        end %event type loop 1:2
    end %mouse loop

    %now loop back though mice again to get event 1 - event 2
    diffKptsMeanX = NaN(numMice, kptDim, numTimesteps);
    diffKptsMeanY = NaN(numMice, kptDim, numTimesteps);
    diffKptsMeanZ = NaN(numMice, kptDim, numTimesteps);
    for mouseIdx = 1:numMice
        for k = 1:kptDim
            diffKptsMeanX(mouseIdx,k,:) = kptMeansX1(mouseIdx,k,:) - kptMeansX2(mouseIdx,k,:);
            diffKptsMeanY(mouseIdx,k,:) = kptMeansY1(mouseIdx,k,:) - kptMeansY2(mouseIdx,k,:);
            diffKptsMeanZ(mouseIdx,k,:) = kptMeansZ1(mouseIdx,k,:) - kptMeansZ2(mouseIdx,k,:);
        end
    end

    % plot mean in period before event as single number for all mice
    cmapX = [0.2 0.2 0.2]; % left-right
    cmapY = [0.0 0 0.8]; % front-back
    cmapZ = [0.8 0.0 0]; % up-down
    f2 = figure; hold on;
    for k = 1:kptDim
        meansX(k,:) = mean(diffKptsMeanX(:,k,1:zeroIdx),3); % mean taken from first sample to zero (i.e. at cue)
        meansY(k,:) = mean(diffKptsMeanY(:,k,1:zeroIdx),3);
        meansZ(k,:) = mean(diffKptsMeanZ(:,k,1:zeroIdx),3);
        plotSpread(meansX(k,:)', 'xValues', k+0.3,'spreadWidth', 0.2, 'spreadFcn', {'lin',[2]}, 'distributionColors', {cmapX});
        bar(k+0.3,median(meansX(k,:)), 0.2,'EdgeColor',cmapX,'FaceColor','none')
        plotSpread(meansY(k,:)', 'xValues', k+0.5,'spreadWidth', 0.2, 'spreadFcn', {'lin',[2]}, 'distributionColors', {cmapY}); 
        bar(k+0.5,median(meansY(k,:)), 0.2,'EdgeColor',cmapY,'FaceColor','none')
        plotSpread(meansZ(k,:)', 'xValues', k+0.7,'spreadWidth', 0.2, 'spreadFcn', {'lin',[2]}, 'distributionColors', {cmapZ});
        bar(k+0.7,median(meansZ(k,:)), 0.2,'EdgeColor',cmapZ,'FaceColor','none')
    end
    zeroYDottedLine
    xticks([1:kptDim]+0.5)
    xticklabels({'tail','R toe','R ankle','R knee','L toe','L ankle','L knee','R finger','L finger','nose','mouth','chest'});
    xtickangle(60)
    ylabel('mean position diff before event (cm)')
    axis([1 13 -1 1])
    hold off;
    f2.Position(3:4) = [1200 300]; %set figure size appropriate for plots
    makepretty
%     allColumns = [meansX' meansY' meansZ']; %NOTE that stats plot below will be separated by X/Y/Z coordinates and not keypoints
%     [p,tbl,stats] = kruskalwallis(allColumns); %overall p-value against null that all columns same
%     pKW = p
%     c = multcompare(stats, 'estimate', 'anova1', 'ctype', 'bonferroni'); %table of p-values across pairs

end

function plotKinAngleHeatmapAcrossTrials(kinAnglesTrialMatrix, timescaleSeg, kinParamName, jointNames, trialsPerMouse)
    % function to plot heatmap of kinematic angle/etc variable across trials/learning
    % input is trials x time struct, timescaleSeg associated

    %names for accessing kinematics struct hierarchy; must match 'getKinAnglesTrials' function definition
    kinNames1 = kinParamName; % jointAnglesNorm, jointAngles, jointVel
    kinNames2 = {'HL'}; 
    kinNames3 = {'L', 'R'}; 
    kinNames4 = jointNames;
%     cLimAll = [-3 3]; % for normalized angles
    cLimAll = [-90 90]; % for raw centered angles
    fontSz = 16;
    f = figure;
    tiledlayout(1,4,'TileSpacing', 'compact') 

    for k1 = kinNames1
        for k2 = kinNames2
            for k3 = kinNames3
                for k4 = kinNames4       
                    kptMatRaw = kinAnglesTrialMatrix.(k1{1}).(k2{1}).(k3{1}).(k4{1});
                    kptMat = kptMatRaw - mean(mean(kptMatRaw)); %subtract grand mean across all mice/trials to use diverging colormap
                    nexttile
                    hold on;
                    imagesc(timescaleSeg, 1:1:size(kptMat,1), kptMat)
                    set(gca, 'YDir', 'reverse');
                    axis tight;
                    zeroXDottedLine;
                    drawMouseSeparatorsHeatmap(trialsPerMouse)
                    colormap(greenMag)
                    clim(cLimAll)
                    cb = colorbar;
                    cb.Label.String = [kinParamName{1}, ' (a.u.)'];
                    xlabel('time (sec)', 'FontSize', fontSz);
                    ylabel('trial #', 'FontSize', fontSz);
                    title([k1{1} ' ' k2{1} ' ' k3{1} ' ' k4{1}], 'FontSize', fontSz+2);
                    hold off;
                end
            end
        end
    end
    f.Position(3:4) = [1200 400]; %set figure size appropriate for plots
    set(gcf,'color','w'); %set figure background white
end
     
function [syncKnee, syncAnkle, syncHLL, syncHLR] = plotKinAnglesCorrs(kinAnglesTrials_allMice, timescaleSeg, kinParamName, trialsPerMouse, accelNotOptoTrials_allMice)
% function to calculate L/R and knee/ankle correlations across trials:
% bilaterality knee, bilaterality ankle, knee/ankle sync L, knee/ankle sync R
% no trial binning here since we care about synchronization on individual trials, but can bin correlations after
% assumes [jointNames = {knee, ankle}] already run
    sampleTime = timescaleSeg(2) - timescaleSeg(1);
    trialSamples = length(timescaleSeg);
    kneeLtrials = kinAnglesTrials_allMice.(kinParamName{1}).HL.L.knee;
    kneeRtrials = kinAnglesTrials_allMice.(kinParamName{1}).HL.R.knee;
    ankleLtrials = kinAnglesTrials_allMice.(kinParamName{1}).HL.L.ankle;
    ankleRtrials = kinAnglesTrials_allMice.(kinParamName{1}).HL.R.ankle;

    % first separate by mice so that we can look at xcorr() and lag rather than just overall corr()
    maxTrials = max(trialsPerMouse);
    kneeLtrialsByMice = NaN(length(trialsPerMouse),maxTrials,trialSamples); % make matrix for a single correlation, for all mice, across all trial bins; pad w/ NaNs for mean calculation later
    kneeRtrialsByMice = NaN(length(trialsPerMouse),maxTrials,trialSamples);
    ankleLtrialsByMice = NaN(length(trialsPerMouse),maxTrials,trialSamples);
    ankleRtrialsByMice = NaN(length(trialsPerMouse),maxTrials,trialSamples);
    kneeXcorrMax = NaN(length(trialsPerMouse),maxTrials);
    kneeLagAtMax = NaN(length(trialsPerMouse),maxTrials);
    kneeCorrZeroLag = NaN(length(trialsPerMouse),maxTrials);
    ankleXcorrMax = NaN(length(trialsPerMouse),maxTrials);
    ankleLagAtMax = NaN(length(trialsPerMouse),maxTrials);
    ankleCorrZeroLag = NaN(length(trialsPerMouse),maxTrials);
    hindlimbLXcorrMax = NaN(length(trialsPerMouse),maxTrials);
    hindlimbLLagAtMax = NaN(length(trialsPerMouse),maxTrials);
    hindlimbLCorrZeroLag = NaN(length(trialsPerMouse),maxTrials);
    hindlimbRXcorrMax = NaN(length(trialsPerMouse),maxTrials);
    hindlimbRLagAtMax = NaN(length(trialsPerMouse),maxTrials);
    hindlimbRCorrZeroLag = NaN(length(trialsPerMouse),maxTrials);

    startTrial = 1;
    mouseNum = 1;
    maxLag = 80; % at 2.5ms frame time, 80*2.5 = 200ms is a reasonable search space for a movement about that long
    for t = trialsPerMouse
        endTrial = startTrial + t - 1;
        kneeLtrialsByMice(mouseNum, 1:t, :) = kneeLtrials(startTrial:endTrial, :); 
        kneeRtrialsByMice(mouseNum, 1:t, :) = kneeRtrials(startTrial:endTrial, :); 
        ankleLtrialsByMice(mouseNum, 1:t, :) = ankleLtrials(startTrial:endTrial, :); 
        ankleRtrialsByMice(mouseNum, 1:t, :) = ankleRtrials(startTrial:endTrial, :); 
        


        % now compute cross-correlation separately for each trial and save max and associated lag (use crosscorr() for Pearson compatibility rather than xcorr() )
        for j = 1:t
            [kneeXcorrs, kneeLags] = crosscorr(squeeze(kneeLtrialsByMice(mouseNum,j,:)),squeeze(kneeRtrialsByMice(mouseNum,j,:)),NumLags=maxLag);
            [kneeXcorrMax(mouseNum,j), lagIdx] = max(kneeXcorrs);
            kneeLagAtMax(mouseNum,j) = kneeLags(lagIdx)*sampleTime; %lag in sec (i.e how much L leads R for positive value, or vice-versa for negative)
%             %optionally look at raw trials %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             figure; subplot(2,1,1); hold on;
%             plot(squeeze(kneeLtrialsByMice(mouseNum,j,:)), 'b-')
%             plot(squeeze(kneeRtrialsByMice(mouseNum,j,:)), 'r-')
%             hold off; subplot(2,1,2);
%             plot(kneeXcorrs, 'k-')
%             keyboard
%             % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [ankleXcorrs, ankleLags] = crosscorr(squeeze(ankleLtrialsByMice(mouseNum,j,:)),squeeze(ankleRtrialsByMice(mouseNum,j,:)),NumLags=maxLag);
            [ankleXcorrMax(mouseNum,j), lagIdx] = max(ankleXcorrs);
            ankleLagAtMax(mouseNum,j) = ankleLags(lagIdx)*sampleTime; %lag in sec (i.e how much L leads R for positive value, or vice-versa for negative)

            [hindlimbLXcorrs, hindlimbLLags] = crosscorr(squeeze(kneeLtrialsByMice(mouseNum,j,:)),squeeze(ankleLtrialsByMice(mouseNum,j,:)),NumLags=maxLag);
            [hindlimbLXcorrMax(mouseNum,j), lagIdx] = max(hindlimbLXcorrs);
            hindlimbLLagAtMax(mouseNum,j) = hindlimbLLags(lagIdx)*sampleTime; %lag in sec (i.e how much knee leads ankle for positive value, or vice-versa for negative)

            [hindlimbRXcorrs, hindlimbRLags] = crosscorr(squeeze(kneeRtrialsByMice(mouseNum,j,:)),squeeze(ankleRtrialsByMice(mouseNum,j,:)),NumLags=maxLag);
            [hindlimbRXcorrMax(mouseNum,j), lagIdx] = max(hindlimbRXcorrs);
            hindlimbRLagAtMax(mouseNum,j) = hindlimbRLags(lagIdx)*sampleTime; %lag in sec (i.e how much knee leads ankle for positive value, or vice-versa for negative)
        end

        % standard Pearson correlation at zero lag for all trials at once:
        kneeCorrZeroLag(mouseNum,1:t) = diag(corr(kneeLtrials(startTrial:endTrial, :)',kneeRtrials(startTrial:endTrial, :)'));
        ankleCorrZeroLag(mouseNum,1:t) = diag(corr(ankleLtrials(startTrial:endTrial, :)',ankleRtrials(startTrial:endTrial, :)'));
        hindlimbLCorrZeroLag(mouseNum,1:t) = diag(corr(kneeLtrials(startTrial:endTrial, :)',ankleLtrials(startTrial:endTrial, :)'));
        hindlimbRCorrZeroLag(mouseNum,1:t) = diag(corr(kneeRtrials(startTrial:endTrial, :)',ankleRtrials(startTrial:endTrial, :)'));

        startTrial = endTrial + 1;
        mouseNum = mouseNum + 1;
    end

    % also calculate overall sum of all 4 measured hindlimb joints, as this should scale with overal force production
    allHindlimbAnglesMat = kneeLtrials+kneeRtrials+ankleLtrials+ankleRtrials;
    maxHindlimbAnglesVsTrials = max(allHindlimbAnglesMat,[],2); % take the max of the angle sums over the time period before extension, which will be highest if all are synced and should correlate with platform accel

    fontSz = 16;
    f = figure;
    tiledlayout(4,3,'TileSpacing', 'compact');
    cmap1 = [0 0 0];
    cmap2 = [0 0 0.5];
    cmap3 = [0.5 0 0];
    maxTrialsPlotted = 51; %standard learning curve 151
    yMin = -0.5;
    lWidth = 3;
    smoothingSamples = 10; % window of trials to smooth over to see trend better
    lineCmaps = linspecer(length(trialsPerMouse));

    for i = 1:length(trialsPerMouse) %overlay plots of for each mouse
        cmapLine = lineCmaps(i,:);
        nexttile(1); hold on;
        plot(smoothdata(kneeCorrZeroLag(i,1:maxTrialsPlotted),2,'gaussian',smoothingSamples), 'Color', cmapLine, 'LineStyle', '-', LineWidth=0.1);
        nexttile(2); hold on;
        plot(smoothdata(kneeXcorrMax(i,1:maxTrialsPlotted),2,'gaussian',smoothingSamples), 'Color', cmap2, 'LineStyle', '--', LineWidth=0.1);
        nexttile(3); hold on;
        plot(smoothdata(kneeLagAtMax(i,1:maxTrialsPlotted).*1000,2,'gaussian',smoothingSamples), 'Color', cmap3, 'LineStyle', '--', LineWidth=0.1);

        nexttile(4); hold on;
        plot(smoothdata(ankleCorrZeroLag(i,1:maxTrialsPlotted),2,'gaussian',smoothingSamples), 'Color', cmapLine, 'LineStyle', '-', LineWidth=0.1);
        nexttile(5); hold on;
        plot(smoothdata(ankleXcorrMax(i,1:maxTrialsPlotted),2,'gaussian',smoothingSamples), 'Color', cmap2, 'LineStyle', '--', LineWidth=0.1);
        nexttile(6); hold on;
        plot(smoothdata(ankleLagAtMax(i,1:maxTrialsPlotted).*1000,2,'gaussian',smoothingSamples), 'Color', cmap3, 'LineStyle', '--', LineWidth=0.1);

        nexttile(7); hold on;
        plot(smoothdata(hindlimbLCorrZeroLag(i,1:maxTrialsPlotted),2,'gaussian',smoothingSamples), 'Color', cmapLine, 'LineStyle', '-', LineWidth=0.1);
        nexttile(8); hold on;
        plot(smoothdata(hindlimbLXcorrMax(i,1:maxTrialsPlotted),2,'gaussian',smoothingSamples), 'Color', cmap2, 'LineStyle', '--', LineWidth=0.1);
        nexttile(9); hold on;
        plot(smoothdata(hindlimbLLagAtMax(i,1:maxTrialsPlotted).*1000,2,'gaussian',smoothingSamples), 'Color', cmap3, 'LineStyle', '--', LineWidth=0.1);

        nexttile(10); hold on;
        plot(smoothdata(hindlimbRCorrZeroLag(i,1:maxTrialsPlotted),2,'gaussian',smoothingSamples), 'Color', cmapLine, 'LineStyle', '-', LineWidth=0.1);
        nexttile(11); hold on;
        plot(smoothdata(hindlimbRXcorrMax(i,1:maxTrialsPlotted),2,'gaussian',smoothingSamples), 'Color', cmap2, 'LineStyle', '--', LineWidth=0.1);
        nexttile(12); hold on;
        plot(smoothdata(hindlimbRLagAtMax(i,1:maxTrialsPlotted).*1000,2,'gaussian',smoothingSamples), 'Color', cmap3, 'LineStyle', '--', LineWidth=0.1);
    end

    % finish each plot separately
    nexttile(1)
    kneeCorrZeroLag = smoothdata(kneeCorrZeroLag, 2, 'gaussian', smoothingSamples); 
    [meankneeCorrZeroLag, semkneeCorrZeroLag] = grpstats(kneeCorrZeroLag,[],{'mean' 'sem'}); %mean across trials
    boundedline(1:maxTrialsPlotted, meankneeCorrZeroLag(1:maxTrialsPlotted), semkneeCorrZeroLag(1:maxTrialsPlotted), 'cmap', cmap1);
%     meanKnee1 = squeeze(mean(kneeCorrZeroLag,1,'omitnan')); %smooth over trials a bit and take mean across mice
%     plot(meanKnee1(1:maxTrialsPlotted), 'Color', cmap1, LineWidth=lWidth); 
    axis([1 maxTrialsPlotted yMin 1])
    ylabel('L-R knee', 'FontSize', fontSz);
    title('zero lag corr', 'FontSize', fontSz);
    hold off;

    nexttile(2)
    kneeXcorrMax = smoothdata(kneeXcorrMax, 2, 'gaussian', smoothingSamples); 
    [meankneeXcorrMax, semkneeXcorrMax] = grpstats(kneeXcorrMax,[],{'mean' 'sem'}); %mean across trials
    boundedline(1:maxTrialsPlotted, meankneeXcorrMax(1:maxTrialsPlotted), semkneeXcorrMax(1:maxTrialsPlotted), 'cmap', cmap2);
%     meanKnee2 = squeeze(mean(kneeXcorrMax,1,'omitnan')); 
%     plot(meanKnee2(1:maxTrialsPlotted), 'Color', cmap2, LineWidth=lWidth); 
    axis([1 maxTrialsPlotted yMin 1])
    title('xcorr max', 'FontSize', fontSz);
    hold off;

    nexttile(3)
    kneeLagAtMax = smoothdata(kneeLagAtMax, 2, 'gaussian', smoothingSamples); 
    [meankneeLagAtMax, semkneeLagAtMax] = grpstats(kneeLagAtMax,[],{'mean' 'sem'}); %mean across trials
    boundedline(1:maxTrialsPlotted, meankneeLagAtMax(1:maxTrialsPlotted), semkneeLagAtMax(1:maxTrialsPlotted), 'cmap', cmap3);
%     meanKnee3 = squeeze(mean(kneeLagAtMax.*1000,1,'omitnan')); 
%     plot(meanKnee3(1:maxTrialsPlotted), 'Color', cmap3, LineWidth=lWidth); 
    axis([1 maxTrialsPlotted -200 200])
    title('lag at max (ms)', 'FontSize', fontSz);
    hold off;

    nexttile(4)
    ankleCorrZeroLag = smoothdata(ankleCorrZeroLag, 2, 'gaussian', smoothingSamples); 
    [meanankleCorrZeroLag, semankleCorrZeroLag] = grpstats(ankleCorrZeroLag,[],{'mean' 'sem'}); %mean across trials
    boundedline(1:maxTrialsPlotted, meanankleCorrZeroLag(1:maxTrialsPlotted), semankleCorrZeroLag(1:maxTrialsPlotted), 'cmap', cmap1);
%     meanAnkle1 = squeeze(mean(ankleCorrZeroLag,1,'omitnan')); %take mean across mice
%     plot(meanAnkle1(1:maxTrialsPlotted), 'Color', cmap1, LineWidth=lWidth); 
    axis([1 maxTrialsPlotted yMin 1])
    ylabel('L-R ankle', 'FontSize', fontSz);
    hold off;

    nexttile(5)
    ankleXcorrMax = smoothdata(ankleXcorrMax, 2, 'gaussian', smoothingSamples); 
    [meanankleXcorrMax, semankleXcorrMax] = grpstats(ankleXcorrMax,[],{'mean' 'sem'}); %mean across trials
    boundedline(1:maxTrialsPlotted, meanankleXcorrMax(1:maxTrialsPlotted), semankleXcorrMax(1:maxTrialsPlotted), 'cmap', cmap2);
%     meanAnkle2 = squeeze(mean(ankleXcorrMax,1,'omitnan')); 
%     plot(meanAnkle2(1:maxTrialsPlotted), 'Color', cmap2, LineWidth=lWidth); 
    axis([1 maxTrialsPlotted yMin 1])
    hold off;

    nexttile(6)
    ankleLagAtMax = smoothdata(ankleLagAtMax, 2, 'gaussian', smoothingSamples); 
    [meanankleLagAtMax, semankleLagAtMax] = grpstats(ankleLagAtMax,[],{'mean' 'sem'}); %mean across trials
    boundedline(1:maxTrialsPlotted, meanankleLagAtMax(1:maxTrialsPlotted), semankleLagAtMax(1:maxTrialsPlotted), 'cmap', cmap3);
%     meanAnkle3 = squeeze(mean(ankleLagAtMax.*1000,1,'omitnan')); 
%     plot(meanAnkle3(1:maxTrialsPlotted), 'Color', cmap3, LineWidth=lWidth); 
    axis([1 maxTrialsPlotted -200 200])
    hold off;

    nexttile(7)
    hindlimbLCorrZeroLag = smoothdata(hindlimbLCorrZeroLag, 2, 'gaussian', smoothingSamples); 
    [meanhindlimbLCorrZeroLag, semhindlimbLCorrZeroLag] = grpstats(hindlimbLCorrZeroLag,[],{'mean' 'sem'}); %mean across trials
    boundedline(1:maxTrialsPlotted, meanhindlimbLCorrZeroLag(1:maxTrialsPlotted), semhindlimbLCorrZeroLag(1:maxTrialsPlotted), 'cmap', cmap1);
%     meanHindlimbL1 = squeeze(mean(hindlimbLCorrZeroLag,1,'omitnan')); %take mean across mice
%     plot(meanHindlimbL1(1:maxTrialsPlotted), 'Color', cmap1, LineWidth=lWidth); 
    axis([1 maxTrialsPlotted yMin 1])
    ylabel('L knee-ankle', 'FontSize', fontSz);
    hold off;

    nexttile(8)
    hindlimbLXcorrMax = smoothdata(hindlimbLXcorrMax, 2, 'gaussian', smoothingSamples); 
    [meanhindlimbLXcorrMax, semhindlimbLXcorrMax] = grpstats(hindlimbLXcorrMax,[],{'mean' 'sem'}); %mean across trials
    boundedline(1:maxTrialsPlotted, meanhindlimbLXcorrMax(1:maxTrialsPlotted), semhindlimbLXcorrMax(1:maxTrialsPlotted), 'cmap', cmap2);
%     meanHindlimbL2 = squeeze(mean(hindlimbLXcorrMax,1,'omitnan')); 
%     plot(meanHindlimbL2(1:maxTrialsPlotted), 'Color', cmap2, LineWidth=lWidth); 
    axis([1 maxTrialsPlotted yMin 1])
    hold off;

    nexttile(9)
    hindlimbLLagAtMax = smoothdata(hindlimbLLagAtMax, 2, 'gaussian', smoothingSamples); 
    [meanhindlimbLLagAtMax, semhindlimbLLagAtMax] = grpstats(hindlimbLLagAtMax,[],{'mean' 'sem'}); %mean across trials
    boundedline(1:maxTrialsPlotted, meanhindlimbLLagAtMax(1:maxTrialsPlotted), semhindlimbLLagAtMax(1:maxTrialsPlotted), 'cmap', cmap3);
%     meanHindlimbL3 = squeeze(mean(hindlimbLLagAtMax.*1000,1,'omitnan')); 
%     plot(meanHindlimbL3(1:maxTrialsPlotted), 'Color', cmap3, LineWidth=lWidth); 
    axis([1 maxTrialsPlotted -200 200])
    hold off;

    nexttile(10)
    hindlimbRCorrZeroLag = smoothdata(hindlimbRCorrZeroLag, 2, 'gaussian', smoothingSamples); 
    [meanhindlimbRCorrZeroLag, semhindlimbRCorrZeroLag] = grpstats(hindlimbRCorrZeroLag,[],{'mean' 'sem'}); %mean across trials
    boundedline(1:maxTrialsPlotted, meanhindlimbRCorrZeroLag(1:maxTrialsPlotted), semhindlimbRCorrZeroLag(1:maxTrialsPlotted), 'cmap', cmap1);
%     meanHindlimbR1 = squeeze(mean(hindlimbRCorrZeroLag,1,'omitnan')); %take mean across mice
%     plot(meanHindlimbR1(1:maxTrialsPlotted), 'Color', cmap1, LineWidth=lWidth); 
    axis([1 maxTrialsPlotted yMin 1])
    xlabel('trials', 'FontSize', fontSz);
    ylabel('R knee-ankle', 'FontSize', fontSz);
    hold off;

    nexttile(11)
    hindlimbRXcorrMax = smoothdata(hindlimbRXcorrMax, 2, 'gaussian', smoothingSamples); 
    [meanhindlimbRXcorrMax, semhindlimbRXcorrMax] = grpstats(hindlimbRXcorrMax,[],{'mean' 'sem'}); %mean across trials
    boundedline(1:maxTrialsPlotted, meanhindlimbRXcorrMax(1:maxTrialsPlotted), semhindlimbRXcorrMax(1:maxTrialsPlotted), 'cmap', cmap2);
%     meanHindlimbR2 = squeeze(mean(hindlimbRXcorrMax,1,'omitnan')); 
%     plot(meanHindlimbR2(1:maxTrialsPlotted), 'Color', cmap2, LineWidth=lWidth); 
    axis([1 maxTrialsPlotted yMin 1])
    xlabel('trials', 'FontSize', fontSz);
    hold off;

    nexttile(12)
    hindlimbRLagAtMax = smoothdata(hindlimbRLagAtMax, 2, 'gaussian', smoothingSamples); 
    [meanhindlimbRLagAtMax, semhindlimbRLagAtMax] = grpstats(hindlimbRLagAtMax,[],{'mean' 'sem'}); %mean across trials
    boundedline(1:maxTrialsPlotted, meanhindlimbRLagAtMax(1:maxTrialsPlotted), semhindlimbRLagAtMax(1:maxTrialsPlotted), 'cmap', cmap3);
%     meanHindlimbR3 = squeeze(mean(hindlimbRLagAtMax.*1000,1,'omitnan')); 
%     plot(meanHindlimbR3(1:maxTrialsPlotted), 'Color', cmap3, LineWidth=lWidth); 
    axis([1 maxTrialsPlotted -200 200])
    xlabel('trials', 'FontSize', fontSz);
    hold off;

    set(gcf,'color','w'); %set figure background white

    % also scatterplot hindlimb angles vs platform accel
%     figure;
%     hold on,
%     scatter(accelNotOptoTrials_allMice, maxHindlimbAnglesVsTrials)
%     axis tight;
%     xlabel('trial peak accel. (cm/sec^2)', 'FontSize', fontSz);
%     ylabel('trial hindlimb angles max. sum', 'FontSize', fontSz);
%     hold off;
%     set(gcf,'color','w'); %set figure background white
%     rAngleVsAccel = corrcoef(accelNotOptoTrials_allMice, maxHindlimbAnglesVsTrials)

    % return correlation measures for comparison with ephys:
    syncKnee = kneeCorrZeroLag;
    syncAnkle = ankleCorrZeroLag;
    syncHLL = hindlimbLCorrZeroLag;
    syncHLR = hindlimbRCorrZeroLag;

end

function trialSimilarity = plotKinCorrMatrixAcrossTrials(kinVariablesTrialMatrix, kptName, trialsPerMouse)
    % function to plot how well trials correlate with each other (i.e. correlation matrix) wrt kinematic variable
    % return matrix diags(mouseNum, EuclideanDimension, TrialBins) for kptName, to plot all together elsewhere

    % first extract only X/Y/Z positions of keypoint of interest
    kptMat = kinVariablesTrialMatrix.(kptName);
    numRows = size(kptMat,1);
    numTrials = numRows/3; %always have X/Y/Z positions
    xIdx = 1:3:numRows-2;
    yIdx = 2:3:numRows-1;
    zIdx = 3:3:numRows;
    trialsBinned = 10; % bin several trials to reduce noise
    plotCorr = 0; %turn this off if only getting trialSimilarity across all keypoints for comparison plot
    maxBinsPlot = 30;

    %use normalized positions for correlation, with diverging colormap to more clearly visualize differences
    kptMatXRaw = kptMat(xIdx,:); %raw position wrt (0,0) mean nose point, in cm
    kptMatYRaw = kptMat(yIdx,:);
    kptMatZRaw = kptMat(zIdx,:);
%     kptMatX = (kptMatXRaw - mean(kptMatXRaw,2)) ./ std(kptMatXRaw,0,2); %take mean & stddev across time
%     kptMatY = (kptMatYRaw - mean(kptMatYRaw,2)) ./ std(kptMatYRaw,0,2);
%     kptMatZ = (kptMatZRaw - mean(kptMatZRaw,2)) ./ std(kptMatZRaw,0,2);
    kptMatX = kptMatXRaw;
    kptMatY = kptMatYRaw;
    kptMatZ = kptMatZRaw;
    numTimepoints = size(kptMat,2);

    % bin across trials:
    binIdx = 1;
    startTrial = 1;
    binsPerMouse = [];
    for i = 1:length(trialsPerMouse) %need to deal with boundaries across mice
        endTrial = startTrial + trialsPerMouse(i) - trialsBinned;
        startBin = binIdx;
        for j = startTrial:trialsBinned:endTrial
%                 kptMatXBinned(binIdx,:) = mean(kptMatX(j:j+trialsBinned-1,:),1); 
%                 kptMatYBinned(binIdx,:) = mean(kptMatY(j:j+trialsBinned-1,:),1); 
%                 kptMatZBinned(binIdx,:) = mean(kptMatZ(j:j+trialsBinned-1,:),1); 
                kptMatXBinned(binIdx,:) = median(kptMatX(j:j+trialsBinned-1,:),1); %option to use median instead of mean, but doesn't seem to matter much
                kptMatYBinned(binIdx,:) = median(kptMatY(j:j+trialsBinned-1,:),1); 
%                 kptMatZBinned(binIdx,:) = median(kptMatZ(j:j+trialsBinned-1,:),1); 

                % OPTION to plot model of motor skill learning instead:
                kptMatZBinned(binIdx,:) = (-0.5:(1/numTimepoints):0.5-(1/numTimepoints)).*(j) + rand(1,numTimepoints).*((endTrial-j).^1.05); 

                binIdx = binIdx + 1;

        end
        currentMouseBins = binIdx - startBin;
        unbinnedTrials = trialsPerMouse(i) - (j+trialsBinned-1);
        startTrial = (j+trialsBinned) + unbinnedTrials;
        binsPerMouse = [binsPerMouse currentMouseBins];
%         binsPerMouse = [binsPerMouse maxBinsPlot]; %option to plot with maxBins for calculation example
    end

    [rhoX, pvalX] = corr(kptMatXBinned');
    [rhoY, pvalY] = corr(kptMatYBinned');
    [rhoZ, pvalZ] = corr(kptMatZBinned');
%     [rhoX, pvalX] = corr(kptMatXBinned(1:maxBinsPlot,:)');
%     [rhoY, pvalY] = corr(kptMatYBinned(1:maxBinsPlot,:)');
%     [rhoZ, pvalZ] = corr(kptMatZBinned(1:maxBinsPlot,:)');
    
    if plotCorr
%         alpha = 0.05; %option to only plot significant values
%         Xsig = find(pvalX>alpha);
%         Ysig = find(pvalY>alpha);
%         Zsig = find(pvalZ>alpha);

        Xsig = 1:length(pvalX); % for OPTION to plot model of motor skill learning instead:
        Ysig = 1:length(pvalY);
        Zsig = 1:length(pvalZ);
        %figure; imagesc(kptMatZBinned)
        

        fontSz = 16;
        climAll = [-1 1];
        fCorr = figure;
        tiledlayout(1,3,'TileSpacing', 'compact')
    
        nexttile
        hold on;
        imx = imagesc(rhoX);
        sigTemp = ones(size(rhoX));
        sigTemp(Xsig) = zeros(length(Xsig),1); 
        set(imx, 'AlphaData', sigTemp);    %set insignificant values white
        set(gca, 'YDir', 'reverse');
        axis image;
        colormap("parula")
        clim(climAll);
        drawMouseSeparatorsCorr(binsPerMouse);
        ylabel('trial bins');
        xlabel('trial bins');
        title([kptName, ' med-lat'], 'FontSize', fontSz+2, 'Interpreter', 'none'); %fix interpreter for underscore
        hold off;
        
        nexttile
        hold on;
        imy = imagesc(rhoY);
        sigTemp = ones(size(rhoY));
        sigTemp(Ysig) = zeros(length(Ysig),1); 
        set(imy, 'AlphaData', sigTemp);    %set insignificant values white
        set(gca, 'YDir', 'reverse');
        axis image;
        colormap("parula")
        clim(climAll);
        drawMouseSeparatorsCorr(binsPerMouse);
        ylabel('trial bins');
        xlabel('trial bins');
        title(['front-back'], 'FontSize', fontSz+2, 'Interpreter', 'none'); %fix interpreter for underscore
        hold off;
    
        nexttile
        hold on;
        imz = imagesc(rhoZ);
        sigTemp = ones(size(rhoZ));
%         sigTemp(Zsig) = zeros(length(Zsig),1); 

        for z = 1:size(rhoZ,1)
            sigTemp(z,z) = 0; % for OPTION to plot model of motor skill learning instead:
        end


        set(imz, 'AlphaData', sigTemp);    %set insignificant values white
        set(gca, 'YDir', 'reverse');
        axis image;
        colormap("parula")
        clim(climAll);
%         clim([0 1]); % for OPTION to plot model of motor skill learning instead
        drawMouseSeparatorsCorr(binsPerMouse);
        cb = colorbar;
        cb.Label.String = 'trial correlation';
        ylabel('trial bins');
        xlabel('trial bins');
        title(['vertical'], 'FontSize', fontSz+2, 'Interpreter', 'none'); %fix interpreter for underscore
        hold off;
        
        fCorr.Position(3:4) = [1200 400]; %set figure size appropriate for plots
        set(gcf,'color','w'); %set figure background white
    end

    % return 1-off diagonals [diag(X,1)] as a measure of each trial bin correlation with next trial bin, 
    % which should increase over time (slope != 0) for behavior that becomes more self similar
    % OR, optionally use the last column that reflects correlation of the last trial bin with all others, but do this 
    % for each mouse individual matrix so compute in loop below (however, this is biased towards tendency of temporally 
    % contigous trials to have similar posture, so correlations always increase towards last bin)
    useDiag = 1; 
    if useDiag
        similarityX = diag(rhoX,1);
        similarityY = diag(rhoY,1);
        similarityZ = diag(rhoZ,1);
    end
    
    maxBins = max(binsPerMouse);
    trialSimilarity = NaN(length(binsPerMouse),3,maxBins); % make matrix with X/Y/Z dimensions for a single keypoint, for each mouse, across all trial bins; pad w/ NaNs for mean calculation later
    startBin = 1;
    mouseNum = 1;
    for b = binsPerMouse
         endBin = startBin + b - 2; % skip last correlation of bins across mice
        if ~useDiag
            similarityX = rhoX(endBin,startBin:endBin);
            similarityY = rhoY(endBin,startBin:endBin);
            similarityZ = rhoZ(endBin,startBin:endBin);
            trialSimilarity(mouseNum, 1, 1:b-1) = similarityX; 
            trialSimilarity(mouseNum, 2, 1:b-1) = similarityY;
            trialSimilarity(mouseNum, 3, 1:b-1) = similarityZ;
        else
            %for use with diagonal measure:
            trialSimilarity(mouseNum, 1, 1:b-1) = similarityX(startBin:endBin); 
            trialSimilarity(mouseNum, 2, 1:b-1) = similarityY(startBin:endBin);
            trialSimilarity(mouseNum, 3, 1:b-1) = similarityZ(startBin:endBin);
        end
        startBin = endBin + 2;
        mouseNum = mouseNum + 1;
    end


end

function plotKinCorrSimilarities(kptTrialSimilarities,kptNames)
% plot all diagonals as a measure of self-similarity of kinematics across trials 
% input trialSimilarity is structure with (mouseNum, EuclideanDimension, TrialBins) for each kptName
    f1 = figure;
    tiledlayout('flow','TileSpacing', 'compact')
    fontSz = 18;
    cmapX = [0.2 0.2 0.2]; % left-right
    cmapY = [0.0 0 0.8]; % front-back
    cmapZ = [0.8 0.0 0]; % up-down
    maxBinsPlotted = 6; %standard learning curve ~15 bins
    kptNum = 1;
    numKpts = length(kptNames);
    

    for kptIdx = kptNames %loop separate subplot for each kpt containing X/Y/Z corrs in different colors
        diagsMat = kptTrialSimilarities.(kptIdx{1}); % single matrix diags(mouseNum, EuclideanDimension, TrialBins) for kptName
        
        nexttile
        hold on;
        for i = 1:size(diagsMat,1) %loop for separate corr lines for each mouse
            currMouseLR = squeeze(diagsMat(i,1,1:maxBinsPlotted));
            currMouseFB = squeeze(diagsMat(i,2,1:maxBinsPlotted));
            currMouseUD = squeeze(diagsMat(i,3,1:maxBinsPlotted));
            plot(currMouseLR,'Color',cmapX, 'LineStyle', '--', LineWidth=0.1) %left-right
            plot(currMouseFB,'Color',cmapY, 'LineStyle', '--', LineWidth=0.1) %front-back
            plot(currMouseUD,'Color',cmapZ, 'LineStyle', '--', LineWidth=0.1) %vertical / up-down
            % store mean value and slope for each mouse
            magsLR{kptNum}(i) = mean(currMouseLR,'omitnan');
            slopesLR{kptNum}(i) = mean(diff(currMouseLR),'omitnan');
            magsFB{kptNum}(i) = mean(currMouseFB,'omitnan');
            slopesFB{kptNum}(i) = mean(diff(currMouseFB),'omitnan');
            magsUD{kptNum}(i) = mean(currMouseUD,'omitnan');
            slopesUD{kptNum}(i) = mean(diff(currMouseUD),'omitnan');
        end

        meanX = squeeze(mean(diagsMat(:,1,:),1,'omitnan'));
        meanY = squeeze(mean(diagsMat(:,2,:),1,'omitnan'));
        meanZ = squeeze(mean(diagsMat(:,3,:),1,'omitnan'));
        psubX = plot(meanX(1:maxBinsPlotted), 'Color', cmapX, LineWidth=5); %medial-lateral
        psubY = plot(meanY(1:maxBinsPlotted), 'Color', cmapY, LineWidth=5); %front-back
        psubZ =plot(meanZ(1:maxBinsPlotted), 'Color', cmapZ, LineWidth=5); %vertical
        axis([1 maxBinsPlotted -1 1])
        ylabel([kptIdx{1} ' corr.'], 'Interpreter','none', 'FontSize', fontSz);
        xlabel('trial bins', 'FontSize', fontSz);
        hold off;
        kptNum = kptNum + 1;
    end
    
    f1.Position(3:4) = [1200 1200]; %set figure size appropriate for plots
    legend([psubX psubY psubZ], {'medial-lateral','front-back','vertical'}, 'location', 'best')
    set(gcf,'color','w'); %set figure background white

    % also make dot plots for mean slopes and correlation maginitudes across keypoints
    %plotSpread() needs a cell array with different column vectors as input
    numMice = size(diagsMat,1);
    f2 = figure; hold on;
    for k = 1:numKpts
        plotSpread(magsLR, 'xValues', [1:numKpts]+0.3,'spreadWidth', 0.2, 'spreadFcn', {'lin',[2]}, 'distributionColors', {cmapX});
        bar([1:numKpts]+0.3,mean(reshape(cell2mat(magsLR),numMice,numKpts),1), 0.2,'EdgeColor',cmapX,'FaceColor','none')
        plotSpread(magsFB, 'xValues', [1:numKpts]+0.5,'spreadWidth', 0.2, 'spreadFcn', {'lin',[2]}, 'distributionColors', {cmapY}); 
        bar([1:numKpts]+0.5,mean(reshape(cell2mat(magsFB),numMice,numKpts),1), 0.2,'EdgeColor',cmapY,'FaceColor','none')
        plotSpread(magsUD, 'xValues', [1:numKpts]+0.7,'spreadWidth', 0.2, 'spreadFcn', {'lin',[2]}, 'distributionColors', {cmapZ});
        bar([1:numKpts]+0.7,mean(reshape(cell2mat(magsUD),numMice,numKpts),1), 0.2,'EdgeColor',cmapZ,'FaceColor','none')
    end
    xticks([1:numKpts]+0.5)
    xticklabels({'tail','R toe','R ankle','R knee','L toe','L ankle','L knee','R finger','L finger','nose','mouth','chest'});
    xtickangle(60)
    ylabel('mean mag. of corr. across trials')
    axis([1 13 -1 1])
    hold off;
    f2.Position(3:4) = [1200 300]; %set figure size appropriate for plots
    set(gcf,'color','w'); %set figure background white
    makepretty
%     allColumns = [reshape(cell2mat(magsLR),numMice,numKpts) reshape(cell2mat(magsFB),numMice,numKpts) reshape(cell2mat(magsUD),numMice,numKpts)];
%     [pKW,tblKW,statsKW] = kruskalwallis(allColumns); %overall p-value against null that all columns same
%     cKW = multcompare(statsKW, 'estimate', 'kruskalwallis', 'ctype', 'dunn-sidak'); %table of p-values across pairs

    f3 = figure; hold on;
    for k = 1:numKpts
        plotSpread(slopesLR, 'xValues', [1:numKpts]+0.3,'spreadWidth', 0.2, 'spreadFcn', {'lin',[2]}, 'distributionColors', {cmapX});
        bar([1:numKpts]+0.3,mean(reshape(cell2mat(slopesLR),numMice,numKpts),1), 0.2,'EdgeColor',cmapX,'FaceColor','none')
        plotSpread(slopesFB, 'xValues', [1:numKpts]+0.5,'spreadWidth', 0.2, 'spreadFcn', {'lin',[2]}, 'distributionColors', {cmapY}); 
        bar([1:numKpts]+0.5,mean(reshape(cell2mat(slopesFB),numMice,numKpts),1), 0.2,'EdgeColor',cmapY,'FaceColor','none')
        plotSpread(slopesUD, 'xValues', [1:numKpts]+0.7,'spreadWidth', 0.2, 'spreadFcn', {'lin',[2]}, 'distributionColors', {cmapZ});
        bar([1:numKpts]+0.7,mean(reshape(cell2mat(slopesUD),numMice,numKpts),1), 0.2,'EdgeColor',cmapZ,'FaceColor','none')
    end
    xticks([1:numKpts]+0.5)
    xticklabels({'tail','R toe','R ankle','R knee','L toe','L ankle','L knee','R finger','L finger','nose','mouth','chest'});
    xtickangle(60)
    ylabel('mean slope of corr. across trials')
    axis([1 13 -1 1])
    hold off;
    f3.Position(3:4) = [1200 300]; %set figure size appropriate for plots
    set(gcf,'color','w'); %set figure background white
    makepretty
%     allColumns = [reshape(cell2mat(slopesLR),numMice,numKpts) reshape(cell2mat(slopesFB),numMice,numKpts) reshape(cell2mat(slopesUD),numMice,numKpts)];
%     [pKW,tblKW,statsKW] = kruskalwallis(allColumns); %overall p-value against null that all columns same
%     cKW = multcompare(statsKW, 'estimate', 'kruskalwallis', 'ctype', 'dunn-sidak'); %table of p-values across pairs
end

%% ephys functions
function probeData = initializeEphyData(anatNames)
% function to initialize and organize arrays for ephys analyses
    probeData.allUnitsZFR = []; % collect all (not just mean) for comprehensive rastermap sort
    probeData.allUnitsMeanZFR = []; % for mean across all VHEx
    probeData.allUnitsMeanZFR_shuf1 = []; % make shuffled versions with randomly selected avoid vs react trials, for control comparions
    probeData.allUnitsMeanZFR_shuf2 = [];
    probeData.allUnitsAtCueMeanZFR_shuf1 = [];
    probeData.allUnitsAtCueMeanZFR_shuf2 = [];
    probeData.diffMeanZFR = []; % avoid - react @ VHEx
    probeData.avoidMeanZFR = [];
    probeData.reactMeanZFR = [];
    probeData.optoMeanZFR = []; % all VHEx opto trials 
    probeData.avoidAtCueMeanZFR = [];
    probeData.reactAtCueMeanZFR = [];
    probeData.avoidAtItiMeanZFR = [];
    probeData.reactAtItiMeanZFR = [];
    probeData.optoAtCueMeanZFR = [];
    probeData.beepOnMeanZFR = [];
    probeData.beepOffMeanZFR = [];
    probeData.optoOnMeanZFR = [];
    probeData.platformReleaseMeanZFR = [];
    probeData.platformCommandMeanZFR = [];
    probeData.itiMeanZFR = [];
    probeData.coldAirOnMeanZFR = [];
    probeData.optoJumpMeanZFR = [];
    probeData.avoidMinusOptoJump = [];
    probeData.reactMinusOptoJump = [];
    probeData.numGoodUnits = []; % collect number of units per probe to be able to search later using seDepths / anatHiers                    
    probeData.stDepths = []; % collect probe depths for each unit in stCell 
    probeData.anatHiers = []; % collect anatomical hierarchy for each unit in stCell for later selecting subsets
    probeData.LFP = [];

    % also make sub-structures for specified anatomical subsets; this uses more memory than indexing later, but simplifies per trial analyses
    for anatIdx = anatNames
        probeData.(anatIdx{1}).allUnitsMeanZFR = []; % for mean across all VHEx
        probeData.(anatIdx{1}).allUnitsMeanZFR_shuf1 = []; % make 2 shuffled versions with randomly selected avoid vs react trials, for control comparions
        probeData.(anatIdx{1}).allUnitsMeanZFR_shuf2 = [];
        probeData.(anatIdx{1}).allUnitsAtCueMeanZFR_shuf1 = [];
        probeData.(anatIdx{1}).allUnitsAtCueMeanZFR_shuf2 = [];
        probeData.(anatIdx{1}).allUnitsAvoidMeanZFR = [];
        probeData.(anatIdx{1}).allUnitsReactMeanZFR = [];
        probeData.(anatIdx{1}).allUnitsOptoMeanZFR = [];
        probeData.(anatIdx{1}).allUnitsPlatformReleaseMeanZFR = [];
        probeData.(anatIdx{1}).allUnitsPlatformCommandMeanZFR = [];
        probeData.(anatIdx{1}).itiMeanZFR = [];
        probeData.(anatIdx{1}).allUnitsAvoidAtCueMeanZFR = [];
        probeData.(anatIdx{1}).allUnitsReactAtCueMeanZFR = [];
        probeData.(anatIdx{1}).allUnitsAvoidAtItiMeanZFR = [];
        probeData.(anatIdx{1}).allUnitsReactAtItiMeanZFR = [];
        probeData.(anatIdx{1}).allUnitsOptoAtCueMeanZFR = [];
        probeData.(anatIdx{1}).allUnitsDiffMeanZFR = []; % avoid - react @ VHEx
        probeData.(anatIdx{1}).allUnitsColdAirOnMeanZFR = [];
        probeData.(anatIdx{1}).allUnitsBeepOnAvoidMeanZFR = [];
        probeData.(anatIdx{1}).allUnitsBeepOnReactMeanZFR = [];
        probeData.(anatIdx{1}).allUnitsOptoOnMeanZFR = [];
        probeData.(anatIdx{1}).numGoodUnits = [];
        probeData.(anatIdx{1}).stDepths = [];
        probeData.(anatIdx{1}).anatHiers = [];
        probeData.(anatIdx{1}).platDistZAvoid = []; %for platform decoders
        probeData.(anatIdx{1}).platDistZReact = [];
        probeData.(anatIdx{1}).platDistZOpto = [];
        probeData.(anatIdx{1}).platDistZIti = [];
        probeData.(anatIdx{1}).allUnitsTrialsAvoid = [];
        probeData.(anatIdx{1}).allUnitsTrialsReact = [];
        probeData.(anatIdx{1}).allUnitsTrialsOpto = [];
        probeData.(anatIdx{1}).allUnitsTrialsIti = [];
        probeData.(anatIdx{1}).allUnitsTrialsAvoidAtCue = []; 
        probeData.(anatIdx{1}).allUnitsTrialsReactAtCue = []; 
        probeData.(anatIdx{1}).allUnitsTrialsOptoAtCue = []; 
        probeData.(anatIdx{1}).allUnitsTrialsColdAir = [];
        probeData.(anatIdx{1}).allUnitsTrialsShuf1 = []; 
        probeData.(anatIdx{1}).allUnitsTrialsShuf2 = []; 
        probeData.(anatIdx{1}).allUnitsTrialsShufAtCue1 = []; 
        probeData.(anatIdx{1}).allUnitsTrialsShufAtCue2 = [];

        if strcmp(anatIdx{1}, 'CTX') %for cortex also make a non-cortex subset for comparison
            probeData.nonCTX.allUnitsMeanZFR = []; % for mean across all VHEx
            probeData.nonCTX.allUnitsMeanZFR_shuf1 = []; % make 2 shuffled versions with randomly selected avoid vs react trials, for control comparions
            probeData.nonCTX.allUnitsMeanZFR_shuf2 = [];
            probeData.nonCTX.allUnitsAtCueMeanZFR_shuf1 = [];
            probeData.nonCTX.allUnitsAtCueMeanZFR_shuf2 = [];
            probeData.nonCTX.allUnitsAvoidMeanZFR = [];
            probeData.nonCTX.allUnitsReactMeanZFR = [];
            probeData.nonCTX.allUnitsOptoMeanZFR = [];
            probeData.nonCTX.allUnitsPlatformReleaseMeanZFR = [];
            probeData.nonCTX.allUnitsPlatformCommandMeanZFR = [];
            probeData.nonCTX.itiMeanZFR = [];
            probeData.nonCTX.allUnitsAvoidAtCueMeanZFR = [];
            probeData.nonCTX.allUnitsReactAtCueMeanZFR = [];
            probeData.nonCTX.allUnitsAvoidAtItiMeanZFR = [];
            probeData.nonCTX.allUnitsReactAtItiMeanZFR = [];
            probeData.nonCTX.allUnitsOptoAtCueMeanZFR = [];
            probeData.nonCTX.allUnitsDiffMeanZFR = []; % avoid - react
            probeData.nonCTX.allUnitsColdAirOnMeanZFR = [];
            probeData.nonCTX.allUnitsBeepOnAvoidMeanZFR = [];
            probeData.nonCTX.allUnitsBeepOnReactMeanZFR = [];
            probeData.nonCTX.allUnitsOptoOnMeanZFR = [];
            probeData.nonCTX.numGoodUnits = [];
            probeData.nonCTX.stDepths = [];
            probeData.nonCTX.anatHiers = [];
            probeData.nonCTX.platDistZAvoid = []; 
            probeData.nonCTX.platDistZReact = [];
            probeData.nonCTX.platDistZOpto = [];
            probeData.nonCTX.platDistZIti = [];
            probeData.nonCTX.allUnitsTrialsAvoid = [];
            probeData.nonCTX.allUnitsTrialsReact = [];
            probeData.nonCTX.allUnitsTrialsOpto = [];
            probeData.nonCTX.allUnitsTrialsIti = [];
            probeData.nonCTX.allUnitsTrialsAvoidAtCue = []; 
            probeData.nonCTX.allUnitsTrialsReactAtCue = []; 
            probeData.nonCTX.allUnitsTrialsOptoAtCue = []; 
            probeData.nonCTX.allUnitsTrialsColdAir = [];
            probeData.nonCTX.allUnitsTrialsShuf1 = []; 
            probeData.nonCTX.allUnitsTrialsShuf2 = []; 
            probeData.nonCTX.allUnitsTrialsShufAtCue1 = []; 
            probeData.nonCTX.allUnitsTrialsShufAtCue2 = [];
        end
    end
end

function probeData = concatEphyData(probeDataBefore, currentProbeData, anatNames)
% concatenate probes across different mice
    probeData.allUnitsZFR = [probeDataBefore.allUnitsZFR; currentProbeData.allUnitsZFR]; % collect all (not just mean) for comprehensive rastermap sort
    probeData.allUnitsMeanZFR = [probeDataBefore.allUnitsMeanZFR; currentProbeData.allUnitsMeanZFR]; % for mean across all VHEx
    probeData.allUnitsMeanZFR_shuf1 = [probeDataBefore.allUnitsMeanZFR_shuf1; currentProbeData.allUnitsMeanZFR_shuf1]; % make 2 shuffled versions with randomly selected avoid vs react trials, for control comparions
    probeData.allUnitsMeanZFR_shuf2 = [probeDataBefore.allUnitsMeanZFR_shuf2; currentProbeData.allUnitsMeanZFR_shuf2];
    probeData.allUnitsAtCueMeanZFR_shuf1 = [probeDataBefore.allUnitsAtCueMeanZFR_shuf1; currentProbeData.allUnitsAtCueMeanZFR_shuf1]; % make 2 shuffled versions with randomly selected avoid vs react trials, for control comparions
    probeData.allUnitsAtCueMeanZFR_shuf2 = [probeDataBefore.allUnitsAtCueMeanZFR_shuf2; currentProbeData.allUnitsAtCueMeanZFR_shuf2];
    probeData.diffMeanZFR = [probeDataBefore.diffMeanZFR; currentProbeData.diffMeanZFR]; % avoid - react
    probeData.avoidMeanZFR = [probeDataBefore.avoidMeanZFR; currentProbeData.avoidMeanZFR];
    probeData.reactMeanZFR = [probeDataBefore.reactMeanZFR; currentProbeData.reactMeanZFR];
    probeData.optoMeanZFR = [probeDataBefore.optoMeanZFR; currentProbeData.optoMeanZFR];
    probeData.avoidAtCueMeanZFR = [probeDataBefore.avoidAtCueMeanZFR; currentProbeData.avoidAtCueMeanZFR];
    probeData.reactAtCueMeanZFR = [probeDataBefore.reactAtCueMeanZFR; currentProbeData.reactAtCueMeanZFR];
    probeData.avoidAtItiMeanZFR = [probeDataBefore.avoidAtItiMeanZFR; currentProbeData.avoidAtItiMeanZFR];
    probeData.reactAtItiMeanZFR = [probeDataBefore.reactAtItiMeanZFR; currentProbeData.reactAtItiMeanZFR];
    probeData.optoAtCueMeanZFR = [probeDataBefore.optoAtCueMeanZFR; currentProbeData.optoAtCueMeanZFR];
    probeData.beepOnMeanZFR = [probeDataBefore.beepOnMeanZFR; currentProbeData.beepOnMeanZFR];
    probeData.beepOffMeanZFR = [probeDataBefore.beepOffMeanZFR; currentProbeData.beepOffMeanZFR];
    probeData.optoOnMeanZFR = [probeDataBefore.optoOnMeanZFR; currentProbeData.optoOnMeanZFR];
    probeData.platformReleaseMeanZFR = [probeDataBefore.platformReleaseMeanZFR; currentProbeData.platformReleaseMeanZFR];
    probeData.platformCommandMeanZFR = [probeDataBefore.platformCommandMeanZFR; currentProbeData.platformCommandMeanZFR];
    probeData.itiMeanZFR = [probeDataBefore.itiMeanZFR; currentProbeData.itiMeanZFR];
    probeData.coldAirOnMeanZFR = [probeDataBefore.coldAirOnMeanZFR; currentProbeData.coldAirOnMeanZFR];
    probeData.optoJumpMeanZFR = [probeDataBefore.optoJumpMeanZFR; currentProbeData.optoJumpMeanZFR];
    probeData.avoidMinusOptoJump = [probeDataBefore.avoidMinusOptoJump; currentProbeData.avoidMinusOptoJump];
    probeData.reactMinusOptoJump = [probeDataBefore.reactMinusOptoJump; currentProbeData.reactMinusOptoJump];
    probeData.numGoodUnits = [probeDataBefore.numGoodUnits currentProbeData.numGoodUnits]; % collect number of units per probe to be able to search later using seDepths / anatHiers                    
    probeData.stDepths = [probeDataBefore.stDepths; currentProbeData.stDepths]; % collect probe depths for each unit in stCell 
    probeData.anatHiers = [probeDataBefore.anatHiers currentProbeData.anatHiers]; % collect anatomical hierarchy for each unit in stCell for later selecting subsets
    probeData.LFP = [probeDataBefore.LFP currentProbeData.LFP];

    % also make sub-structures for specified anatomical subsets; this uses more memory than indexing later, but simplifies per trial analyses
    for anatIdx = anatNames
        probeData.(anatIdx{1}).allUnitsMeanZFR = [probeDataBefore.(anatIdx{1}).allUnitsMeanZFR; currentProbeData.(anatIdx{1}).allUnitsMeanZFR]; % for mean across all VHEx

        %for EqualTrials sampling repeats, concatenate over units dimension
        probeData.(anatIdx{1}).allUnitsMeanZFR_shuf1 = cat(2, probeDataBefore.(anatIdx{1}).allUnitsMeanZFR_shuf1, currentProbeData.(anatIdx{1}).allUnitsMeanZFR_shuf1); % 
        probeData.(anatIdx{1}).allUnitsMeanZFR_shuf2 = cat(2, probeDataBefore.(anatIdx{1}).allUnitsMeanZFR_shuf2, currentProbeData.(anatIdx{1}).allUnitsMeanZFR_shuf2);
        probeData.(anatIdx{1}).allUnitsAtCueMeanZFR_shuf1 = cat(2, probeDataBefore.(anatIdx{1}).allUnitsAtCueMeanZFR_shuf1, currentProbeData.(anatIdx{1}).allUnitsAtCueMeanZFR_shuf1); % 
        probeData.(anatIdx{1}).allUnitsAtCueMeanZFR_shuf2 = cat(2, probeDataBefore.(anatIdx{1}).allUnitsAtCueMeanZFR_shuf2, currentProbeData.(anatIdx{1}).allUnitsAtCueMeanZFR_shuf2);
        probeData.(anatIdx{1}).allUnitsAvoidMeanZFR = cat(2, probeDataBefore.(anatIdx{1}).allUnitsAvoidMeanZFR, currentProbeData.(anatIdx{1}).allUnitsAvoidMeanZFR);
        probeData.(anatIdx{1}).allUnitsReactMeanZFR = cat(2, probeDataBefore.(anatIdx{1}).allUnitsReactMeanZFR, currentProbeData.(anatIdx{1}).allUnitsReactMeanZFR);
        probeData.(anatIdx{1}).allUnitsOptoMeanZFR = cat(2, probeDataBefore.(anatIdx{1}).allUnitsOptoMeanZFR, currentProbeData.(anatIdx{1}).allUnitsOptoMeanZFR);
        probeData.(anatIdx{1}).allUnitsAvoidAtCueMeanZFR = cat(2, probeDataBefore.(anatIdx{1}).allUnitsAvoidAtCueMeanZFR, currentProbeData.(anatIdx{1}).allUnitsAvoidAtCueMeanZFR);
        probeData.(anatIdx{1}).allUnitsReactAtCueMeanZFR = cat(2, probeDataBefore.(anatIdx{1}).allUnitsReactAtCueMeanZFR, currentProbeData.(anatIdx{1}).allUnitsReactAtCueMeanZFR);
        probeData.(anatIdx{1}).allUnitsOptoAtCueMeanZFR = cat(2, probeDataBefore.(anatIdx{1}).allUnitsOptoAtCueMeanZFR, currentProbeData.(anatIdx{1}).allUnitsOptoAtCueMeanZFR);

        probeData.(anatIdx{1}).allUnitsPlatformReleaseMeanZFR = [probeDataBefore.(anatIdx{1}).allUnitsPlatformReleaseMeanZFR; currentProbeData.(anatIdx{1}).allUnitsPlatformReleaseMeanZFR];
        probeData.(anatIdx{1}).allUnitsPlatformCommandMeanZFR = [probeDataBefore.(anatIdx{1}).allUnitsPlatformCommandMeanZFR; currentProbeData.(anatIdx{1}).allUnitsPlatformCommandMeanZFR];
        probeData.(anatIdx{1}).itiMeanZFR = [probeDataBefore.(anatIdx{1}).itiMeanZFR; currentProbeData.(anatIdx{1}).itiMeanZFR];
        probeData.(anatIdx{1}).allUnitsAvoidAtItiMeanZFR = [probeDataBefore.(anatIdx{1}).allUnitsAvoidAtItiMeanZFR; currentProbeData.(anatIdx{1}).allUnitsAvoidAtItiMeanZFR];
        probeData.(anatIdx{1}).allUnitsReactAtItiMeanZFR = [probeDataBefore.(anatIdx{1}).allUnitsReactAtItiMeanZFR; currentProbeData.(anatIdx{1}).allUnitsReactAtItiMeanZFR];
        probeData.(anatIdx{1}).allUnitsDiffMeanZFR = [probeDataBefore.(anatIdx{1}).allUnitsDiffMeanZFR; currentProbeData.(anatIdx{1}).allUnitsDiffMeanZFR]; % avoid - react
        probeData.(anatIdx{1}).allUnitsColdAirOnMeanZFR = [probeDataBefore.(anatIdx{1}).allUnitsColdAirOnMeanZFR; currentProbeData.(anatIdx{1}).allUnitsColdAirOnMeanZFR];
        probeData.(anatIdx{1}).allUnitsBeepOnAvoidMeanZFR = [probeDataBefore.(anatIdx{1}).allUnitsBeepOnAvoidMeanZFR; currentProbeData.(anatIdx{1}).allUnitsBeepOnAvoidMeanZFR];
        probeData.(anatIdx{1}).allUnitsBeepOnReactMeanZFR = [probeDataBefore.(anatIdx{1}).allUnitsBeepOnReactMeanZFR; currentProbeData.(anatIdx{1}).allUnitsBeepOnReactMeanZFR];
        probeData.(anatIdx{1}).allUnitsOptoOnMeanZFR = [probeDataBefore.(anatIdx{1}).allUnitsOptoOnMeanZFR; currentProbeData.(anatIdx{1}).allUnitsOptoOnMeanZFR];
        probeData.(anatIdx{1}).numGoodUnits = [probeDataBefore.(anatIdx{1}).numGoodUnits currentProbeData.(anatIdx{1}).numGoodUnits];
        probeData.(anatIdx{1}).stDepths = [probeDataBefore.(anatIdx{1}).stDepths; currentProbeData.(anatIdx{1}).stDepths];
        probeData.(anatIdx{1}).anatHiers = [probeDataBefore.(anatIdx{1}).anatHiers currentProbeData.(anatIdx{1}).anatHiers];

        probeData.(anatIdx{1}).platDistZAvoid = [probeDataBefore.(anatIdx{1}).platDistZAvoid {currentProbeData.(anatIdx{1}).platDistZAvoid}]; 
        probeData.(anatIdx{1}).platDistZReact = [probeDataBefore.(anatIdx{1}).platDistZReact {currentProbeData.(anatIdx{1}).platDistZReact}]; 
        probeData.(anatIdx{1}).platDistZOpto = [probeDataBefore.(anatIdx{1}).platDistZOpto {currentProbeData.(anatIdx{1}).platDistZOpto}]; 
        probeData.(anatIdx{1}).platDistZIti = [probeDataBefore.(anatIdx{1}).platDistZIti {currentProbeData.(anatIdx{1}).platDistZIti}]; 
        probeData.(anatIdx{1}).allUnitsTrialsAvoid = [probeDataBefore.(anatIdx{1}).allUnitsTrialsAvoid {currentProbeData.(anatIdx{1}).allUnitsTrialsAvoid}]; 
        probeData.(anatIdx{1}).allUnitsTrialsReact = [probeDataBefore.(anatIdx{1}).allUnitsTrialsReact {currentProbeData.(anatIdx{1}).allUnitsTrialsReact}]; 
        probeData.(anatIdx{1}).allUnitsTrialsOpto = [probeDataBefore.(anatIdx{1}).allUnitsTrialsOpto {currentProbeData.(anatIdx{1}).allUnitsTrialsOpto}]; 
        probeData.(anatIdx{1}).allUnitsTrialsIti = [probeDataBefore.(anatIdx{1}).allUnitsTrialsIti {currentProbeData.(anatIdx{1}).allUnitsTrialsIti}]; 
        probeData.(anatIdx{1}).allUnitsTrialsAvoidAtCue = [probeDataBefore.(anatIdx{1}).allUnitsTrialsAvoidAtCue {currentProbeData.(anatIdx{1}).allUnitsTrialsAvoidAtCue}]; 
        probeData.(anatIdx{1}).allUnitsTrialsReactAtCue = [probeDataBefore.(anatIdx{1}).allUnitsTrialsReactAtCue {currentProbeData.(anatIdx{1}).allUnitsTrialsReactAtCue}];  
        probeData.(anatIdx{1}).allUnitsTrialsOptoAtCue = [probeDataBefore.(anatIdx{1}).allUnitsTrialsOptoAtCue {currentProbeData.(anatIdx{1}).allUnitsTrialsOptoAtCue}]; 
        probeData.(anatIdx{1}).allUnitsTrialsColdAir = [probeDataBefore.(anatIdx{1}).allUnitsTrialsColdAir {currentProbeData.(anatIdx{1}).allUnitsTrialsColdAir}];
        probeData.(anatIdx{1}).allUnitsTrialsShuf1 = [probeDataBefore.(anatIdx{1}).allUnitsTrialsShuf1 {currentProbeData.(anatIdx{1}).allUnitsTrialsShuf1}];  
        probeData.(anatIdx{1}).allUnitsTrialsShuf2 = [probeDataBefore.(anatIdx{1}).allUnitsTrialsShuf2 {currentProbeData.(anatIdx{1}).allUnitsTrialsShuf2}];  
        probeData.(anatIdx{1}).allUnitsTrialsShufAtCue1 = [probeDataBefore.(anatIdx{1}).allUnitsTrialsShufAtCue1 {currentProbeData.(anatIdx{1}).allUnitsTrialsShufAtCue1}]; 
        probeData.(anatIdx{1}).allUnitsTrialsShufAtCue2 = [probeDataBefore.(anatIdx{1}).allUnitsTrialsShufAtCue2 {currentProbeData.(anatIdx{1}).allUnitsTrialsShufAtCue2}]; 

        if strcmp(anatIdx{1}, 'CTX') %for cortex also make a non-cortex subset for comparison
            probeData.nonCTX.allUnitsMeanZFR = [probeDataBefore.nonCTX.allUnitsMeanZFR; currentProbeData.nonCTX.allUnitsMeanZFR]; % for mean across all VHEx

            %for EqualTrials sampling repeats, concatenate over units dimension
            probeData.nonCTX.allUnitsMeanZFR_shuf1 = cat(2, probeDataBefore.nonCTX.allUnitsMeanZFR_shuf1, currentProbeData.nonCTX.allUnitsMeanZFR_shuf1); % 
            probeData.nonCTX.allUnitsMeanZFR_shuf2 = cat(2, probeDataBefore.nonCTX.allUnitsMeanZFR_shuf2, currentProbeData.nonCTX.allUnitsMeanZFR_shuf2);
            probeData.nonCTX.allUnitsAtCueMeanZFR_shuf1 = cat(2, probeDataBefore.nonCTX.allUnitsAtCueMeanZFR_shuf1, currentProbeData.nonCTX.allUnitsAtCueMeanZFR_shuf1); % 
            probeData.nonCTX.allUnitsAtCueMeanZFR_shuf2 = cat(2, probeDataBefore.nonCTX.allUnitsAtCueMeanZFR_shuf2, currentProbeData.nonCTX.allUnitsAtCueMeanZFR_shuf2);
            probeData.nonCTX.allUnitsAvoidMeanZFR = cat(2, probeDataBefore.nonCTX.allUnitsAvoidMeanZFR, currentProbeData.nonCTX.allUnitsAvoidMeanZFR);
            probeData.nonCTX.allUnitsReactMeanZFR = cat(2, probeDataBefore.nonCTX.allUnitsReactMeanZFR, currentProbeData.nonCTX.allUnitsReactMeanZFR);
            probeData.nonCTX.allUnitsOptoMeanZFR = cat(2, probeDataBefore.nonCTX.allUnitsOptoMeanZFR, currentProbeData.nonCTX.allUnitsOptoMeanZFR);
            probeData.nonCTX.allUnitsAvoidAtCueMeanZFR = cat(2, probeDataBefore.nonCTX.allUnitsAvoidAtCueMeanZFR, currentProbeData.nonCTX.allUnitsAvoidAtCueMeanZFR);
            probeData.nonCTX.allUnitsReactAtCueMeanZFR = cat(2, probeDataBefore.nonCTX.allUnitsReactAtCueMeanZFR, currentProbeData.nonCTX.allUnitsReactAtCueMeanZFR);
            probeData.nonCTX.allUnitsOptoAtCueMeanZFR = cat(2, probeDataBefore.nonCTX.allUnitsOptoAtCueMeanZFR, currentProbeData.nonCTX.allUnitsOptoAtCueMeanZFR);
            
            probeData.nonCTX.allUnitsPlatformReleaseMeanZFR = [probeDataBefore.nonCTX.allUnitsPlatformReleaseMeanZFR; currentProbeData.nonCTX.allUnitsPlatformReleaseMeanZFR];
            probeData.nonCTX.allUnitsPlatformCommandMeanZFR = [probeDataBefore.nonCTX.allUnitsPlatformCommandMeanZFR; currentProbeData.nonCTX.allUnitsPlatformCommandMeanZFR];
            probeData.nonCTX.itiMeanZFR = [probeDataBefore.nonCTX.itiMeanZFR; currentProbeData.nonCTX.itiMeanZFR];
            probeData.nonCTX.allUnitsAvoidAtItiMeanZFR = [probeDataBefore.nonCTX.allUnitsAvoidAtItiMeanZFR; currentProbeData.nonCTX.allUnitsAvoidAtItiMeanZFR];
            probeData.nonCTX.allUnitsReactAtItiMeanZFR = [probeDataBefore.nonCTX.allUnitsReactAtItiMeanZFR; currentProbeData.nonCTX.allUnitsReactAtItiMeanZFR];
            probeData.nonCTX.allUnitsDiffMeanZFR = [probeDataBefore.nonCTX.allUnitsDiffMeanZFR; currentProbeData.nonCTX.allUnitsDiffMeanZFR]; % avoid - react
            probeData.nonCTX.allUnitsColdAirOnMeanZFR = [probeDataBefore.nonCTX.allUnitsColdAirOnMeanZFR; currentProbeData.nonCTX.allUnitsColdAirOnMeanZFR];
            probeData.nonCTX.allUnitsBeepOnAvoidMeanZFR = [probeDataBefore.nonCTX.allUnitsBeepOnAvoidMeanZFR; currentProbeData.nonCTX.allUnitsBeepOnAvoidMeanZFR];
            probeData.nonCTX.allUnitsBeepOnReactMeanZFR = [probeDataBefore.nonCTX.allUnitsBeepOnReactMeanZFR; currentProbeData.nonCTX.allUnitsBeepOnReactMeanZFR];
            probeData.nonCTX.allUnitsOptoOnMeanZFR = [probeDataBefore.nonCTX.allUnitsOptoOnMeanZFR; currentProbeData.nonCTX.allUnitsOptoOnMeanZFR];
            probeData.nonCTX.numGoodUnits = [probeDataBefore.nonCTX.numGoodUnits currentProbeData.nonCTX.numGoodUnits];
            probeData.nonCTX.stDepths = [probeDataBefore.nonCTX.stDepths; currentProbeData.nonCTX.stDepths];
            probeData.nonCTX.anatHiers = [probeDataBefore.nonCTX.anatHiers currentProbeData.nonCTX.anatHiers];

            probeData.nonCTX.platDistZAvoid = [probeDataBefore.nonCTX.platDistZAvoid {currentProbeData.nonCTX.platDistZAvoid}]; 
            probeData.nonCTX.platDistZReact = [probeDataBefore.nonCTX.platDistZReact {currentProbeData.nonCTX.platDistZReact}]; 
            probeData.nonCTX.platDistZOpto = [probeDataBefore.nonCTX.platDistZOpto {currentProbeData.nonCTX.platDistZOpto}]; 
            probeData.nonCTX.platDistZIti = [probeDataBefore.nonCTX.platDistZIti {currentProbeData.nonCTX.platDistZIti}]; 
            probeData.nonCTX.allUnitsTrialsAvoid = [probeDataBefore.nonCTX.allUnitsTrialsAvoid {currentProbeData.nonCTX.allUnitsTrialsAvoid}]; 
            probeData.nonCTX.allUnitsTrialsReact = [probeDataBefore.nonCTX.allUnitsTrialsReact {currentProbeData.nonCTX.allUnitsTrialsReact}]; 
            probeData.nonCTX.allUnitsTrialsOpto = [probeDataBefore.nonCTX.allUnitsTrialsOpto {currentProbeData.nonCTX.allUnitsTrialsOpto}]; 
            probeData.nonCTX.allUnitsTrialsIti = [probeDataBefore.nonCTX.allUnitsTrialsIti {currentProbeData.nonCTX.allUnitsTrialsIti}]; 
            probeData.nonCTX.allUnitsTrialsAvoidAtCue = [probeDataBefore.nonCTX.allUnitsTrialsAvoidAtCue {currentProbeData.nonCTX.allUnitsTrialsAvoidAtCue}]; 
            probeData.nonCTX.allUnitsTrialsReactAtCue = [probeDataBefore.nonCTX.allUnitsTrialsReactAtCue {currentProbeData.nonCTX.allUnitsTrialsReactAtCue}];  
            probeData.nonCTX.allUnitsTrialsOptoAtCue = [probeDataBefore.nonCTX.allUnitsTrialsOptoAtCue {currentProbeData.nonCTX.allUnitsTrialsOptoAtCue}]; 
            probeData.nonCTX.allUnitsTrialsColdAir = [probeDataBefore.nonCTX.allUnitsTrialsColdAir {currentProbeData.nonCTX.allUnitsTrialsColdAir}]; 
            probeData.nonCTX.allUnitsTrialsShuf1 = [probeDataBefore.nonCTX.allUnitsTrialsShuf1 {currentProbeData.nonCTX.allUnitsTrialsShuf1}];  
            probeData.nonCTX.allUnitsTrialsShuf2 = [probeDataBefore.nonCTX.allUnitsTrialsShuf2 {currentProbeData.nonCTX.allUnitsTrialsShuf2}];  
            probeData.nonCTX.allUnitsTrialsShufAtCue1 = [probeDataBefore.nonCTX.allUnitsTrialsShufAtCue1 {currentProbeData.nonCTX.allUnitsTrialsShufAtCue1}];  
            probeData.nonCTX.allUnitsTrialsShufAtCue2 = [probeDataBefore.nonCTX.allUnitsTrialsShufAtCue2 {currentProbeData.nonCTX.allUnitsTrialsShufAtCue2}]; 
        end
    end
end

function probeSubsetData = getProbeSubset(subsetIdx, shankData_test, tsData_test, idxData_test, times, analyzeWindow, runDecoders, minTrials, decodeWindow)
% get probe subset data defined by subsetIdx
    goodUnitsThresh = 10; %minimum units in subset to analyze; must be >3 to visualize 3 PCs, but can increase to match decoder threshold (i.e. 10)
    shankDataSubset.stCell = shankData_test.stCell(subsetIdx);

    % good units in subset
    if ~isempty(length(shankDataSubset.stCell))
        probeSubsetData.numGoodUnits = length(shankDataSubset.stCell);
        probeSubsetData.stDepths = shankData_test.stDepths(subsetIdx);
        probeSubsetData.anatHiers = shankData_test.anatHiers(subsetIdx);

        shankDataSubset.allLocalFR = shankData_test.allLocalFR(:,subsetIdx);
        shankDataSubset.allDeltaFR = shankData_test.allDeltaFR(:,subsetIdx);
        shankDataSubset.allDeltaFRzScore = shankData_test.allDeltaFRzScore(:,subsetIdx);
        shankDataSubset.psthTimes = shankData_test.psthTimes;
    else
        probeSubsetData.numGoodUnits = 0;
    end

    if  probeSubsetData.numGoodUnits > goodUnitsThresh 
        % actually collect data at event times here:
        [probeSubsetData.allUnitsMeanZFR, ~] = getEphysMeanZscoreFREvent(shankDataSubset, times.allVhexNoOpto, analyzeWindow);
        [probeSubsetData.allUnitsAvoidMeanZFR, ~] = getEphysMeanZscoreFREventEqualTrialsInclOpto(shankDataSubset, times.avoid, minTrials, analyzeWindow);
        [probeSubsetData.allUnitsReactMeanZFR, ~] = getEphysMeanZscoreFREventEqualTrialsInclOpto(shankDataSubset, times.react, minTrials, analyzeWindow);
        [probeSubsetData.allUnitsOptoMeanZFR, ~] = getEphysMeanZscoreFREventEqualTrialsInclOpto(shankDataSubset, times.optoVhex, minTrials, analyzeWindow);
   
        [probeSubsetData.allUnitsPlatformReleaseMeanZFR, ~] = getEphysMeanZscoreFREvent(shankDataSubset, times.platformRelease, analyzeWindow);
        [probeSubsetData.allUnitsPlatformCommandMeanZFR, ~] = getEphysMeanZscoreFREvent(shankDataSubset, times.platformCommand, analyzeWindow);
        [probeSubsetData.itiMeanZFR, ~] = getEphysMeanZscoreFREvent(shankDataSubset, times.iti, analyzeWindow);
        [probeSubsetData.allUnitsAvoidAtCueMeanZFR, ~] = getEphysMeanZscoreFREventEqualTrialsInclOpto(shankDataSubset, times.avoidAtCue, minTrials, analyzeWindow);
        [probeSubsetData.allUnitsReactAtCueMeanZFR, ~] = getEphysMeanZscoreFREventEqualTrialsInclOpto(shankDataSubset, times.reactAtCue, minTrials, analyzeWindow);
        [probeSubsetData.allUnitsOptoAtCueMeanZFR, ~] = getEphysMeanZscoreFREventEqualTrialsInclOpto(shankDataSubset, times.optoAtCue, minTrials, analyzeWindow);

        [probeSubsetData.allUnitsAvoidAtItiMeanZFR, ~] = getEphysMeanZscoreFREvent(shankDataSubset, times.avoidAtIti, analyzeWindow);
        [probeSubsetData.allUnitsReactAtItiMeanZFR, ~] = getEphysMeanZscoreFREvent(shankDataSubset, times.reactAtIti, analyzeWindow); 
        [probeSubsetData.allUnitsMeanZFR_shuf1, ~] = getEphysMeanZscoreFRShuffleEqualTrialsInclOpto(shankDataSubset, times.avoid, times.react, minTrials, analyzeWindow); %
        [probeSubsetData.allUnitsMeanZFR_shuf2, ~] = getEphysMeanZscoreFRShuffleEqualTrialsInclOpto(shankDataSubset, times.avoid, times.react, minTrials, analyzeWindow); %

        [probeSubsetData.allUnitsAtCueMeanZFR_shuf1, ~] = getEphysMeanZscoreFRShuffleEqualTrialsInclOpto(shankDataSubset, times.avoidAtCue, times.reactAtCue, minTrials, analyzeWindow); %
        [probeSubsetData.allUnitsAtCueMeanZFR_shuf2, ~] = getEphysMeanZscoreFRShuffleEqualTrialsInclOpto(shankDataSubset, times.avoidAtCue, times.reactAtCue, minTrials, analyzeWindow); %

        [probeSubsetData.allUnitsColdAirOnMeanZFR, ~] = getEphysMeanZscoreFREvent(shankDataSubset, times.coldAirOn, analyzeWindow);
        [probeSubsetData.allUnitsBeepOnAvoidMeanZFR, probeSubsetData.allUnitsBeepOnReactMeanZFR, ~] = getAvoidVsReactMeanZscoreFREvent(shankDataSubset, tsData_test.timescale, idxData_test, analyzeWindow);
        [probeSubsetData.allUnitsOptoOnMeanZFR, timescaleSeg] = getEphysMeanZscoreFREvent(shankDataSubset, times.optoOn, analyzeWindow);
        probeSubsetData.allUnitsDiffMeanZFR = []; %not used for now probeSubsetData.allUnitsAvoidMeanZFR - probeSubsetData.allUnitsReactMeanZFR;

        if runDecoders
            numDecoderSampRepeats = 10;
        
            for decRpt = 1:numDecoderSampRepeats
                disp(['sampling repeat ' num2str(decRpt) ' of ' num2str(numDecoderSampRepeats)])
                % note: allUnitsTrialsAvoid acts as movement null condition for @cue timing
                [platDistZAvoid{decRpt}, allUnitsTrialsAvoid{decRpt}, timescaleSegPlatDec] = getBehaviorEphysDecEventsBinned(shankDataSubset.allDeltaFRzScore, shankDataSubset.psthTimes, tsData_test.jumpDistance, tsData_test.timescale, times.avoid, decodeWindow, minTrials); %#ok<*AGROW> 
                [platDistZReact{decRpt}, allUnitsTrialsReact{decRpt}, ~] = getBehaviorEphysDecEventsBinned(shankDataSubset.allDeltaFRzScore, shankDataSubset.psthTimes, tsData_test.jumpDistance, tsData_test.timescale, times.react, decodeWindow, minTrials);
                [platDistZOpto{decRpt}, allUnitsTrialsOpto{decRpt}, ~] = getBehaviorEphysDecEventsBinned(shankDataSubset.allDeltaFRzScore, shankDataSubset.psthTimes, tsData_test.jumpDistance, tsData_test.timescale, times.optoVhex, decodeWindow, minTrials);
                [platDistZIti{decRpt}, allUnitsTrialsIti{decRpt}, ~] = getBehaviorEphysDecEventsBinned(shankDataSubset.allDeltaFRzScore, shankDataSubset.psthTimes, tsData_test.jumpDistance, tsData_test.timescale, times.iti, decodeWindow, minTrials);
                % also collect raw trial data at cue to facilitate downstream analyses: 
                [~, allUnitsTrialsAvoidAtCue{decRpt}, ~] = getBehaviorEphysDecEventsBinned(shankDataSubset.allDeltaFRzScore, shankDataSubset.psthTimes, tsData_test.jumpDistance, tsData_test.timescale, times.avoidAtCue, decodeWindow, minTrials);
                [~, allUnitsTrialsReactAtCue{decRpt}, ~] = getBehaviorEphysDecEventsBinned(shankDataSubset.allDeltaFRzScore, shankDataSubset.psthTimes, tsData_test.jumpDistance, tsData_test.timescale, times.reactAtCue, decodeWindow, minTrials);
                [~, allUnitsTrialsOptoAtCue{decRpt}, ~] = getBehaviorEphysDecEventsBinned(shankDataSubset.allDeltaFRzScore, shankDataSubset.psthTimes, tsData_test.jumpDistance, tsData_test.timescale, times.optoAtCue, decodeWindow, minTrials);
                % and save raw shuffle data for control comparisons:
                [allUnitsTrialsShuf1{decRpt}, ~] = getBehaviorEphysDecEventsBinnedShuffle(shankDataSubset.allDeltaFRzScore, shankDataSubset.psthTimes, times.avoid, times.react, decodeWindow, minTrials);
                [allUnitsTrialsShuf2{decRpt}, ~] = getBehaviorEphysDecEventsBinnedShuffle(shankDataSubset.allDeltaFRzScore, shankDataSubset.psthTimes, times.avoid, times.react, decodeWindow, minTrials);
                [allUnitsTrialsShufAtCue1{decRpt}, ~] = getBehaviorEphysDecEventsBinnedShuffle(shankDataSubset.allDeltaFRzScore, shankDataSubset.psthTimes, times.avoidAtCue, times.reactAtCue, decodeWindow, minTrials);
                [allUnitsTrialsShufAtCue2{decRpt}, ~] = getBehaviorEphysDecEventsBinnedShuffle(shankDataSubset.allDeltaFRzScore, shankDataSubset.psthTimes, times.avoidAtCue, times.reactAtCue, decodeWindow, minTrials);
            end

            probeSubsetData.platDistZAvoid = platDistZAvoid; 
            probeSubsetData.platDistZReact = platDistZReact; 
            probeSubsetData.platDistZOpto = platDistZOpto; 
            probeSubsetData.platDistZIti = platDistZIti; 
            probeSubsetData.allUnitsTrialsAvoid = allUnitsTrialsAvoid; 
            probeSubsetData.allUnitsTrialsReact = allUnitsTrialsReact; 
            probeSubsetData.allUnitsTrialsOpto = allUnitsTrialsOpto; 
            probeSubsetData.allUnitsTrialsIti = allUnitsTrialsIti;
            probeSubsetData.allUnitsTrialsAvoidAtCue = allUnitsTrialsAvoidAtCue; 
            probeSubsetData.allUnitsTrialsReactAtCue = allUnitsTrialsReactAtCue; 
            probeSubsetData.allUnitsTrialsOptoAtCue = allUnitsTrialsOptoAtCue; 
            probeSubsetData.allUnitsTrialsShuf1 = allUnitsTrialsShuf1;
            probeSubsetData.allUnitsTrialsShuf2 = allUnitsTrialsShuf2; 
            probeSubsetData.allUnitsTrialsShufAtCue1 = allUnitsTrialsShufAtCue1; 
            probeSubsetData.allUnitsTrialsShufAtCue2 = allUnitsTrialsShufAtCue2; 
            % cold air on is null condition for @ VHEx timing, so collect all trials here:
            [probeSubsetData.allUnitsTrialsColdAir, ~] = getBehaviorEphysDecEventsBinnedAll(shankDataSubset.allDeltaFRzScore, shankDataSubset.psthTimes, times.coldAirOn, decodeWindow);
            
        else
            probeSubsetData.platDistZAvoid = []; 
            probeSubsetData.platDistZReact = []; 
            probeSubsetData.platDistZOpto = []; 
            probeSubsetData.platDistZIti = []; 
            probeSubsetData.allUnitsTrialsAvoid = []; 
            probeSubsetData.allUnitsTrialsReact = []; 
            probeSubsetData.allUnitsTrialsOpto = []; 
            probeSubsetData.allUnitsTrialsIti = []; 
            probeSubsetData.allUnitsTrialsAvoidAtCue = []; 
            probeSubsetData.allUnitsTrialsReactAtCue = []; 
            probeSubsetData.allUnitsTrialsOptoAtCue = []; 
            probeSubsetData.allUnitsTrialsColdAir = [];
            probeSubsetData.allUnitsTrialsShuf1 = []; 
            probeSubsetData.allUnitsTrialsShuf2 = []; 
            probeSubsetData.allUnitsTrialsShufAtCue1 = []; 
            probeSubsetData.allUnitsTrialsShufAtCue2 = [];
        end

    else
        % assign default empty arrays in case not enough units found on this probe
        probeSubsetData.allUnitsMeanZFR = []; % for mean across all VHEx
        probeSubsetData.allUnitsMeanZFR_shuf1 = []; % make 2 shuffled versions with randomly selected avoid vs react trials, for control comparions
        probeSubsetData.allUnitsMeanZFR_shuf2 = [];
        probeSubsetData.allUnitsAtCueMeanZFR_shuf1 = []; % make 2 shuffled versions with randomly selected avoid vs react trials, for control comparions
        probeSubsetData.allUnitsAtCueMeanZFR_shuf2 = [];
        probeSubsetData.allUnitsAvoidMeanZFR = [];
        probeSubsetData.allUnitsReactMeanZFR = [];
        probeSubsetData.allUnitsOptoMeanZFR = [];
        probeSubsetData.allUnitsPlatformReleaseMeanZFR = [];
        probeSubsetData.allUnitsPlatformCommandMeanZFR = [];
        probeSubsetData.itiMeanZFR = [];
        probeSubsetData.allUnitsAvoidAtCueMeanZFR = [];
        probeSubsetData.allUnitsReactAtCueMeanZFR = [];
        probeSubsetData.allUnitsAvoidAtItiMeanZFR = [];
        probeSubsetData.allUnitsReactAtItiMeanZFR = [];
        probeSubsetData.allUnitsOptoAtCueMeanZFR = [];
        probeSubsetData.allUnitsDiffMeanZFR = []; % avoid - react
        probeSubsetData.allUnitsColdAirOnMeanZFR = [];
        probeSubsetData.allUnitsBeepOnAvoidMeanZFR = [];
        probeSubsetData.allUnitsBeepOnReactMeanZFR = [];
        probeSubsetData.allUnitsOptoOnMeanZFR = [];
        probeSubsetData.numGoodUnits = [];
        probeSubsetData.stDepths = [];
        probeSubsetData.anatHiers = [];
        probeSubsetData.platDistZAvoid = []; 
        probeSubsetData.platDistZReact = []; 
        probeSubsetData.platDistZOpto = []; 
        probeSubsetData.platDistZIti = []; 
        probeSubsetData.allUnitsTrialsAvoid = []; 
        probeSubsetData.allUnitsTrialsReact = []; 
        probeSubsetData.allUnitsTrialsOpto = []; 
        probeSubsetData.allUnitsTrialsIti = []; 
        probeSubsetData.allUnitsTrialsAvoidAtCue = []; 
        probeSubsetData.allUnitsTrialsReactAtCue = []; 
        probeSubsetData.allUnitsTrialsOptoAtCue = []; 
        probeSubsetData.allUnitsTrialsColdAir = [];
        probeSubsetData.allUnitsTrialsShuf1 = []; 
        probeSubsetData.allUnitsTrialsShuf2 = []; 
        probeSubsetData.allUnitsTrialsShufAtCue1 = []; 
        probeSubsetData.allUnitsTrialsShufAtCue2 = [];
    end %if numGoodUnits > thresh
end

function allMiceAllUnitMeans = calculateMeanAcrossTrialsRepeats(anatData, numMiceProbes, removeNull)
% calculate means across trials and repeats with standard concatenation across time & conditions for subspace calculations
% return matrix rows as units, each column having concatenated timing across 4 conditions, with mean responses across all trials and repeats
    numRepeats = 10; %should match number in other EqualTrials functions
    numTimepoints = size(anatData.allUnitsTrialsAvoid{1}{1}{1},2);
    numConditions = 4;

    allMiceAllUnitMeans = [];
    for i = 1:numMiceProbes 
        % concatenate over 4 conditions that we care about variance across: avoid/react trial types & cue/VHEx timing
        if removeNull
            thisMouseData1 = anatData.allUnitsTrialsAvoidAtCueNull{i}; %unpack to get per-mouse cells that include all sample repeats
            thisMouseData2 = anatData.allUnitsTrialsAvoidNull{i};
            thisMouseData3 = anatData.allUnitsTrialsReactAtCueNull{i};
            thisMouseData4 = anatData.allUnitsTrialsReactNull{i};
        else
            thisMouseData1 = anatData.allUnitsTrialsAvoidAtCue{i}; %unpack to get per-mouse cells that include all sample repeats
            thisMouseData2 = anatData.allUnitsTrialsAvoid{i};
            thisMouseData3 = anatData.allUnitsTrialsReactAtCue{i};
            thisMouseData4 = anatData.allUnitsTrialsReact{i};
        end
        numTrials = size(thisMouseData1{1},2);
        numUnits = size(thisMouseData1{1}{1},1);
        thisMouseConcatRptsTrials = NaN(numUnits, numTrials*numRepeats, numTimepoints*numConditions);
        endTrialRpt = 1;
        for rpt = 1:numRepeats
            thisRptData1 = thisMouseData1{rpt}; %unpack one more time to get each repeat cell - now this cell for each trial, each having units x timepoints matrix
            thisRptData2 = thisMouseData2{rpt};
            thisRptData3 = thisMouseData3{rpt};
            thisRptData4 = thisMouseData4{rpt};
            for j = 1:numTrials
                for k = 1:numUnits
                    thisMouseConcatRptsTrials(k,endTrialRpt,:) = [thisRptData1{j}(k,:) thisRptData2{j}(k,:) thisRptData3{j}(k,:) thisRptData4{j}(k,:)]; %avoid/react @ cue/VHEx concatenated
                end
                endTrialRpt = endTrialRpt + 1;
            end %end trials
        end %end sampling repeats
        thisMouseMeanRptsTrials = squeeze(mean(thisMouseConcatRptsTrials,2)); %take mean across trials and repeats
        allMiceAllUnitMeans = [allMiceAllUnitMeans; thisMouseMeanRptsTrials];
    end %end numMice
end

function dataNullRemoved = getNullRemovedData(anatData, timescaleSeg)
% function to take trial data corresponding to a particular region and remove information (e.g. cold air) using nullspace projection
% calculate null for each mouse & probe separately

    condNamesVhex = {'allUnitsTrialsAvoid','allUnitsTrialsReact','allUnitsTrialsOpto','allUnitsTrialsShuf1','allUnitsTrialsShuf2'};
    condNamesCue = {'allUnitsTrialsAvoidAtCue','allUnitsTrialsReactAtCue','allUnitsTrialsOptoAtCue','allUnitsTrialsShufAtCue1','allUnitsTrialsShufAtCue2'};
    
    removeNull = true;
    varCutoff = 0.999; %how much null to remove, in terms of varExplained

    nullCondVhex = 'allUnitsTrialsColdAir'; %remove cold air info @ VHEx
    nullCondCue = 'allUnitsTrialsAvoid'; %remove avoid movement info @ cue
%     nullCondCue = 'allUnitsTrialsReact'; %option to remove react movement info @ cue %for Fig. 4f removing react movement from opto trials (see below)
%     nullCondCue = 'allUnitsTrialsIti'; %control
    nullWindowVhex = [-0.5 0.2];  %for cold air on
%     nullWindowVhex = [-0.5 -0.2];  %control
    nullWindowCue = [-1 0]; %for movement
    nullSampsVhex = [find(nullWindowVhex(1)>=timescaleSeg,1,'last') find(nullWindowVhex(2)>=timescaleSeg,1,'last')];
    nullSampsCue = [find(nullWindowCue(1)>=timescaleSeg,1,'last') find(nullWindowCue(2)>=timescaleSeg,1,'last')];
    
    numMiceProbes = length(anatData.numGoodUnits);
    numRepeats = 10; %should match number in other EqualTrials functions
    smoothingSamples = [5 0]; %ex. smooth window 5 is 5*20ms PSTH bins = 100ms; for causal half-kernel use [X 0]
    
    % first get means across trials & repeats for null conditions, and calculate mean across trials (per repeat) for other conditions to facilitate PC distance calculations later
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for cond = condNamesVhex
        condNameMean = [cond{1} 'Mean'];
%         dataNullRemoved.(cond{1}) = anatData.(cond{1}); %also re-save original data in new structure
        for i = 1:numMiceProbes %loop over mice/probes now (those with enough units) to apply PC projection and concatenate output
            thisMouseCondData = anatData.(cond{1}){i}; %unpack to get per-mouse cells that include all sample repeats
            numTrials = size(thisMouseCondData{1},2);
            numUnits = size(thisMouseCondData{1}{1},1);
            numTimepoints = size(thisMouseCondData{1}{1},2);
            for rpt = 1:numRepeats
                thisRptCondData = thisMouseCondData{rpt}; %unpack one more time to get each repeat cell - now this cell for each trial, each having units x timepoints matrix
                dataForMean = zeros(numTrials, numUnits, numTimepoints);
                for j = 1:numTrials
                        dataForMean(j,:,:) = thisRptCondData{j};
                        trialSmoothed = smoothdata(thisRptCondData{j}, 2, 'gaussian', smoothingSamples); %add some smoothing to try and match PCA-denoised subspace decoder
                        dataNullRemoved.(cond{1}){i}{rpt}{j} = trialSmoothed; %also re-save original data in new structure (with option to smooth)
                end %end trials
                dataNullRemoved.(condNameMean){i}{rpt} = squeeze(mean(dataForMean,1)); % save mean across trials for each unit
            end %end sampling repeats
            
            if strcmp(cond,condNamesVhex{1}) %since number of trials is different for cold air null + no repeats, compute in separate loop, only once per mouse/probe
                thisMouseNullCondVhex = anatData.(nullCondVhex){i}; %per-mouse cell containing all trial cells 
                numNullTrials = size(thisMouseNullCondVhex,2);
                numNullTimepoints = nullSampsVhex(2) - nullSampsVhex(1) + 1;
                thisMouseNullVhexConcatTrials = NaN(numNullTrials, numUnits, numNullTimepoints);
                for jn = 1:numNullTrials
                    thisNullTrial = thisMouseNullCondVhex{jn}(:,nullSampsVhex(1):nullSampsVhex(2));
                    thisMouseNullVhexConcatTrials(jn,:,:) = thisNullTrial;
                end
                thisMouseMeanNullVhexTrials = squeeze(mean(thisMouseNullVhexConcatTrials,1)); %take mean across trials
                allMiceAllTrialsNullDataVhex{i} = thisMouseMeanNullVhexTrials;
            end
        end %end numMice
    end %end condition
    
    for cond = condNamesCue
        condNameMean = [cond{1} 'Mean'];
        for i = 1:numMiceProbes %loop over mice/probes now (those with enough units) to apply PC projection and concatenate output
            thisMouseCondData = anatData.(cond{1}){i}; %unpack to get per-mouse cells that include all sample repeats
            numTrials = size(thisMouseCondData{1},2);
            numUnits = size(thisMouseCondData{1}{1},1);
            numTimepoints = size(thisMouseCondData{1}{1},2);

            if strcmp(cond,condNamesCue{1}) %only compute null once, but here null has repeats so need to take mean across trials and repeats
                thisMouseNullCondCue = anatData.(nullCondCue){i}; %per-mouse cell 
                numTrials = size(thisMouseCondData{1},2);
                numTimepointsNull = nullSampsCue(2) - nullSampsCue(1) + 1;
                thisMouseNullCueConcatRptsTrials = NaN(numTrials*numRepeats, numUnits, numTimepointsNull);
                endTrialRpt = 1;
            end

            for rpt = 1:numRepeats
                thisRptCondData = thisMouseCondData{rpt}; %unpack one more time to get each repeat cell - now this cell for each trial, each having units x timepoints matrix
                thisRptNullData = thisMouseNullCondCue{rpt};
                dataForMean = zeros(numTrials, numUnits, numTimepoints);
                for j = 1:numTrials
                    if strcmp(cond,condNamesCue{1})
                        thisNullTrial = thisRptNullData{j}(:,nullSampsCue(1):nullSampsCue(2));
                        thisMouseNullCueConcatRptsTrials(endTrialRpt,:,:) = thisNullTrial;
                        endTrialRpt = endTrialRpt + 1;
                    end
                    dataForMean(j,:,:) = thisRptCondData{j};
                    trialSmoothed = smoothdata(thisRptCondData{j}, 2, 'gaussian', smoothingSamples); %add some smoothing to try and match PCA-denoised subspace decoder
                    dataNullRemoved.(cond{1}){i}{rpt}{j} = trialSmoothed; %also re-save original data in new structure (with option to smooth)
                end %end trials
                dataNullRemoved.(condNameMean){i}{rpt} = squeeze(mean(dataForMean,1)); % save mean across trials for each unit
            end %end sampling repeats

            if strcmp(cond,condNamesCue{1})  %only compute null once
                thisMouseMeanNullCueRptsTrials = squeeze(mean(thisMouseNullCueConcatRptsTrials,1)); %take mean across trials and repeats
                allMiceAllTrialsNullDataCue{i} = thisMouseMeanNullCueRptsTrials;
            end
        end %end numMice
    end %end condition

    % now calculate nulls separately for each mouse
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:numMiceProbes

        thisMouseMeanNullDataVhex = allMiceAllTrialsNullDataVhex{i}; 
        thisMouseMeanNullDataVhex = thisMouseMeanNullDataVhex - mean(thisMouseMeanNullDataVhex,2);
        thisMouseMeanNullDataVhexNorm = thisMouseMeanNullDataVhex./norm(thisMouseMeanNullDataVhex,'fro');
    
        thisMouseMeanNullDataCue = allMiceAllTrialsNullDataCue{i};
        thisMouseMeanNullDataCue = thisMouseMeanNullDataCue - mean(thisMouseMeanNullDataCue,2);
        thisMouseMeanNullDataCueNorm = thisMouseMeanNullDataCue./norm(thisMouseMeanNullDataCue,'fro');
    
        % since the mean null condition matrices are generally full rank, just set a varExplained threshold for null dimensions:
        [U,S,V_nullCue] = svd(thisMouseMeanNullDataCueNorm');
        varExplainedNullCue = cumsum(diag(S).^2)./sum(diag(S).^2);
        varExplCueIdx = find(varExplainedNullCue>varCutoff,1,'first'); % find cutoff where variance explained reaches threshold
        %figure; plot(varExplainedNullCue); xlabel('dimensions'); ylabel('varExplained'); title('cue null'); makepretty;
        nullCue = V_nullCue(:,varExplCueIdx+1:end);
        nullCueProj{i} = nullCue*nullCue'; %make projection matrix to get back to original activity space

        [U,S,V_nullVhex] = svd(thisMouseMeanNullDataVhexNorm');
        varExplainedNullVhex = cumsum(diag(S).^2)./sum(diag(S).^2);
        varExplVhexIdx = find(varExplainedNullVhex>varCutoff,1,'first');
        %figure; plot(varExplainedNullCue); xlabel('dimensions'); ylabel('varExplained'); title('VHEx null'); makepretty;
        nullVhex = V_nullVhex(:,varExplVhexIdx+1:end);
        nullVhexProj{i} = nullVhex*nullVhex'; 

%         nullVhexProj{i} = nullCueProj{i}; %optionally use the VHEx null as control; ex. for Fig. 4f removing avoid/react movement activity on opto trials
    end

    %finally, go back through trial data and project into the nullspace for each condition:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for cond = condNamesVhex
        nullCond = [cond{1} 'Null'];
        nullCondMean = [cond{1} 'NullMean'];
        for i = 1:numMiceProbes %loop over mice/probes now (those with enough units) to apply PC projection and concatenate output
            thisMouseCondData = anatData.(cond{1}){i}; %unpack to get per-mouse cells that include all sample repeats
            for rpt = 1:numRepeats
                thisRptCondData = thisMouseCondData{rpt}; %unpack one more time to get each repeat cell - now this cell for each trial, each having units x timepoints matrix
                numTrials = size(thisRptCondData,2);
                dataForNullMean = [];
                for j = 1:numTrials
                    trialProjCond = nullVhexProj{i}*thisRptCondData{j};
                    trialProjCond = smoothdata(trialProjCond, 2, 'gaussian', smoothingSamples); %add some smoothing to try and match PCA-denoised subspace decoder
                    dataNullRemoved.(nullCond){i}{rpt}{j} = trialProjCond;
                    dataForNullMean(j,:,:) = trialProjCond;
%                         figure; imagesc(thisRptCondData{j}); clim([-1 1]);
%                         figure; imagesc(trialProjCond); clim([-1 1]);
                end %end trials
                meanData = squeeze(mean(dataForNullMean,1)); %mean across trials
                dataNullRemoved.(nullCondMean){i}{rpt} = meanData;
            end %end sampling repeats
        end %end numMice
    end %end condition
    
    for cond = condNamesCue
        nullCond = [cond{1} 'Null'];
        nullCondMean = [cond{1} 'NullMean'];
        for i = 1:numMiceProbes %loop over mice/probes now (those with enough units) to apply PC projection and concatenate output
            thisMouseCondData = anatData.(cond{1}){i}; %unpack to get per-mouse cells that include all sample repeats
            for rpt = 1:numRepeats
                thisRptCondData = thisMouseCondData{rpt}; %unpack one more time to get each repeat cell - now this cell for each trial, each having units x timepoints matrix
                numTrials = size(thisRptCondData,2);
                dataForNullMean = [];
                for j = 1:numTrials
%                         trialProjCond = nullCue'*thisRptCondData{j}; %apply PCs to individual trials; output is #PCsToUse x #timepoints
                    trialProjCond = nullCueProj{i}*thisRptCondData{j};
                    trialProjCond = smoothdata(trialProjCond, 2, 'gaussian', smoothingSamples);
                    dataNullRemoved.(nullCond){i}{rpt}{j} = trialProjCond;
                    dataForNullMean(j,:,:) = trialProjCond;
%                         figure; imagesc(thisRptCondData{j}); clim([-1 1]);
%                         figure; imagesc(trialProjCond); clim([-1 1]);
                end %end trials
                meanData = squeeze(mean(dataForNullMean,1)); %mean across trials
                dataNullRemoved.(nullCondMean){i}{rpt} = meanData;
            end %end sampling repeats
%                 temp = squeeze(mean(dataForNullMean,1)); figure; imagesc(temp(1,:)); clim([-1 1]);
        end %end numMice
    end %end condition

    % for each mouse/probe and condition, optionally compare original mean activity and null-removed mean activity (for 1 trial sampling repeat)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     climAll = [-2 2];
%     for cond = condNamesVhex
%         condMean = [cond{1} 'Mean'];
%         nullCondMean = [cond{1} 'NullMean'];
%         f = figure; t = tiledlayout(2,numMiceProbes,'TileSpacing', 'tight');
%         for i = 1:numMiceProbes
%             meanData = dataNullRemoved.(condMean){i}{rpt};
%             meanDataNull = dataNullRemoved.(nullCondMean){i}{rpt};
%             tileNum = tilenum(t,1,i);
%             nexttile(tileNum); imagesc(meanData); clim(climAll); if i==1; title(condMean); end
%             tileNum = tilenum(t,2,i);
%             nexttile(tileNum); imagesc(meanDataNull); clim(climAll); if i==1; title(nullCondMean); end
%         end
%         f.Position(3:4) = [1800 600]; %set figure size appropriate for plots
%     end
% 
%     for cond = condNamesCue
%         condMean = [cond{1} 'Mean'];
%         nullCondMean = [cond{1} 'NullMean'];
%         f = figure; t = tiledlayout(2,numMiceProbes,'TileSpacing', 'tight');
%         for i = 1:numMiceProbes
%             meanData = dataNullRemoved.(condMean){i}{rpt};
%             meanDataNull = dataNullRemoved.(nullCondMean){i}{rpt};
%             tileNum = tilenum(t,1,i);
%             nexttile(tileNum); imagesc(meanData); clim(climAll); if i==1; title(condMean); end
%             tileNum = tilenum(t,2,i);
%             nexttile(tileNum); imagesc(meanDataNull); clim(climAll); if i==1; title(nullCondMean); end
%         end
%         f.Position(3:4) = [1800 600]; %set figure size appropriate for plots
%     end

end

function subspaceData = getSubspace(anatData, timescaleSeg)
% function to take trial data corresponding to a particular region and put it into a common subspace using PCA,
% optionally removing information (e.g. cold air) using nullspace projection

    condNamesVhex = {'allUnitsTrialsAvoid','allUnitsTrialsReact','allUnitsTrialsOpto','allUnitsTrialsShuf1','allUnitsTrialsShuf2'};
    condNamesCue = {'allUnitsTrialsAvoidAtCue','allUnitsTrialsReactAtCue','allUnitsTrialsOptoAtCue','allUnitsTrialsShufAtCue1','allUnitsTrialsShufAtCue2'};
    
    removeNull = true;

    alignPCs = true;
    nullCondVhex = 'allUnitsTrialsColdAir'; %remove cold air info @ VHEx
    nullCondCue = 'allUnitsTrialsAvoid'; %remove movement info @ cue
%     nullCondCue = 'allUnitsTrialsIti'; %control
    nullWindowVhex = [-0.5 0.2];  %for cold air on
%     nullWindowVhex = [-0.5 -0.2];  %control
    nullWindowCue = [-1 0]; %for movement
    nullSampsVhex = [find(nullWindowVhex(1)>=timescaleSeg,1,'last') find(nullWindowVhex(2)>=timescaleSeg,1,'last')];
    nullSampsCue = [find(nullWindowCue(1)>=timescaleSeg,1,'last') find(nullWindowCue(2)>=timescaleSeg,1,'last')];
    
    pcaWindowVhex = [-2 2]; %seconds around extension to calculate PCs for subspace
    pcaWindowCue = [-2 2]; %seconds around CUE to calculate PCs for subspace; include some ITI time to allow "readiness" as part of dimensionality
    pcaSampsVhex = [find(pcaWindowVhex(1)>=timescaleSeg,1,'last') find(pcaWindowVhex(2)>=timescaleSeg,1,'last')];
    pcaSampsCue = [find(pcaWindowCue(1)>=timescaleSeg,1,'last') find(pcaWindowCue(2)>=timescaleSeg,1,'last')];
    
    numMiceProbes = length(anatData.numGoodUnits);
    numPcsToUse = 10;
    smoothingSamples = [25 0]; %ex. smooth window 5 is 5*20ms PSTH bins = 100ms; for causal half-kernel use [X 0]
    numRepeats = 10; %should match number in other EqualTrials functions
    
    % first GET PCs (based on mean activity across all conditions, aligned at both cue & VHEx & concatenated)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    numTimepoints = (pcaSampsCue(2)-pcaSampsCue(1)+1)*2 + (pcaSampsVhex(2)-pcaSampsVhex(1)+1)*2;
    numConditions = 4;
    
    % NOTE: this loop is same as calculateMeanAcrossTrialsRepeats() but using pcsSamps to restrict timing
    allMiceAllUnitMeans = [];
    for i = 1:numMiceProbes 
        % concatenate over 4 conditions that we care about variance across: avoid/react trial types & cue/VHEx timing
        thisMouseData1 = anatData.allUnitsTrialsAvoidAtCue{i}; %unpack to get per-mouse cells that include all sample repeats
        thisMouseData2 = anatData.allUnitsTrialsAvoid{i};
        thisMouseData3 = anatData.allUnitsTrialsReactAtCue{i};
        thisMouseData4 = anatData.allUnitsTrialsReact{i};
        numTrials = size(thisMouseData1{1},2);
        numUnits = size(thisMouseData1{1}{1},1);
        thisMouseConcatRptsTrials = NaN(numUnits, numTrials*numRepeats, numTimepoints);
        endTrialRpt = 1;
        for rpt = 1:numRepeats
            thisRptData1 = thisMouseData1{rpt}; %unpack one more time to get each repeat cell - now this cell for each trial, each having units x timepoints matrix
            thisRptData2 = thisMouseData2{rpt};
            thisRptData3 = thisMouseData3{rpt};
            thisRptData4 = thisMouseData4{rpt};
            for j = 1:numTrials
                for k = 1:numUnits
                    thisMouseConcatRptsTrials(k,endTrialRpt,:) = [thisRptData1{j}(k,pcaSampsCue(1):pcaSampsCue(2)) thisRptData2{j}(k,pcaSampsVhex(1):pcaSampsVhex(2)) thisRptData3{j}(k,pcaSampsCue(1):pcaSampsCue(2)) thisRptData4{j}(k,pcaSampsVhex(1):pcaSampsVhex(2))]; %avoid/react @ cue/VHEx concatenated
                end
                endTrialRpt = endTrialRpt + 1;
            end %end trials
        end %end sampling repeats
        thisMouseMeanRptsTrials = squeeze(mean(thisMouseConcatRptsTrials,2)); %take mean across trials and repeats
        allMiceAllUnitMeans = [allMiceAllUnitMeans; thisMouseMeanRptsTrials];
    end %end numMice

    allMiceAllUnitMeansFullTime = calculateMeanAcrossTrialsRepeats(anatData, numMiceProbes, false); % also calculate means across entire timespan for comparison with below
%     figure; imagesc(allMiceAllUnitMeansFullTime); clim([-0.5 0.5]);

    % mean subtract and normalize again so that some mice don't dominate PCA
    allMiceAllUnitMeans = allMiceAllUnitMeans - mean(allMiceAllUnitMeans,2);
    allMiceAllUnitMeansNorm = allMiceAllUnitMeans./norm(allMiceAllUnitMeans,'fro');
%     figure; imagesc(allMiceAllUnitMeansNorm); clim([-0.005 0.005]);

    [pcSeg, varExplained] = getEphysPCsFromData(allMiceAllUnitMeansNorm, anatData.numGoodUnits, numMiceProbes, numPcsToUse, alignPCs); %

    % now project trial data for each condition onto common aligned & re-orthogonalized PCs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for cond = condNamesVhex
        condNameMean = [cond{1} 'Mean'];
        for i = 1:numMiceProbes %loop over mice/probes now (those with enough units) to apply PC projection and concatenate output
            thisMouseCondData = anatData.(cond{1}){i}; %unpack to get per-mouse cells that include all sample repeats
            if strcmp(cond,condNamesVhex{1}) %only compute null once
                thisMouseNullCondVhex = anatData.(nullCondVhex){i}; %per-mouse cell containing all trial cells  
            end
            for rpt = 1:numRepeats
                thisRptCondData = thisMouseCondData{rpt}; %unpack one more time to get each repeat cell - now this cell for each trial, each having units x timepoints matrix
                numTrials = size(thisRptCondData,2);
                dataForMean = [];
                for j = 1:numTrials
                    trialProjCond = pcSeg{i}'*thisRptCondData{j}; %apply PCs to individual trials; output is #PCsToUse x #timepoints
                    trialProjCond = smoothdata(trialProjCond, 2, 'gaussian', smoothingSamples); %ex. smooth window 25 is 25*20ms PSTH bins = 500ms
                    subspaceData.(cond{1}){i}{rpt}{j} = trialProjCond;
                    dataForMean(j,:,:) = trialProjCond; %collect data for mean across trials as well, to make  PC distance calculation easier later
                end %end trials
                subspaceData.(condNameMean){i}{rpt} = squeeze(mean(dataForMean,1));
            end %end sampling repeats
            if strcmp(cond,condNamesVhex{1})
                numNullTrials = size(thisMouseNullCondVhex,2);
                numTimepoints = nullSampsVhex(2) - nullSampsVhex(1) + 1;
                thisMouseNullVhexConcatRptsTrials = NaN(numNullTrials, numPcsToUse, numTimepoints);
                for jn = 1:numNullTrials
                    thisNullTrial = thisMouseNullCondVhex{jn}(:,nullSampsVhex(1):nullSampsVhex(2));
                    nullTrialProjCond = pcSeg{i}'*thisNullTrial; %put null data into supbspace as well to prepare for PCA below
                    subspaceData.(nullCondVhex){i}{jn} = nullTrialProjCond;
                    thisMouseNullVhexConcatRptsTrials(jn,:,:) = nullTrialProjCond;
                end
                thisMouseMeanNullVhexRptsTrials = squeeze(mean(thisMouseNullVhexConcatRptsTrials,1)); %take mean across trials
                allMiceAllTrialsNullDataVhex(i,:,:) = thisMouseMeanNullVhexRptsTrials;
            end
        end %end numMice
    end %end condition
    
    for cond = condNamesCue
        condNameMean = [cond{1} 'Mean'];
        for i = 1:numMiceProbes %loop over mice/probes now (those with enough units) to apply PC projection and concatenate output
            thisMouseCondData = anatData.(cond{1}){i}; %unpack to get per-mouse cells that include all sample repeats
            if strcmp(cond,condNamesCue{1}) %only compute null once
                thisMouseNullCondCue = anatData.(nullCondCue){i}; %per-mouse cell 
                numTrials = size(thisMouseCondData{1},2);
                numTimepoints = nullSampsCue(2) - nullSampsCue(1) + 1;
                thisMouseNullCueConcatRptsTrials = NaN(numTrials*numRepeats, numPcsToUse, numTimepoints);
                endTrialRpt = 1;
            end
            for rpt = 1:numRepeats
                thisRptCondData = thisMouseCondData{rpt}; %unpack one more time to get each repeat cell - now this cell for each trial, each having units x timepoints matrix
                thisRptNullData = thisMouseNullCondCue{rpt};
                numTrials = size(thisRptCondData,2);
                dataForMean = [];
%                 rawDataForMean = [];
                for j = 1:numTrials
                    trialProjCond = pcSeg{i}'*thisRptCondData{j}; %apply PCs to individual trials; output is #PCsToUse x #timepoints
                    trialProjCond = smoothdata(trialProjCond, 2, 'gaussian', smoothingSamples); %ex. smooth window 25 is 25*20ms PSTH bins = 500ms
                    subspaceData.(cond{1}){i}{rpt}{j} = trialProjCond;
                    dataForMean(j,:,:) = trialProjCond;
%                     rawDataForMean(j,:,:) = thisRptCondData{j};
                    if strcmp(cond,condNamesCue{1})
                        thisNullTrial = thisRptNullData{j}(:,nullSampsCue(1):nullSampsCue(2));
                        nullTrialProjCond = pcSeg{i}'*thisNullTrial; %put null data into supbspace as well to prepare for PCA below
                        thisMouseNullCueConcatRptsTrials(endTrialRpt,:,:) = nullTrialProjCond;
                        endTrialRpt = endTrialRpt + 1;
                    end
                end %end trials
                meanData = squeeze(mean(dataForMean,1)); %mean across trials
%                 meanRawData = squeeze(mean(rawDataForMean,1)); %mean across trials
%                 figure; imagesc(meanData); caxis([-1 1]);
%                 figure; imagesc(meanRawData); caxis([-1 1]);
                subspaceData.(condNameMean){i}{rpt} = meanData; %mean across repeats
            end %end sampling repeats
            if strcmp(cond,condNamesCue{1})
                thisMouseMeanNullCueRptsTrials = squeeze(mean(thisMouseNullCueConcatRptsTrials,1)); %take mean across trials and repeats
                allMiceAllTrialsNullDataCue(i,:,:) = thisMouseMeanNullCueRptsTrials;
            end
        end %end numMice
    end %end condition

    subspaceData.numGoodUnits = anatData.numGoodUnits;
    allMiceAllUnitMeansProjected = calculateMeanAcrossTrialsRepeats(subspaceData, numMiceProbes, false);
%     figure; tiledlayout(3,3,'TileSpacing', 'tight');
%     for dim=1:9
%         nexttile;
%         imagesc(allMiceAllUnitMeansProjected(dim:10:end,:)); clim([-2 2]); %for plotting mean projection of first dimension, across mice...use (2:10:end) for 2nd dimension, etc
%     end

    % finally, calculate null based on subspace dimensions and remove, if necessary
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if removeNull
        allMiceMeanNullDataVhex = squeeze(mean(allMiceAllTrialsNullDataVhex,1)); %mean across mice
        allMiceMeanNullDataVhex = allMiceMeanNullDataVhex - mean(allMiceMeanNullDataVhex,2);
        allMiceMeanNullDataVhexNorm = allMiceMeanNullDataVhex./norm(allMiceMeanNullDataVhex,'fro');

        allMiceMeanNullDataCue = squeeze(mean(allMiceAllTrialsNullDataCue,1)); %mean across mice
        allMiceMeanNullDataCue = allMiceMeanNullDataCue - mean(allMiceMeanNullDataCue,2);
        allMiceMeanNullDataCueNorm = allMiceMeanNullDataCue./norm(allMiceMeanNullDataCue,'fro');

        % calculate null space as the null of the top X PCs for each condition; more PCs removes more info
        [U,S,V_nullCue] = svd(allMiceMeanNullDataCueNorm');
        varExplainedNullCue = cumsum(diag(S).^2)./sum(diag(S).^2);
%         figure; plot(varExplainedNullCue); xlabel('dimensions'); ylabel('varExplained'); title('cue null'); makepretty;
        nullCue = V_nullCue(:,end); %take last 1-2 dimension with least information remaining, to construct projection matrix
        nullCueProj = nullCue*nullCue';
        [U,S,V_nullVhex] = svd(allMiceMeanNullDataVhexNorm');
        varExplainedNullVhex = cumsum(diag(S).^2)./sum(diag(S).^2);
%         figure; plot(varExplainedNullCue); xlabel('dimensions'); ylabel('varExplained'); title('VHEx null'); makepretty;
        nullVhex = V_nullVhex(:,end);
        nullVhexProj = nullVhex*nullVhex';
%         nullVhexProj = nullCueProj; % control

%     nullForData = null(covForNull,0.0001); %use reasonable tolerance of 0.0001 (otherwise can end up with empty null matrix for cov matrix that is barely full rank); this will give a nullspace [# neurons x # null dimensions], where # null dimensions is less than the number of cells
    % ex.
%     [U, S, V] = svd(covForNull);
%     tol = 0.0001; %max(size(covForNull)) * eps(norm(S, 'fro'));
%     r = sum(diag(S) > tol); %rank
%     nullForData = V(:, r+1:end);

        %finally, go back through trial data and project into the nullspace:
        for cond = condNamesVhex
            nullCond = [cond{1} 'Null'];
            nullCondMean = [cond{1} 'NullMean'];
            for i = 1:numMiceProbes %loop over mice/probes now (those with enough units) to apply PC projection and concatenate output
                thisMouseCondData = subspaceData.(cond{1}){i}; %unpack to get per-mouse cells that include all sample repeats
                for rpt = 1:numRepeats
                    thisRptCondData = thisMouseCondData{rpt}; %unpack one more time to get each repeat cell - now this cell for each trial, each having units x timepoints matrix
                    numTrials = size(thisRptCondData,2);
                    dataForNullMean = [];
                    for j = 1:numTrials
%                         trialProjCond = nullVhex'*thisRptCondData{j}; %apply PCs to individual trials; output is #PCsToUse x #timepoints
                        trialProjCond = nullVhexProj*thisRptCondData{j};
%                         trialProjCond = smoothdata(trialProjCond, 2, 'gaussian', smoothingSamples); %ex. smooth window 25 is 25*20ms PSTH bins = 500ms
                        subspaceData.(nullCond){i}{rpt}{j} = trialProjCond;
                        dataForNullMean(j,:,:) = trialProjCond;
%                         figure; imagesc(thisRptCondData{j}); clim([-1 1]);
%                         figure; imagesc(trialProjCond); clim([-1 1]);
                    end %end trials
                    meanData = squeeze(mean(dataForNullMean,1)); %mean across trials
                    subspaceData.(nullCondMean){i}{rpt} = meanData;
                end %end sampling repeats
            end %end numMice
        end %end condition
        
        for cond = condNamesCue
            nullCond = [cond{1} 'Null'];
            nullCondMean = [cond{1} 'NullMean'];
            for i = 1:numMiceProbes %loop over mice/probes now (those with enough units) to apply PC projection and concatenate output
                thisMouseCondData = subspaceData.(cond{1}){i}; %unpack to get per-mouse cells that include all sample repeats
                for rpt = 1:numRepeats
                    thisRptCondData = thisMouseCondData{rpt}; %unpack one more time to get each repeat cell - now this cell for each trial, each having units x timepoints matrix
                    numTrials = size(thisRptCondData,2);
                    dataForNullMean = [];
                    for j = 1:numTrials
%                         trialProjCond = nullCue'*thisRptCondData{j}; %apply PCs to individual trials; output is #PCsToUse x #timepoints
                        trialProjCond = nullCueProj*thisRptCondData{j};
%                         trialProjCond = smoothdata(trialProjCond, 2, 'gaussian', smoothingSamples); %ex. smooth window 25 is 25*20ms PSTH bins = 500ms
                        subspaceData.(nullCond){i}{rpt}{j} = trialProjCond;
                        dataForNullMean(j,:,:) = trialProjCond;
%                         figure; imagesc(thisRptCondData{j}); clim([-1 1]);
%                         figure; imagesc(trialProjCond); clim([-1 1]);
                    end %end trials
                    meanData = squeeze(mean(dataForNullMean,1)); %mean across trials
                    subspaceData.(nullCondMean){i}{rpt} = meanData;
                end %end sampling repeats
%                 temp = squeeze(mean(dataForNullMean,1)); figure; imagesc(temp(1,:)); clim([-1 1]);
            end %end numMice
        end %end condition
    end %end if removeNull

    allMiceAllUnitMeansNullProjected = calculateMeanAcrossTrialsRepeats(subspaceData, numMiceProbes, true);
%     figure; tiledlayout(3,3,'TileSpacing', 'tight');
%     for dim=1:9
%         nexttile;
%         imagesc(allMiceAllUnitMeansNullProjected(dim:10:end,:)); clim([-2 2]); %for plotting mean projection of first dimension, across mice...use (2:10:end) for 2nd dimension, etc
%     end

end

function [behTrialsZ, ephysTrials, timescaleSeg] = getBehaviorEphysDecEventsBinned(subsetZFR, subsetTimescale, behaviorData, behaviorTimescale, eventTimes, eventWindow, minTrials)
% function to collect z-scored behavior measure in matrix of windows trials, but binned to match decoder predictor timescale (here, ephys bins)
% also collect cell array of [units/predictors x ZFR] for each trial (use cell since # predictors can change for each mouse)
    if isempty(eventTimes)
        behTrialsZ = [];
        ephysTrials = [];
        timescaleSeg = [];
        return
    end
    
    % randomly choose minTrials to include
    numEvents = length(eventTimes);
    s = RandStream('dsfmt19937','Seed','shuffle');
    randIdx = randperm(s, numEvents, minTrials); %draw randomly # minimum of avoid/react/opto trials for this mouse
    eventTimesShuf = eventTimes(randIdx);
    numUnits = size(subsetZFR,2);
    smoothWindow = [5 0]; % how many PSTH bins in window for Gaussian smoothing below; ex. 20ms bins x 10 bins = 200ms; causal half kernel use [X 0]; keep on the order of neurotransmitter time constants

    % first normalize behavior data to z-scored version
    behZ = normalize(behaviorData);

    sampTimeEphys = (subsetTimescale(2) - subsetTimescale(1));
    timescaleSeg = eventWindow(1):sampTimeEphys:eventWindow(2)-sampTimeEphys;
    samplesBefore = round(eventWindow(1)/sampTimeEphys);
    samplesAfter = round(eventWindow(2)/sampTimeEphys);
 
    % now collect across events 
    ephysTrials = [];
    behTrialsZ = zeros(length(eventTimesShuf),length(timescaleSeg));
    thisEphysTrial = zeros(numUnits,length(timescaleSeg)); %#ok<PREALL> 
    for j = 1:length(eventTimesShuf)
        currentEphysIdx = find(subsetTimescale>eventTimesShuf(j), 1, 'first');
        thisEphysTrial = subsetZFR(currentEphysIdx+samplesBefore:currentEphysIdx+samplesAfter-1,:);
        thisEphysTrial = smoothdata(thisEphysTrial,1,'gaussian',smoothWindow); %smooth Z-score over time to reduce noise for decoders
        thisSegTimes = subsetTimescale(currentEphysIdx+samplesBefore:currentEphysIdx+samplesAfter-1);
        ephysTrials = [ephysTrials {thisEphysTrial'}];
        % for platform distance behavior, downsample in time to match ephys PSTH bin size; since niSampRate and imSampRate are not same 
        % and not synced, go through SLOW process of finding indeces where they match best; first get approximate reduced timescale to search
        % to reduce the search time considerably:
        firstSearchSample =  find(behaviorTimescale>thisSegTimes(1), 1, 'first') - 10;
        lastSearchSample =  find(behaviorTimescale>thisSegTimes(end), 1, 'first') + 10;
        behaviorTimescaleToSearch = behaviorTimescale(firstSearchSample:lastSearchSample);
        for t = 1:length(thisSegTimes)
            currentBehTimeIdx = find(behaviorTimescaleToSearch>thisSegTimes(t), 1, 'first'); %this is index into reduced timescale, so add back on firstSearchSample
            behTrialsZ(j,t) = behZ(currentBehTimeIdx + firstSearchSample - 1);
        end
    end
end

function [ephysTrials, timescaleSeg] = getBehaviorEphysDecEventsBinnedShuffle(subsetZFR, subsetTimescale, eventTimes1, eventTimes2, eventWindow, minTrials)
% function to collect shuffled z-scored matrix of windowed trials

    %use min of the two # of events to pick from each distribution, to match normal EqualTrials sampling procedure
    numEvents1 = length(eventTimes1);
    numEvents2 = length(eventTimes2);
    numShufEvents = min([numEvents1 numEvents2 minTrials]);

    if mod(numShufEvents,2)==0 %if even
        numEventsToChoose1 = numShufEvents/2;
        numEventsToChoose2 = numShufEvents/2;
    else %if odd
        numEventsToChoose1 = (numShufEvents+1)/2;
        numEventsToChoose2 = (numShufEvents-1)/2;
    end

    s = RandStream('dsfmt19937','Seed','shuffle');
    randIdx1 = randperm(s, numEvents1, numEventsToChoose1); %take equal # random trials from each type of event
    randIdx2 = randperm(s, numEvents2, numEventsToChoose2);
    eventTimesShuf = [eventTimes1(randIdx1) eventTimes2(randIdx2)]; %1/2 times from event1, 1/2 times from event2

    numUnits = size(subsetZFR,2);
    smoothWindow = [5 0]; % how many PSTH bins in window for Gaussian smoothing below; ex. 20ms bins x 10 bins = 200ms; causal half kernel use [X 0]; keep on the order of neurotransmitter time constants

    sampTimeEphys = (subsetTimescale(2) - subsetTimescale(1));
    timescaleSeg = eventWindow(1):sampTimeEphys:eventWindow(2)-sampTimeEphys;
    samplesBefore = round(eventWindow(1)/sampTimeEphys);
    samplesAfter = round(eventWindow(2)/sampTimeEphys);
 
    % now collect across events 
    ephysTrials = [];
    thisEphysTrial = zeros(numUnits,length(timescaleSeg)); %#ok<PREALL> 
    for j = 1:length(eventTimesShuf)
        currentEphysIdx = find(subsetTimescale>eventTimesShuf(j), 1, 'first');
        thisEphysTrial = subsetZFR(currentEphysIdx+samplesBefore:currentEphysIdx+samplesAfter-1,:);
        thisEphysTrial = smoothdata(thisEphysTrial,1,'gaussian',smoothWindow); %smooth Z-score over time to reduce noise for decoders
        ephysTrials = [ephysTrials {thisEphysTrial'}];
    end

end

function [ephysTrials, timescaleSeg] = getBehaviorEphysDecEventsBinnedAll(subsetZFR, subsetTimescale, eventTimes, eventWindow)
% function to collect z-scored behavior measure in matrix of windows trials, but binned to match decoder predictor timescale (here, ephys bins)
% also collect cell array of [units/predictors x ZFR] for each trial (use cell since # predictors can change for each mouse)
% this version does NOT select randomly from min # avoid/react/opto trials, but takes all trials (to be used for null condition computation)
    if isempty(eventTimes)
        ephysTrials = [];
        timescaleSeg = [];
        return
    end
    
    numUnits = size(subsetZFR,2);
    smoothWindow = [5 0]; % how many PSTH bins in window for Gaussian smoothing below; ex. 20ms bins x 10 bins = 200ms; causal half kernel use [X 0]; keep on the order of neurotransmitter time constants

    sampTimeEphys = (subsetTimescale(2) - subsetTimescale(1));
    timescaleSeg = eventWindow(1):sampTimeEphys:eventWindow(2)-sampTimeEphys;
    samplesBefore = round(eventWindow(1)/sampTimeEphys);
    samplesAfter = round(eventWindow(2)/sampTimeEphys);
 
    % now collect across events 
    ephysTrials = [];
    thisEphysTrial = zeros(numUnits,length(timescaleSeg)); %#ok<PREALL> 
    for j = 1:length(eventTimes)
        currentEphysIdx = find(subsetTimescale>eventTimes(j), 1, 'first');
        thisEphysTrial = subsetZFR(currentEphysIdx+samplesBefore:currentEphysIdx+samplesAfter-1,:);
        thisEphysTrial = smoothdata(thisEphysTrial,1,'gaussian',smoothWindow); %smooth Z-score over time to reduce noise for decoders
        ephysTrials = [ephysTrials {thisEphysTrial'}];
    end
end

function [sortIdx] = sortEphysDeltaAroundEvent(shankDataToSort, eventTimes, window)
    %function to sort ephys data by magnitude of deltaFR (from average, signed) after some event(s)
    % window is [-secondsBeforeEvent secondsAfterEvent]
    numUnits = size(shankDataToSort.allDeltaFR,2);
    eventWindowBefore = [window(1) window(2)]; %windows around events to calculate deltaFR for sorting
    eventWindowAfter = [window(3) window(4)]; %
    psthBinSize = shankDataToSort.psthTimes(2) - shankDataToSort.psthTimes(1);
    
    allBeforeMinusAfterEventFR = zeros(1, numUnits);
    k = 1;
    
    % NOTE: should find latency to peak modulation and multiply mean by that to get sort order (or use covariance diagonal or more complicated RasterMap algorithm)
    while(k<numUnits+1) %go through units to get mean before and after around event times
        % use zscore so easier to compare across units
        [meanFrBefore, ~, ~] = meanSemFiringRateAcrossEvents(shankDataToSort.allDeltaFRzScore(:,k), eventTimes, psthBinSize, eventWindowBefore);
        [meanFrAfter, ~, ~] = meanSemFiringRateAcrossEvents(shankDataToSort.allDeltaFRzScore(:,k), eventTimes, psthBinSize, eventWindowAfter);
        allBeforeMinusAfterEventFR(k) = mean(meanFrBefore) - mean(meanFrAfter); %positive # for decrease after event, negative for increase after event
        k = k + 1;
    end
    
    [~,sortIdx] = sort(allBeforeMinusAfterEventFR, 'descend'); %clusters with largest decrease after event first
end

function [sortIdx] = sortEphysDeltaAroundEventNegOnly(shankDataToSort, eventTimes, window)
    %function to sort ephys data by magnitude of deltaFR (from average, signed) after some event(s)
    % ONLY KEEP UNITS WITH NEGATIVE MODULATION (ex. throw out putative inhibtory neurons for vgat-Chr2 effect quantification)
    % window is [-secondsBeforeEvent secondsAfterEvent]
    numUnits = size(shankDataToSort.allDeltaFR,2);
    eventWindowBefore = [window(1) window(2)]; %windows around events to calculate deltaFR for sorting
    eventWindowAfter = [window(3) window(4)]; %
    psthBinSize = shankDataToSort.psthTimes(2) - shankDataToSort.psthTimes(1);
    
    allBeforeMinusAfterEventFR = zeros(1, numUnits);
    k = 1;
    
    % NOTE: should find latency to peak modulation and multiply mean by that to get sort order (or use covariance diagonal or more complicated RasterMap algorithm)
    while(k<numUnits+1) %go through units to get mean before and after around event times
        % use zscore so easier to compare across units
        [meanFrBefore, ~, ~] = meanSemFiringRateAcrossEvents(shankDataToSort.allLocalFR(:,k), eventTimes, psthBinSize, eventWindowBefore);
        [meanFrAfter, ~, ~] = meanSemFiringRateAcrossEvents(shankDataToSort.allLocalFR(:,k), eventTimes, psthBinSize, eventWindowAfter);
        allBeforeMinusAfterEventFR(k) = mean(meanFrBefore) - mean(meanFrAfter); %positive # for decrease after event, negative for increase after event
        k = k + 1;
    end
    
    [~,sortIdxAll] = sort(allBeforeMinusAfterEventFR, 'descend'); %clusters with largest decrease after event first

    %go through units one more time and throw out those indeces with positive modulated units
    logicalKeepIdx = true(1, length(sortIdxAll));
    for k = 1:length(sortIdxAll)
        currIdx = sortIdxAll(k);
        if allBeforeMinusAfterEventFR(currIdx) < 0
            logicalKeepIdx(k) = false;
        end
    end
    sortIdx = sortIdxAll(logicalKeepIdx);
end

function [sortIdx] = sortEphysZMeansAroundEvent(allUnitsZFR, timescaleSeg, windows)
    %function to sort ephys data which is already organized by mean
    %z-scored FR around an event, respose at window(3:4) minus response at window(1:2) 
    numUnits = size(allUnitsZFR,1);
    resortSampsBefore = [find(windows(1)>=timescaleSeg,1,'last') find(windows(2)>=timescaleSeg,1,'last')];
    resortSampsAfter = [find(windows(3)>=timescaleSeg,1,'last') find(windows(4)>=timescaleSeg,1,'last')];
    
    allBeforeMinusAfterEventFR = zeros(1, numUnits);
    % NOTE: can also find latency to peak modulation and multiply mean by that to get sort order (or use covariance diagonal or more complicated RasterMap algorithm)
    for i = 1:numUnits
        meanZFrBefore = mean(allUnitsZFR(i,resortSampsBefore(1):resortSampsBefore(2)));
        meanZFrAfter = mean(allUnitsZFR(i,resortSampsAfter(1):resortSampsAfter(2)));
        allBeforeMinusAfterEventFR(i) = meanZFrBefore - meanZFrAfter; %positive # for decrease after event, negative for increase after event
    end

    [~,sortIdx] = sort(allBeforeMinusAfterEventFR, 'descend'); %clusters with largest decrease after event first
end

function [sortIdx] = rastermapSort(zScoredFrData)
% wrapper to convert to numpy array, call Rastermap, and then convert back to MATLAB array
% this assumes Python interpreter already setup with something like:
% pyenv('Version','C:\Users\labadmin\.conda\envs\rastermap\pythonw.exe', 'ExecutionMode', 'OutOfProcess')
% also set "PYTHONHOME" environmental variable to to rastermap environment root folder first
% see also:
% https://www.mathworks.com/help/matlab/matlab_external/install-supported-python-implementation.html

    % first check if python ev already loaded, and load if not
    pe = pyenv; 
    if pe.Status ~= "Loaded" 
        pyenv('Version','C:\Users\labadmin\.conda\envs\rastermap\pythonw.exe', 'ExecutionMode', 'OutOfProcess')
    end

    data = zScoredFrData';
    dataNdArray = py.numpy.array(data);
    pyrun("from rastermap import Rastermap") %load interpreter, import main function
    rmModel = pyrun("model = Rastermap(locality=0.0, time_lag_window=100).fit(spks)", "model", spks=dataNdArray);
    sortIdx = int16(py.memoryview(rmModel.isort.data)) + 1; %back to MATLAB array, 1-indexing

end

function plotEphysClassicRaster(tsData, idxData, stCell, sortIdx) 
    % rasterize from cell array of spike times and plot; pass sortIdx = 0 if not sorting (use default depth sorting)
    if sortIdx ~= 0 
        rasterIdx = sortIdx;
    else
        rasterIdx = 1:1:length(stCell);
    end
    currentIdx = 1;
    figure;
    hold on;
    %colored patches when opto on
    y = [1 length(stCell) length(stCell) 1];
    for j = 1:length(idxData.optoOnIndeces)
        x = [tsData.timescale(idxData.optoOnIndeces(j)) tsData.timescale(idxData.optoOnIndeces(j)) tsData.timescale(idxData.optoOffIndeces(j)) tsData.timescale(idxData.optoOffIndeces(j))];
        patch(x, y, [0 0.5 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); %blue box for stimulus
    end
    %colored patches when cold air on
    for j = 1:length(idxData.solenoidOnIndeces)
        x = [tsData.timescale(idxData.solenoidOnIndeces(j)) tsData.timescale(idxData.solenoidOnIndeces(j)) tsData.timescale(idxData.solenoidOffIndeces(j)) tsData.timescale(idxData.solenoidOffIndeces(j))];
        patch(x, y, [0.7 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); %blue box for stimulus
    end
    %movement overlay scaled across Y-axis
    jumpDistanceScaled = (tsData.jumpDistance./max(tsData.jumpDistance)).*length(stCell);
    plot(tsData.timescale, jumpDistanceScaled, 'Color',[0 1 0 0.6], 'LineWidth', 1) %4 element color vector for transparency
%     %opto voltage overlay scaled across Y-axis
%     optoVoltageScaled = (tsData.optoVoltage./max(tsData.optoVoltage)).*length(stCell);
%     plot(tsData.timescale, optoVoltageScaled, 'Color',[0 0 1 0.8], 'LineWidth', 1) %4 element color vector for transparency

    for k = rasterIdx %now go through each row,unit, rasterize & plot
        rasterX = nan(1,length(stCell{k})*3); %use 3rd NaN column so plot is discontinuous
        rasterX(1:3:end) = stCell{k}; %every first & second point is the spike time, but third point is NaN
        rasterX(2:3:end) = stCell{k};
        rasterY = nan(1,length(stCell{k})*3);
        minValues = ones(1,length(stCell{k})).*currentIdx; %for each unit, the raster ticks are 1 a.u.
        rasterY(1:3:end) = minValues;
        rasterY(2:3:end) = minValues+1;
        plot(rasterX,rasterY,'Color',[0 0 0 0.2], 'Marker', 'none', 'LineWidth', 0.01); %4 element color vector for transparency
        currentIdx = currentIdx + 1;
    end
    ylabel('cluster/unit #');  
    xlabel('time (sec)')

    hold off;
    axis tight;
    makepretty;
    set(gcf,'color','w'); %set figure background white
end

function plotEphysClassicRasterByAnatomy(tsData, idxData, stCell, anatHiers, movementIndex, frameTimes) 
    % rasterize from cell array of spike times and plot; pass sortIdx = 0 if not sorting (use default depth sorting)
    rasterIdx = 1:1:length(stCell);
%     rasterIdx = 120:330; %option to plot a subset of cells   
    platformToBeepTime = 0.5;
    cueTimes = tsData.timescale(idxData.beepOnIndeces) - platformToBeepTime; %just use platfrom command time here

    figure;
    hold on;

    %restrict to trial times (defined as 0.25sec before platform command [cue] to 3.0sec after VHEx thresh) to make comparisons easier
    trialIdxDiff = diff(tsData.trialTimes); % 1 when trial on, 0 otherwise, so get edges using diff
    trialStartIdx = find(trialIdxDiff==1);
    trialEndIdx = find(trialIdxDiff==-1) - 22500; %10k samp rate and 3 sec added before, so subtract to make end at extension + 750ms
    trialStartTimes = tsData.timescale(trialStartIdx); %250ms before platform command
    trialEndTimes = tsData.timescale(trialEndIdx);
    timeStep = tsData.timescale(2) - tsData.timescale(1);
    [unitCmaps, unitCmapsIdx, anatCmapsPlot, anatNamesPlot] = anatHiersToCmaps(anatHiers); %get anatomy colormaps
%     jumpDistanceScaled = (tsData.jumpDistance./max(tsData.jumpDistance)).*length(rasterIdx);
    jumpDistanceScaled = tsData.jumpDistance.*10; %non-scaled option for examples with actual distance
    movementIndexScaled = movementIndex.*length(rasterIdx);
    currentTrialPlotStartTime = 0;

    optoIdx = find(idxData.jumpOptoIndeces);
    avoidIdx = find(idxData.jumpAntIndeces);
    nonOptoIdx = find(~idxData.jumpOptoIndeces);
    allIdxToPlot = sort([avoidIdx optoIdx],'ascend');

%     for t = 1:length(trialStartTimes) % use multiple or custom list here to restrict plot to certain trials
    for t = nonOptoIdx(2:6)
%     for t = [13:15,25:29,46:50]
        trialTimesteps = trialEndIdx(t) - trialStartIdx(t) + 1;
        currentTrialPlotEndTime = currentTrialPlotStartTime + trialTimesteps*timeStep;
        trialTimescale = currentTrialPlotStartTime:timeStep:(currentTrialPlotEndTime-timeStep);

        currentIdx = 1; 

        for k = rasterIdx %now go through each row,unit, rasterize & plot; for every trial need to re-rasterize :-(, so takes a minute or two

            rasterX = stCell{k}; % all spike times for this unit
            firstSpike = find(rasterX>trialStartTimes(t),1,'first'); % index of first spike time in trial
            lastSpike = find(rasterX<trialEndTimes(t),1,'last');
            rasterXSeg = rasterX(firstSpike:lastSpike) - trialStartTimes(t) + currentTrialPlotStartTime; %make relative to trial start time, then add current concatenation time
    
            rasterXslash = nan(1,length(rasterXSeg)*3); %use 3rd NaN column so plot is discontinuous
            rasterXslash(1:3:end) = rasterXSeg; %every first & second point is the spike time, but third point is NaN
            rasterXslash(2:3:end) = rasterXSeg;
            rasterY = nan(1,length(rasterXSeg)*3);
            minValues = ones(1,length(rasterXSeg)).*currentIdx; %for each unit, the raster ticks are 1 a.u.
            rasterY(1:3:end) = minValues;
            rasterY(2:3:end) = minValues+1;
            plot(rasterXslash,rasterY,'Color',unitCmaps(k,:), 'Marker', 'none', 'LineWidth', 2); %increase LineWidth for zoom capability
            currentIdx = currentIdx + 1;
        end
        
        % thick black lines for trial boundaries
        y = [1 length(rasterIdx)];
        x = [currentTrialPlotStartTime currentTrialPlotStartTime];
        plot(x, y, 'k', 'LineWidth',3);

        % green dotted line for cue timing (which is constant offset for each trial)
        x = [currentTrialPlotStartTime+0.25 currentTrialPlotStartTime+0.25];
        plot(x, y, 'g--', 'LineWidth',1);

        % magenta dotted line for extension threshold
        lastVhexBeforeTrialEnd = find(idxData.jumpIndeces<trialEndIdx(t),1,'last'); %index of last VHEx before the trial ended
        lastVhexTime = tsData.timescale(idxData.jumpIndeces(lastVhexBeforeTrialEnd));
        lastVhexTime = lastVhexTime - trialStartTimes(t) + currentTrialPlotStartTime; %make relative to trial start time, then add current concatenation time
        x = [lastVhexTime lastVhexTime];
        plot(x, y, 'm--', 'LineWidth',1);
        
        % blue colored patches when opto on
        lastOptoOnBeforeTrialEnd = find(idxData.optoOnIndeces<trialEndIdx(t),1,'last'); %index of last opto on before the trial ended
        lastOptoOnTime = tsData.timescale(idxData.optoOnIndeces(lastOptoOnBeforeTrialEnd));
        if (lastOptoOnTime > trialStartTimes(t)) %test if opto trial
            lastOptoOffTime = tsData.timescale(idxData.optoOffIndeces(lastOptoOnBeforeTrialEnd));
            lastOptoOnTime = lastOptoOnTime - trialStartTimes(t) + currentTrialPlotStartTime; %make relative to trial start time, then add current concatenation time
            lastOptoOffTime = lastOptoOffTime - trialStartTimes(t) + currentTrialPlotStartTime; 
            y = [1 length(rasterIdx) length(rasterIdx) 1];
            x = [lastOptoOnTime lastOptoOnTime lastOptoOffTime lastOptoOffTime];
            patch(x, y, [0 0.5 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); %
        end

        % red colored patches when cold air on
        lastColdOnBeforeTrialEnd = find(idxData.solenoidOnIndeces<trialEndIdx(t),1,'last'); %index of last cold on before the trial ended
        lastColdOnTime = tsData.timescale(idxData.solenoidOnIndeces(lastColdOnBeforeTrialEnd));
        if (lastColdOnTime > trialStartTimes(t)) %test if cold / react trial
            lastColdOffTime = tsData.timescale(idxData.solenoidOffIndeces(lastColdOnBeforeTrialEnd));
            lastColdOnTime = lastColdOnTime - trialStartTimes(t) + currentTrialPlotStartTime; %make relative to trial start time, then add current concatenation time
            lastColdOffTime = lastColdOffTime - trialStartTimes(t) + currentTrialPlotStartTime; 
            y = [1 length(rasterIdx) length(rasterIdx) 1];
            x = [lastColdOnTime lastColdOnTime lastColdOffTime lastColdOffTime];
            patch(x, y, [0.7 0 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
        end

        % green platform movement overlay scaled across Y-axis
        jumDistanceSeg = jumpDistanceScaled(trialStartIdx(t):trialEndIdx(t));
        plot(trialTimescale, jumDistanceSeg, 'Color',[0 0.7 0 0.8], 'LineWidth', 2) %4 element color vector for transparency

        % blue movement index overlay scaled across Y-axis
        firstFrame = find(frameTimes>trialStartTimes(t),1,'first'); % index of first frame of trial
        lastFrame = find(frameTimes<trialEndTimes(t),1,'last');
        trialTimescaleFrames = frameTimes(firstFrame):(1/400):frameTimes(lastFrame);
        trialTimescaleFramesAdj = trialTimescaleFrames - trialStartTimes(t) + currentTrialPlotStartTime;
        movementIndexSeg = movementIndexScaled(firstFrame:lastFrame);
        plot(trialTimescaleFramesAdj, movementIndexSeg, 'Color',[0.1 0.1 0.7 0.8], 'LineWidth', 2) %4 element color vector for transparency

        currentTrialPlotStartTime = currentTrialPlotEndTime + timeStep;

    end %end trials

    ylabel('cluster/unit #');  
    xlabel('time (sec)')
    drawAnatColormapLegend(anatCmapsPlot, anatNamesPlot)
    hold off;
    axis tight;
    makepretty;
    set(gcf,'color','w'); %set figure background white
end

function plotEphysTrialRaster(tsData, idxData, numData, stCell, anatHiers) 
    % plot a single neuron raster over trials, with various annotations
    % rasterize from cell array of spike times and plot; pass sortIdx = 0 if not sorting (use default depth sorting)
    [rasterIdx, ~] = getAnatSubset(anatHiers, 'Isocortex');
%     rasterIdx = 301:360; % indeces of cells to plot (single plots in a loop)

    platformToBeepTime = 0.5;
    cueTimes = tsData.timescale(idxData.beepOnIndeces) - platformToBeepTime; % use platfrom command time here for cue
    blueCmap = [0 0.4470 0.7410];
    orangeCmap = [0.8500 0.3250 0.0980];

    for k = rasterIdx %now go through each row,unit, rasterize & plot; for every trial need to re-rasterize :-(, so takes a minute or two

        figure;
        hold on;
    
        %restrict to trial times (defined as 0.25sec before platform command [cue] to 3.0sec after VHEx thresh) to make comparisons easier
        trialIdxDiff = diff(tsData.trialTimes); % 1 when trial on, 0 otherwise, so get edges using diff
        trialStartIdx = find(trialIdxDiff==1);
        trialStartTimes = tsData.timescale(trialStartIdx)-0.25; %250ms before platform command, so add another 250ms to make it 500ms before cue
        trialEndTimes = trialStartTimes + 6.5; % here trial end times are a constant 6 sec after cue
        timeStep = tsData.timescale(2) - tsData.timescale(1); % ephys timescale
        trialTimescale = -0.5:timeStep:6;
        trialTimesteps = length(trialTimescale);

        optoIdx = find(idxData.jumpOptoIndeces);
        avoidIdx = find(idxData.jumpAntIndeces);
        nonOptoIdx = find(~idxData.jumpOptoIndeces(1:end-1)); %don't use last since it can differ from trialStart index
        allIdxToPlot = sort([avoidIdx optoIdx],'ascend');
        allLatencies = numData.jumpLatencies + 3; % add back in 3 seconds since latency relative to cold air (beepToColdAir from readBehaviorSaveMat.m)
        nonOptoLatencies = allLatencies(nonOptoIdx);

        % sort trials by latency, and filter for < 6sec for plotting 
        [~,sortIdx] = sort(nonOptoLatencies, 'ascend');
        nonOptoIdx = nonOptoIdx(sortIdx);
        nonOptoLatencies = nonOptoLatencies(sortIdx);
        lastIdx = find(nonOptoLatencies>6, 1, 'first');
        if ~isempty(lastIdx)
            nonOptoIdx = nonOptoIdx(1:lastIdx-1);
            nonOptoLatencies = nonOptoLatencies(1:lastIdx-1);
        end
        
        currentTrialIdx = 1;
%     for t = 1:length(trialStartTimes) % use multiple or custom list here to restrict plot to certain trials
        for tr = nonOptoIdx % all non-opto trials
    
            rasterX = stCell{k}; % all spike times for this unit
            firstSpike = find(rasterX>trialStartTimes(tr),1,'first'); % index of first spike time in trial
            lastSpike = find(rasterX<trialEndTimes(tr),1,'last');
            rasterXSeg = rasterX(firstSpike:lastSpike) - trialStartTimes(tr) - 0.5; %make relative to cue at t=0: trial start time - 0.5
    
            rasterXslash = nan(1,length(rasterXSeg)*3); %use 3rd NaN column so plot is discontinuous
            rasterXslash(1:3:end) = rasterXSeg; %every first & second point is the spike time, but third point is NaN
            rasterXslash(2:3:end) = rasterXSeg;
            rasterY = nan(1,length(rasterXSeg)*3);
            minValues = ones(1,length(rasterXSeg)).*currentTrialIdx; %for each unit, the raster ticks are 1 a.u.
            rasterY(1:3:end) = minValues;
            rasterY(2:3:end) = minValues+1;
            plot(rasterXslash,rasterY,'Color','k', 'Marker', 'none', 'LineWidth', 1); %increase LineWidth for zoom capability

            % orange (avoid) or blue (react) dot for extension threshold timing for each trial
%             lastVhexBeforeTrialEnd = find(idxData.jumpIndeces<trialEndIdx(t),1,'last'); %index of last VHEx before the trial ended
%             lastVhexTime = tsData.timescale(idxData.jumpIndeces(lastVhexBeforeTrialEnd));
%             lastVhexTime = lastVhexTime - trialStartTimes(t); %make relative to trial start time, then add current concatenation time
            lastVhexTime = allLatencies(tr); % after correction above, this is relative to time = 0 at platform command cue
            if lastVhexTime > 3 %react if latency more than 3 seconds
                trialVhexCmap = blueCmap;
            else
                trialVhexCmap = orangeCmap;
            end
            plot(lastVhexTime, currentTrialIdx+0.5, 'Marker', "diamond", 'MarkerEdgeColor', trialVhexCmap, 'MarkerFaceColor', trialVhexCmap, 'MarkerSize',5);

            currentTrialIdx = currentTrialIdx + 1;
        end %end trials

    % green dotted line for cue timing
    x = [0 0];
    y = [1 length(nonOptoIdx)];
    plot(x, y, 'g--', 'LineWidth',3);

%     % cyan dotted line when cold air on
%     x = [3 3];
%     plot(x, y, 'c--', 'LineWidth',1);

    set(gca, 'YDir','reverse')
    yticks([])
    xlim([-0.5 6.1])
    xticks([0 1 2 3 4 5 6])
    ylabel('trials');  
    xlabel('time (sec)')
    hold off;
    axis tight;
    makepretty;
    set(gcf,'color','w'); %set figure background white

    end %end k units

end

function plotMeanOverUnits(allUnitsCtxAvoidMeanZFR_Probe0, allUnitsCtxReactMeanZFR_Probe0, allUnitsCtxAvoidMeanZFR_Probe1, allUnitsCtxReactMeanZFR_Probe1, allUnitsCtxAvoidMeanZFR_Probe2, allUnitsCtxReactMeanZFR_Probe2, timescaleSeg, blueCmap, orangeCmap)
 %combined ctx or non-ctx units avoid vs react mean responses
    [meanLatencyAvoidCtx_Probe0, semLatencyAvoidCtx_Probe0] = grpstats(allUnitsCtxAvoidMeanZFR_Probe0,[],{'mean' 'sem'}); %mean across columns/timepoints
    [meanLatencyAvoidCtx_Probe1, semLatencyAvoidCtx_Probe1] = grpstats(allUnitsCtxAvoidMeanZFR_Probe1,[],{'mean' 'sem'}); 
    [meanLatencyAvoidCtx_Probe2, semLatencyAvoidCtx_Probe2] = grpstats(allUnitsCtxAvoidMeanZFR_Probe2,[],{'mean' 'sem'}); 
    [meanLatencyReactCtx_Probe0, semLatencyReactCtx_Probe0] = grpstats(allUnitsCtxReactMeanZFR_Probe0,[],{'mean' 'sem'}); %mean across columns/timepoints
    [meanLatencyReactCtx_Probe1, semLatencyReactCtx_Probe1] = grpstats(allUnitsCtxReactMeanZFR_Probe1,[],{'mean' 'sem'}); 
    [meanLatencyReactCtx_Probe2, semLatencyReactCtx_Probe2] = grpstats(allUnitsCtxReactMeanZFR_Probe2,[],{'mean' 'sem'}); 

    figure;
    tiledlayout(1,3,'TileSpacing', 'tight');

    nexttile; hold on;
    boundedline(timescaleSeg, meanLatencyReactCtx_Probe0, semLatencyReactCtx_Probe0, 'cmap', blueCmap)
    boundedline(timescaleSeg, meanLatencyAvoidCtx_Probe0, semLatencyAvoidCtx_Probe0, 'cmap', orangeCmap)
    zeroXDottedLine;
    hold off;
    ylabel('mean z-scored FR');
    xlabel('time (sec)');
    title('Probe 0 mean')

    nexttile; hold on;
    boundedline(timescaleSeg, meanLatencyReactCtx_Probe1, semLatencyReactCtx_Probe1, 'cmap', blueCmap)
    boundedline(timescaleSeg, meanLatencyAvoidCtx_Probe1, semLatencyAvoidCtx_Probe1, 'cmap', orangeCmap)
    zeroXDottedLine;
    hold off;
    xlabel('time (sec)');
    title('Probe 1 mean')

    nexttile; hold on;
    [hl1,~] = boundedline(timescaleSeg, meanLatencyReactCtx_Probe2, semLatencyReactCtx_Probe2, 'cmap', blueCmap);
    [hl2,~] =boundedline(timescaleSeg, meanLatencyAvoidCtx_Probe2, semLatencyAvoidCtx_Probe2, 'cmap', orangeCmap);
    zeroXDottedLine;
    hold off;
    xlabel('time (sec)');
    title('Probe 2 mean')
    legend([hl1, hl2], {'react', 'avoid'});
    legend('boxoff')
    set(gcf,'color','w'); %set figure background white

end

function plotMeanOverUnitsOpto(allUnitsMeanZFR_optoTest, timescaleSeg)
 %combined ctx or non-ctx units avoid vs react mean responses
    [meanUnits, semUnits] = grpstats(allUnitsMeanZFR_optoTest,[],{'mean' 'sem'}); %mean across columns/timepoints
    figure;
    hold on;
    boundedline(timescaleSeg, meanUnits, semUnits, 'cmap', [0 0 0])
    ylabel('mean z-scored FR');
    xlabel('time (sec)');
    title('mean response at opto, across units')

    %colored rectangle patch when opto on, separate offset patch for ramp
    %down to apply gradient color elsewhere (or just leave off since matlab
    %groups patches on SVG export, & use graphics program to create)
    optoDuration = 2.0;
%     optoRampDown = 0.5;
    yLims = get(gca, 'YLim');
    y = [yLims(1) yLims(2) yLims(2) yLims(1)];
    x = [0 0 optoDuration optoDuration];
    patch(x, y, [0 0.5 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
%     x = [optoDuration+0.01 optoDuration+0.01 optoDuration+optoRampDown optoDuration+optoRampDown];
%     patch(x, y, [0 0.5 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
    hold off;
    set(gcf,'color','w'); %set figure background white

end

function plotMeanOverUnitsFromStruct(allData, removeNull)
% function to plot mean firing rate aligned to events, across units in a given region

    nullRemovalYaxisMultiplier = 0.5; %since null will remove activity, optionally magnify
    % make structure to iterate brain regions in a loop to plot all at once
    params(1).probeName = 'Probe1'; params(1).anatName = 'CTX'; params(1).labelName = 'PFC'; params(1).atCue = true; params(1).cmapIdx = 2; params(1).ylim = [-0.2 0.5]; %[-0.2 0.5] without null removal
    params(2).probeName = 'Probe1'; params(2).anatName = 'CTX'; params(2).labelName = 'PFC'; params(2).atCue = false; params(2).cmapIdx = 2; params(2).ylim = [-0.2 0.5];
    params(3).probeName = 'Probe0'; params(3).anatName = 'CTX'; params(3).labelName = 'M1'; params(3).atCue = true; params(3).cmapIdx = 1; params(3).ylim = [-0.2 1.2]; %[-0.2 1.2]
    params(4).probeName = 'Probe0'; params(4).anatName = 'CTX'; params(4).labelName = 'M1'; params(4).atCue = false; params(4).cmapIdx = 1; params(4).ylim = [-0.2 1.2];
    params(5).probeName = 'Probe0'; params(5).anatName = 'BS'; params(5).labelName = 'TH/HY'; params(5).atCue = true; params(5).cmapIdx = 3; params(5).ylim = [-0.3 0.3]; %[-0.3 0.3]
    params(6).probeName = 'Probe0'; params(6).anatName = 'BS'; params(6).labelName = 'TH/HY'; params(6).atCue = false; params(6).cmapIdx = 3; params(6).ylim = [-0.3 0.3];
    params(7).probeName = 'Probe2'; params(7).anatName = 'HB'; params(7).labelName = 'HB'; params(7).atCue = true; params(7).cmapIdx = 4; params(7).ylim = [-0.3 1.0]; %[-0.3 1.0]
    params(8).probeName = 'Probe2'; params(8).anatName = 'HB'; params(8).labelName = 'HB'; params(8).atCue = false; params(8).cmapIdx = 4; params(8).ylim = [-0.3 1.0];

    % option to only plot particular region (ex. PFC for Fig. 4f)    
%     params(1).probeName = 'Probe1'; params(1).anatName = 'CTX'; params(1).labelName = 'PFC'; params(1).atCue = false; params(1).cmapIdx = 2; params(1).ylim = [-0.05 0.1]*nullRemovalYaxisMultiplier; %[-0.2 0.5] without null removal

    f = figure; tiledlayout(4,2,"TileSpacing","compact");
    for p = 1:length(params)
    
        atCueOrVhex = params(p).atCue;
        probeName = params(p).probeName;
        anatName = params(p).anatName;
        labelName = params(p).labelName;
        cmapIdx = params(p).cmapIdx;
    
        if atCueOrVhex %true for cue timing
            if removeNull
                condAvoid = 'allUnitsTrialsAvoidAtCueNullMean'; %#ok<*NASGU> 
                condReact = 'allUnitsTrialsReactAtCueNullMean'; 
                condShuf1 = 'allUnitsTrialsShufAtCue1NullMean';
                condOpto = 'allUnitsTrialsOptoAtCueNullMean';
                titleStr = [labelName '(-null) at cue'];
                params(p).ylim = params(p).ylim.*nullRemovalYaxisMultiplier;
            else
                condAvoid = 'allUnitsTrialsAvoidAtCueMean'; 
                condReact = 'allUnitsTrialsReactAtCueMean'; 
                condShuf1 = 'allUnitsTrialsShufAtCue1Mean';
                condOpto = 'allUnitsTrialsOptoAtCueMean';
                titleStr = [labelName ' at cue'];
            end
        else %else @ VHEx timing
            if removeNull
                condAvoid = 'allUnitsTrialsAvoidNullMean';
                condReact = 'allUnitsTrialsReactNullMean';
                condShuf1 = 'allUnitsTrialsShuf1NullMean';
                condOpto = 'allUnitsTrialsOptoNullMean';
                titleStr = [labelName '(-null) at VHEx'];
                params(p).ylim = params(p).ylim.*nullRemovalYaxisMultiplier;
            else
                condAvoid = 'allUnitsTrialsAvoidMean';
                condReact = 'allUnitsTrialsReactMean';
                condShuf1 = 'allUnitsTrialsShuf1Mean';
                condOpto = 'allUnitsTrialsOptoMean';
                titleStr = [labelName ' at VHEx'];
            end
        end

        numMiceProbes = size(allData.(probeName).(anatName).(condAvoid),2);
        plotWindow = [-2 1]; % [-1 2] for cue, [-2 1] for vhex
        plotSamps = [find(plotWindow(1)>=allData.timescaleSeg,1,'last') find(plotWindow(2)>=allData.timescaleSeg,1,'last')];
        numPlotSamps = plotSamps(2) - plotSamps(1) + 1;
        timescalePlot = allData.timescaleSeg(plotSamps(1):plotSamps(2));
        %other constants
        fontSz = 16;
        blueCmap = [0 0.4470 0.7410];
        orangeCmap = [0.8500 0.3250 0.0980];
        cmapLine = allData.lineCmaps(cmapIdx,:);
        numRepeats = 10;
        smoothingSamples = [5 0]; %ex. smooth window 10 is 10*20ms PSTH bins = 200ms; for causal half-kernel use [X 0]
    
        allMiceUnitsMeanAvoidData = []; 
        allMiceUnitsMeanReactData = [];
        allMiceUnitsMeanShuf1Data = [];
        allMiceUnitsMeanOptoData = [];

        % for each mouse/repeat, compute Euclidean distance from subspace trial means, then average over trial sampling repeats
        miceToPlot = 1:3;
        for i = 1:numMiceProbes 
%         for i = miceToPlot %     
            thisMouseAvoidData = allData.(probeName).(anatName).(condAvoid){i}; %unpack to get per-mouse cells that include all sample repeats
            thisMouseReactData = allData.(probeName).(anatName).(condReact){i};
            thisMouseShuf1Data = allData.(probeName).(anatName).(condShuf1){i};
            thisMouseOptoData = allData.(probeName).(anatName).(condOpto){i};
            numUnits = size(thisMouseAvoidData{1},1);

            % optionally remove VGAT positively-modulated opto neurons for analysis of residual activity in PFC (Fig. 4f) - also uncomment below
%             optoRawData = allData.(probeName).(anatName).allUnitsTrialsOptoAtCueMean{i}; % use cue-aligned data since opto on at same time
%             sortWindowBefore = [-0.5 -0.1];
%             sortWindowAfter = [0.1 0.5];
%             sortSampsBefore = [find(sortWindowBefore(1)>=allData.timescaleSeg,1,'last') find(sortWindowBefore(2)>=allData.timescaleSeg,1,'last')];
%             sortSampsAfter = [find(sortWindowAfter(1)>=allData.timescaleSeg,1,'last') find(sortWindowAfter(2)>=allData.timescaleSeg,1,'last')];
%             dataBefore = mean(optoRawData{1}(:,sortSampsBefore(1):sortSampsBefore(2)),2); %just use first opto on
%             dataAfter = mean(optoRawData{1}(:,sortSampsAfter(1):sortSampsAfter(2)),2); 
%             optoModulation = dataAfter - dataBefore;
%             posOptoModIdx = optoModulation>1.5; %threshold is empirical
%             vgatRemovedIdx = ~posOptoModIdx; %make logical index for remaining non-VGAT neurons 
%             numUnits = length(find(vgatRemovedIdx));
            
            thisRptAvoidData = zeros(numRepeats, numUnits, numPlotSamps);
            thisRptReactData = zeros(numRepeats, numUnits, numPlotSamps);
            thisRptShuf1Data = zeros(numRepeats, numUnits, numPlotSamps);
            thisRptOptoData = zeros(numRepeats, numUnits, numPlotSamps);

            for rpt = 1:numRepeats
                thisRptAvoidData(rpt,:,:) = thisMouseAvoidData{rpt}(:,plotSamps(1):plotSamps(2)); %unpack one more time to get each repeat cell - now this cell for each trial mean
                thisRptReactData(rpt,:,:) = thisMouseReactData{rpt}(:,plotSamps(1):plotSamps(2));
                thisRptShuf1Data(rpt,:,:) = thisMouseShuf1Data{rpt}(:,plotSamps(1):plotSamps(2));
                thisRptOptoData(rpt,:,:) = thisMouseOptoData{rpt}(:,plotSamps(1):plotSamps(2));
                % for Fig. 4f residual activity:
%                 thisRptAvoidData(rpt,:,:) = thisMouseAvoidData{rpt}(vgatRemovedIdx,plotSamps(1):plotSamps(2)); %unpack one more time to get each repeat cell - now this cell for each trial mean
%                 thisRptReactData(rpt,:,:) = thisMouseReactData{rpt}(vgatRemovedIdx,plotSamps(1):plotSamps(2));
%                 thisRptShuf1Data(rpt,:,:) = thisMouseShuf1Data{rpt}(vgatRemovedIdx,plotSamps(1):plotSamps(2));
%                 thisRptOptoData(rpt,:,:) = thisMouseOptoData{rpt}(vgatRemovedIdx,plotSamps(1):plotSamps(2));

            end %end sampling repeats

            % now take mean across repeats for each mouse & smooth:
            thisMouseMeanAvoidData = squeeze(mean(thisRptAvoidData,1)); 
            thisMouseMeanReactData = squeeze(mean(thisRptReactData,1));
            thisMouseMeanShuf1Data = squeeze(mean(thisRptShuf1Data,1));
            thisMouseMeanOptoData = squeeze(mean(thisRptOptoData,1));
            thisMouseMeanAvoidData =  smoothdata(thisMouseMeanAvoidData, 2, 'gaussian', smoothingSamples); 
            thisMouseMeanReactData =  smoothdata(thisMouseMeanReactData, 2, 'gaussian', smoothingSamples); 
            thisMouseMeanShuf1Data =  smoothdata(thisMouseMeanShuf1Data, 2, 'gaussian', smoothingSamples); 
            thisMouseMeanOptoData =  smoothdata(thisMouseMeanOptoData, 2, 'gaussian', smoothingSamples); 

            % and finally take mean across all units to get 1-D array for each mouse (for this iteration of params)
            thisMouseAvoidAcrossUnits = squeeze(mean(thisMouseMeanAvoidData,1));
            thisMouseReactAcrossUnits = squeeze(mean(thisMouseMeanReactData,1));
            thisMouseShuf1AcrossUnits = squeeze(mean(thisMouseMeanShuf1Data,1));
            thisMouseOptoAcrossUnits = squeeze(mean(thisMouseMeanOptoData,1));
            allMiceUnitsMeanAvoidData = [allMiceUnitsMeanAvoidData; thisMouseAvoidAcrossUnits]; 
            allMiceUnitsMeanReactData = [allMiceUnitsMeanReactData; thisMouseReactAcrossUnits];
            allMiceUnitsMeanShuf1Data = [allMiceUnitsMeanShuf1Data; thisMouseShuf1AcrossUnits];
            allMiceUnitsMeanOptoData = [allMiceUnitsMeanOptoData; thisMouseOptoAcrossUnits];
            %optionally subtract median at beginning for comparison of residual activity across Vgat & Vglut1 (Fig. 4f)
%             allMiceUnitsMeanAvoidData = [allMiceUnitsMeanAvoidData; thisMouseAvoidAcrossUnits - median(thisMouseAvoidAcrossUnits(1:20))]; 
%             allMiceUnitsMeanReactData = [allMiceUnitsMeanReactData; thisMouseReactAcrossUnits - median(thisMouseReactAcrossUnits(1:20))]; 
%             allMiceUnitsMeanShuf1Data = [allMiceUnitsMeanShuf1Data; thisMouseShuf1AcrossUnits - median(thisMouseShuf1AcrossUnits(1:20))]; 
%             allMiceUnitsMeanOptoData = [allMiceUnitsMeanOptoData; thisMouseOptoAcrossUnits - median(thisMouseOptoAcrossUnits(1:20))]; 

        end %end numMice
 
       nexttile; hold on;
        % optionally plot individual mice to see variability:
        for i = 1:size(allMiceUnitsMeanOptoData,1)
            plot(timescalePlot, squeeze(allMiceUnitsMeanOptoData(i,:)), 'Color', [0 0 0], 'LineWidth', 0.1)
        end

        [mean1, sem1] = grpstats(allMiceUnitsMeanAvoidData,[],{'mean' 'sem'}); %mean across columns/timepoints
        [h1(p),~] = boundedline(timescalePlot, mean1, sem1, 'cmap', orangeCmap, 'alpha','y', 'transparency', 0.5, 'LineWidth', 3);
        [mean2, sem2] = grpstats(allMiceUnitsMeanReactData,[],{'mean' 'sem'}); %mean across columns/timepoints
        [h2(p),~] = boundedline(timescalePlot, mean2, sem2, 'cmap', blueCmap, 'alpha','y', 'transparency', 0.5, 'LineWidth', 3);
        [mean3, sem3] = grpstats(allMiceUnitsMeanOptoData,[],{'mean' 'sem'}); %mean across columns/timepoints
        [h3(p),~] = boundedline(timescalePlot, mean3, sem3, 'cmap', [0 0 0], 'alpha','y', 'transparency', 0.5, 'LineWidth', 3);
        ylim(params(p).ylim)
%         ylim([-0.05 0.1]) %for VGAT-ChR2 residual opto activity Fig. 4f
%         xlim([timescalePlot(1)+0.5 timescalePlot(end)]) %adjust for smoothing edge effects
        zeroXDottedLine;
        makepretty;
        hold off;
        title(titleStr)
        if p==5
            ylabel('mean z-score FR across neurons', 'FontSize', fontSz)
        elseif p==8
            xlabel('time (sec)', 'FontSize', fontSz)
        end

    end %end params loop

    f.Position(3:4) = [1200 900]; %set figure size appropriate for plots
    movegui(f,'center')
%     legend([h1(1), h1(3), h1(5), h1(7), h2], {params(1).labelName, params(3).labelName, params(5).labelName, params(7).labelName,'shuffle'}, 'Box', 'off', 'Location','northeast');

end

function [meanFR, semFR, timescale] = meanSemFiringRateAcrossEvents(firingRates, eventTimes, psthBinSize, eventWindow)
% function to return mean & SEM of firing rates across neurons, around multiple instances of an event
% NOTE: this discards the variance across neurons in a single trial
% eventWindow should be a pair of how many seconds to use around the event time; negative is before, positive is after; ex [-1 0.5]
%     numUnits = size(allLocalFR,2);
    timeBefore = eventWindow(1); timeAfter = eventWindow(2);
    timescale = timeBefore:psthBinSize:timeAfter-psthBinSize;
    psthTimes = 0:psthBinSize:psthBinSize*size(firingRates,1);
    psthSamplesBefore = round(timeBefore/psthBinSize);
    psthSamplesAfter = round(timeAfter/psthBinSize);
    allLocalFrEvents = zeros(length(eventTimes), psthSamplesAfter-psthSamplesBefore); %put time in column dimension
    for i = 1:length(eventTimes)
        psthSample = find(psthTimes>eventTimes(i), 1, 'first');
        allLocalFrEvents(i,:) = mean(firingRates(psthSample+psthSamplesBefore:psthSample+psthSamplesAfter-1,:),2);
    end
    [meanFR,semFR] = grpstats(allLocalFrEvents,[],{'mean' 'sem'}); %stats are across columns/timepoints
end

function [behTrialsZ, timescaleSeg] = getBehaviorZscoreEvents(behaviorData, behaviorTimescale, eventTimes, eventWindow)
% function to collect z-scored behavior measure in matrix of windowed trials
% input is full behavior timeseries + timescale in seconds, and event times of interest + window in seconds
% output is matrix of [trials x timepointsInWindow]
    if isempty(eventTimes)
        behTrialsZ = [];
        timescaleSeg = [];
        return
    end

    sampTime = behaviorTimescale(2) - behaviorTimescale(1);
    timescaleSeg = eventWindow(1):sampTime:eventWindow(2)-sampTime;
    samplesBefore = round(eventWindow(1)/sampTime);
    samplesAfter = round(eventWindow(2)/sampTime);
    behTrialsZ = zeros(length(eventTimes),length(timescaleSeg));

    % first normalize data to z-scored version
    behZ = normalize(behaviorData);
    
    % now collect across events
    for j = 1:length(eventTimes)
        currentSample = find(behaviorTimescale>eventTimes(j), 1, 'first');
        behTrialsZ(j,:) = behZ(currentSample+samplesBefore:currentSample+samplesAfter-1);
    end

end

function [ephysTrialsZ, timescaleSeg] = getEphysZscoreEvents(shankData, eventTimes, window)
% function to collect z-scored ephys in matrix of windowed trials
% input is full z-scored ephys timeseries and event times of interest + window in seconds
% output is matrix of [neurons x trials x timepointsInWindow]

    numUnits = length(shankData.stCell);
    psthBinSize = shankData.psthTimes(2) - shankData.psthTimes(1);
    psthTimes = 0:psthBinSize:psthBinSize*size(shankData.allLocalFR,1);
    timescaleSeg = window(1):psthBinSize:window(2)-psthBinSize;
    psthSamplesBefore = round(window(1)/psthBinSize);
    psthSamplesAfter = round(window(2)/psthBinSize);
    ephysTrialsZ = zeros(numUnits, length(eventTimes),length(timescaleSeg)); %neurons x trials x time

    for i = 1:numUnits
        for j = 1:length(eventTimes)
            psthSample = find(psthTimes>eventTimes(j), 1, 'first');
            ephysTrialsZ(i,j,:) = shankData.allDeltaFRzScore(psthSample+psthSamplesBefore:psthSample+psthSamplesAfter-1,i);
        end
    end
end

function [allUnitsMeanZFR, timescaleSeg] = getEphysMeanZscoreFREvent(shankData, eventTimes, window)
% get mean zscored FR around some events; zscore is across entire session % assume shankData already in desired sort order
% NOTE that this assumes window(1) is signed appropriately
    numUnits = length(shankData.stCell);
    psthBinSize = shankData.psthTimes(2) - shankData.psthTimes(1);
    psthTimes = 0:psthBinSize:psthBinSize*size(shankData.allLocalFR,1);
    timescaleSeg = window(1):psthBinSize:window(2)-psthBinSize;
    psthSamplesBefore = round(window(1)/psthBinSize);
    psthSamplesAfter = round(window(2)/psthBinSize);
    thisUnitAllEvents = zeros(length(eventTimes),length(timescaleSeg));
    allUnitsMeanZFR = zeros(numUnits,length(timescaleSeg));
    smoothWindow = 5; % how many PSTH bins in window for Gaussian smoothing below; ex. 20ms bins x 5 bins = 100ms; causal half kernel use [X 0]; keep on the order of neurotransmitter time constants
    for i = 1:numUnits %go in sort order
        for j = 1:length(eventTimes)
            psthSample = find(psthTimes>eventTimes(j), 1, 'first');
            thisUnitAllEvents(j,:) = shankData.allDeltaFRzScore(psthSample+psthSamplesBefore:psthSample+psthSamplesAfter-1,i);
            thisUnitAllEvents(j,:) = smoothdata(thisUnitAllEvents(j,:),2,'gaussian',smoothWindow); %smooth Z-score over time
        end
        allUnitsMeanZFR(i,:) = mean(thisUnitAllEvents,1); %mean across columns/timepoints
    end
end

function [allUnitsMeanZFR, timescaleSeg] = getEphysMeanZscoreFREventEqualTrials(shankData, eventTimes1, eventTimes2, window)
% get mean zscored FR around some events, but match # trials between event types
% NOTE that this assumes window(1) is signed appropriately
    numRepeats = 10; %since we do random sampling, repeat the processto allow averaging out sampling bias later
    for rpt = 1:numRepeats
        %use min of the two # of events to pick from each distribution
        numEvents1 = length(eventTimes1);
        numEvents2 = length(eventTimes2);
        if numEvents1 < numEvents2
            eventTimes = eventTimes1; %if less of current event types, just use all of them
        else
            s = RandStream('dsfmt19937','Seed','shuffle');
            randIdx = randperm(s, numEvents1, numEvents2); %else take randowm subset euqal to smaller # of trials
            eventTimes = eventTimes1(randIdx);
        end
    
        % now call original function to compute trial means:
        [allUnitsMeanZFR(rpt,:,:), timescaleSeg] = getEphysMeanZscoreFREvent(shankData, eventTimes, window);
    end
end

function [allUnitsMeanZFR, timescaleSeg] = getEphysMeanZscoreFREventEqualTrialsInclOpto(shankData, eventTimes1, minTrials, window)
% get mean zscored FR around some events, but match # trials between event types
% this function is different from above since it is based on minTrials betwenn opto/react/avoid, not just react/avoid
% NOTE that this assumes window(1) is signed appropriately
    numRepeats = 10; %since we do random sampling, repeat the processto allow averaging out sampling bias later
    for rpt = 1:numRepeats
        %use min of the two # of events to pick from each distribution
        numEvents1 = length(eventTimes1);
        if numEvents1 <= minTrials
            eventTimes = eventTimes1; %if less of current event types, just use all of them
        else
            s = RandStream('dsfmt19937','Seed','shuffle');
            randIdx = randperm(s, numEvents1, minTrials); %else take randowm subset euqal to smaller # of trials
            eventTimes = eventTimes1(randIdx);
        end
    
        % now call original function to compute trial means:
        [allUnitsMeanZFR(rpt,:,:), timescaleSeg] = getEphysMeanZscoreFREvent(shankData, eventTimes, window);
    end
end

function [allUnitsMeanFR, timescaleSeg] = getEphysMeanFREvent(shankData, eventTimes, window)
% get mean FR (Hz) around some events; % assume shankData already in desired sort order
    numUnits = length(shankData.stCell);
    psthBinSize = shankData.psthTimes(2) - shankData.psthTimes(1);
    psthTimes = 0:psthBinSize:psthBinSize*size(shankData.allLocalFR,1);
    timescaleSeg = window(1):psthBinSize:window(2)-psthBinSize;
    psthSamplesBefore = round(window(1)/psthBinSize);
    psthSamplesAfter = round(window(2)/psthBinSize);
    thisUnitAllEvents = zeros(length(eventTimes),length(timescaleSeg));
    allUnitsMeanFR = zeros(numUnits,length(timescaleSeg));
    for i = 1:numUnits %go in sort order
        for j = 1:length(eventTimes)
            psthSample = find(psthTimes>eventTimes(j), 1, 'first');
            thisUnitAllEvents(j,:) = shankData.allLocalFR(psthSample+psthSamplesBefore:psthSample+psthSamplesAfter-1,i);
        end
        allUnitsMeanFR(i,:) = mean(thisUnitAllEvents,1); %mean across columns/timepoints
    end
end

function [allUnitsMeanPercentFR, timescaleSeg] = getEphysMeanPercentFREvent(shankData, eventTimes, window, sortIdxNegOnly)
% get mean FR (Hz) as a percentage of baseline before some events; use only negatively modulated neurons
    sortedLocalFR = shankData.allLocalFR(:,sortIdxNegOnly);
    numUnits = size(sortedLocalFR,2);
    psthBinSize = shankData.psthTimes(2) - shankData.psthTimes(1);
    psthTimes = 0:psthBinSize:psthBinSize*size(sortedLocalFR,1);
    timescaleSeg = window(1):psthBinSize:window(2)-psthBinSize;
    psthSamplesBefore = round(window(1)/psthBinSize);
    psthSamplesAfter = round(window(2)/psthBinSize);
    thisUnitBeforeAllEvents = zeros(length(eventTimes),length(timescaleSeg)/2);
    thisUnitAfterAllEvents = zeros(length(eventTimes),length(timescaleSeg)/2);
    allUnitsMeanPercentFR = zeros(numUnits,1);
    for i = 1:numUnits %go in sort order
        for j = 1:length(eventTimes)
            psthSample = find(psthTimes>eventTimes(j), 1, 'first');
            thisUnitBeforeAllEvents(j,:) = sortedLocalFR(psthSample+psthSamplesBefore:psthSample-1,i);
            thisUnitAfterAllEvents(j,:) = sortedLocalFR(psthSample+1:psthSample+psthSamplesAfter,i);
        end
        currentPercentFR = mean(mean(thisUnitAfterAllEvents,1)) ./ mean(mean(thisUnitBeforeAllEvents,1)); %mean across columns/timepoints
        if currentPercentFR == Inf
            allUnitsMeanPercentFR(i) = NaN; %if not enough spikes before, ignore this unit (divide by 0 = Inf) by setting NaN
        else
            allUnitsMeanPercentFR(i) = currentPercentFR;
        end
    end
end

function [avoidMeanZFR, reactMeanZFR, timescaleSeg] = getAvoidVsReactMeanZscoreFREvent(shankData, behaviorTimescale, idxData, window)
% plot meanFR around some events either zscored or not; zscore is across entire session % assume shankData already in desired sort order
% separate trials by avoid/react, but align to another event hardcoded below (e.g. beep on)

    numUnits = length(shankData.stCell);
    psthBinSize = shankData.psthTimes(2) - shankData.psthTimes(1);
    psthTimes = 0:psthBinSize:psthBinSize*size(shankData.allLocalFR,1);
    timescaleSeg = window(1):psthBinSize:window(2)-psthBinSize;
    psthSamplesBefore = round(window(1)/psthBinSize);
    psthSamplesAfter = round(window(2)/psthBinSize);

    %first get all event times of interest, plus an index for whether the time is for avoid or react trial
    nonOptoCurrentTrial = 1;
    avoidIdx = [];
    beepTimes = [];
    for i = 1:length(idxData.jumpIndeces) %loop over all extensions for the session
        if ~idxData.jumpOptoIndeces(i) %only non-opto trials for now
            nearestBeepIdx = idxData.beepOnIndeces(find(idxData.beepOnIndeces < idxData.jumpIndeces(i), 1, 'last')); %need to find nearest beep on
            beepTimes(nonOptoCurrentTrial) = behaviorTimescale(nearestBeepIdx); 
            nonOptoCurrentTrial = nonOptoCurrentTrial + 1;
            if idxData.jumpAntIndeces(i)
                avoidIdx(nonOptoCurrentTrial) = 1; %construct avoid/react indeces (avoid=1) for only non-opto extensions
            else
                avoidIdx(nonOptoCurrentTrial) = 0;
            end
        end
    end

    thisUnitAvoidEvents = zeros(length(find(avoidIdx)), length(timescaleSeg));
    thisUnitReactEvents = zeros(length(find(~avoidIdx)), length(timescaleSeg));
    avoidMeanZFR = zeros(numUnits,length(timescaleSeg));
    reactMeanZFR = zeros(numUnits,length(timescaleSeg));

    for i = 1:numUnits %go in sort order
        avoidEventIdx = 1;
        reactEventIdx = 1;
        for j = 1:length(beepTimes)
            psthSample = find(psthTimes>beepTimes(j), 1, 'first');
            if avoidIdx(j)
                thisUnitAvoidEvents(avoidEventIdx,:) = shankData.allDeltaFRzScore(psthSample+psthSamplesBefore:psthSample+psthSamplesAfter-1,i);
                avoidEventIdx = avoidEventIdx + 1;
            else
                thisUnitReactEvents(reactEventIdx,:) = shankData.allDeltaFRzScore(psthSample+psthSamplesBefore:psthSample+psthSamplesAfter-1,i);
                reactEventIdx = reactEventIdx + 1;
            end
        end
        avoidMeanZFR(i,:) = mean(thisUnitAvoidEvents,1); %mean across columns/timepoints
        reactMeanZFR(i,:) = mean(thisUnitReactEvents,1);
    end
end

function [allUnitsMeanZFRShuf, timescaleSeg] = getEphysMeanZscoreFRShuffleEven(shankData, eventTimes1, eventTimes2, window)
% get mean from shuffled trials for control; zscore is across entire session % assume shankData already in desired sort order
% this 'Even' version uses equal numbers of events from Times 1 & 2, which may mask differences due to trial #s (see alternative in function below) 

    %use min of the two # of events to pick from each distribution
    numEvents1 = length(eventTimes1);
    numEvents2 = length(eventTimes2);
    numShufEvents = min(numEvents1, numEvents2)*2;

    s = RandStream('dsfmt19937','Seed','shuffle');
    randIdx1 = randperm(s, numEvents1, numShufEvents/2); %take equal # random trials from each type of event
    randIdx2 = randperm(s, numEvents2, numShufEvents/2);
    eventTimesShuf1 = eventTimes1(randIdx1);
    eventTimesShuf2 = eventTimes2(randIdx2);


    numUnits = length(shankData.stCell);
    psthBinSize = shankData.psthTimes(2) - shankData.psthTimes(1);
    psthTimes = 0:psthBinSize:psthBinSize*size(shankData.allLocalFR,1);
    timescaleSeg = window(1):psthBinSize:window(2)-psthBinSize;
    psthSamplesBefore = round(window(1)/psthBinSize);
    psthSamplesAfter = round(window(2)/psthBinSize);

    thisUnitAllEvents1 = zeros(length(eventTimesShuf1),length(timescaleSeg));
    thisUnitAllEvents2 = zeros(length(eventTimesShuf2),length(timescaleSeg));

    allUnitsMeanZFRShuf = zeros(numUnits,length(timescaleSeg));
    for i = 1:numUnits %go in sort order through all units
        for j = 1:length(eventTimesShuf1)
            psthSample = find(psthTimes>eventTimesShuf1(j), 1, 'first');
            thisUnitAllEvents1(j,:) = shankData.allDeltaFRzScore(psthSample+psthSamplesBefore:psthSample+psthSamplesAfter-1,i);
        end
        for j = 1:length(eventTimesShuf2)
            psthSample = find(psthTimes>eventTimesShuf2(j), 1, 'first');
            thisUnitAllEvents2(j,:) = shankData.allDeltaFRzScore(psthSample+psthSamplesBefore:psthSample+psthSamplesAfter-1,i);
        end
        thisUnitAllEvents = cat(1, thisUnitAllEvents1, thisUnitAllEvents2); %concatentate random trials across both event types
        allUnitsMeanZFRShuf(i,:) = mean(thisUnitAllEvents,1); %mean across columns/timepoints
    end
end

function [allUnitsMeanZFRShuf, timescaleSeg] = getEphysMeanZscoreFRShuffleEqualTrials(shankData, eventTimes1, eventTimes2, window)
% get mean from shuffled trials for control
% this 'EqualTrials' version uses equal numbers of events from Times 1 & 2 to make a shuffle distribution with # trials equal to the lesser # trials 

    numShuffles = 10; 
    for sh = 1:numShuffles
        %use min of the two # of events to pick from each distribution, to match normal EqualTrials sampling procedure
        numEvents1 = length(eventTimes1);
        numEvents2 = length(eventTimes2);
        numShufEvents = min(numEvents1, numEvents2);
    
        if mod(numShufEvents,2)==0 %if even
            numEventsToChoose1 = numShufEvents/2;
            numEventsToChoose2 = numShufEvents/2;
        else %if odd
            numEventsToChoose1 = (numShufEvents+1)/2;
            numEventsToChoose2 = (numShufEvents-1)/2;
        end
    
        s = RandStream('dsfmt19937','Seed','shuffle');
        randIdx1 = randperm(s, numEvents1, numEventsToChoose1); %take equal # random trials from each type of event
        randIdx2 = randperm(s, numEvents2, numEventsToChoose2);
        eventTimesShuf = [eventTimes1(randIdx1) eventTimes2(randIdx2)]; %1/2 times from event1, 1/2 times from event2
    
%         % option to choose randomly from all trials rather than 1/2 of each type, which may unfairly mask trial varaibility in one type but doesn't seem to make a difference:
%         concatEventTimes = [eventTimes1 eventTimes2];
%         randIdx = randperm(s, length(concatEventTimes), numShufEvents);
%         eventTimesShuf = concatEventTimes(randIdx);

        % now call original function to compute trial means:
        [allUnitsMeanZFRShuf(sh,:,:), timescaleSeg] = getEphysMeanZscoreFREvent(shankData, eventTimesShuf, window);
    end
end

function [allUnitsMeanZFRShuf, timescaleSeg] = getEphysMeanZscoreFRShuffleEqualTrialsInclOpto(shankData, eventTimes1, eventTimes2, minTrials, window)
% get mean from shuffled trials for control; this version includes minTrials that takes optoTrials into account
% this 'EqualTrials' version uses equal numbers of events from Times 1 & 2 to make a shuffle distribution with # trials equal to the lesser # trials 

    numShuffles = 10; 
    for sh = 1:numShuffles
        %use min of the two # of events to pick from each distribution, to match normal EqualTrials sampling procedure
        numEvents1 = length(eventTimes1);
        numEvents2 = length(eventTimes2);
        numShufEvents = min([numEvents1 numEvents2 minTrials]);
    
        if mod(numShufEvents,2)==0 %if even
            numEventsToChoose1 = numShufEvents/2;
            numEventsToChoose2 = numShufEvents/2;
        else %if odd
            numEventsToChoose1 = (numShufEvents+1)/2;
            numEventsToChoose2 = (numShufEvents-1)/2;
        end
    
        s = RandStream('dsfmt19937','Seed','shuffle');
        randIdx1 = randperm(s, numEvents1, numEventsToChoose1); %take equal # random trials from each type of event
        randIdx2 = randperm(s, numEvents2, numEventsToChoose2);
        eventTimesShuf = [eventTimes1(randIdx1) eventTimes2(randIdx2)]; %1/2 times from event1, 1/2 times from event2
    
%         % option to choose randomly from all trials rather than 1/2 of each type, which may unfairly mask trial varaibility in one type but doesn't seem to make a difference:
%         concatEventTimes = [eventTimes1 eventTimes2];
%         randIdx = randperm(s, length(concatEventTimes), numShufEvents);
%         eventTimesShuf = concatEventTimes(randIdx);

        % now call original function to compute trial means:
        [allUnitsMeanZFRShuf(sh,:,:), timescaleSeg] = getEphysMeanZscoreFREvent(shankData, eventTimesShuf, window);
    end
end

function [allUnitsMeanZFRShuf, timescaleSeg] = getEphysMeanZscoreFRShuffle(shankData, eventTimes1, eventTimes2, useEvents1, window)
% get mean from shuffled trials for control, but need to match unequal # of avoid/react trials, by choosing an equal number from each distribution to get to the target # of trials
% useEvents1 decides whether to sample the # of events from eventTime1 vs eventTimes2; this allows making 2 distributions that match respective trial numbers

%NOTE: must use same smoothWindow as getEphysMeanZscoreFREvent() for fair comparison
smoothWindow = 5; % how many PSTH bins in window for Gaussian smoothing below; ex. 20ms bins x 5 bins = 100ms; causal half kernel use [X 0]; keep on the order of neurotransmitter time constants

    numShuffles = 1; 
    for sh = 1:numShuffles
        if useEvents1
            numShufEvents = length(eventTimes1); 
        else
            numShufEvents = length(eventTimes2); 
        end
        
        if mod(numShufEvents,2)==0 %if even
            numEventsToChoose1 = numShufEvents/2;
            numEventsToChoose2 = numShufEvents/2;
        else %if odd
            numEventsToChoose1 = (numShufEvents+1)/2;
            numEventsToChoose2 = (numShufEvents-1)/2;
        end

        s = RandStream('dsfmt19937','Seed','shuffle');
        ei1 = datasample(s, eventTimes1, numEventsToChoose1); %draw randomly, with replacement, from eventTimes1
        ei2 = datasample(s, eventTimes2, numEventsToChoose2);
        eventTimesShuf = [ei1 ei2];
    
        numUnits = length(shankData.stCell);
        psthBinSize = shankData.psthTimes(2) - shankData.psthTimes(1);
        psthTimes = 0:psthBinSize:psthBinSize*size(shankData.allLocalFR,1);
        timescaleSeg = window(1):psthBinSize:window(2)-psthBinSize;
        psthSamplesBefore = round(window(1)/psthBinSize);
        psthSamplesAfter = round(window(2)/psthBinSize);
    
        thisUnitAllEvents = zeros(length(eventTimesShuf),length(timescaleSeg));
        allUnitsMeanZFRThisShuf = zeros(numShuffles,numUnits,length(timescaleSeg));
        for i = 1:numUnits %go in sort order through all units
            for j = 1:length(eventTimesShuf)
                psthSample = find(psthTimes>eventTimesShuf(j), 1, 'first');
                thisUnitAllEvents(j,:) = shankData.allDeltaFRzScore(psthSample+psthSamplesBefore:psthSample+psthSamplesAfter-1,i);
                thisUnitAllEvents(j,:) = smoothdata(thisUnitAllEvents(j,:),2,'gaussian',smoothWindow); %smooth Z-score over time
            end
            allUnitsMeanZFRThisShuf(sh,i,:) = mean(thisUnitAllEvents,1); %mean across columns/trials
        end
    end
    allUnitsMeanZFRShuf = squeeze(mean(allUnitsMeanZFRThisShuf,1)); %mean across shuffles
end

function [specificAvoidIndeces, nonSpecificAvoidIndeces] = findAvoidAtVhexNeurons(anatData, timescaleSeg, label)
    %find "avoid neurons" that are active before avoid VHEx versus react
    %inputs are rpt x units x timeSeg 
    % specificAvoidIndeces are neurons that are more active before avoid VHEx AND relatively inactive before react (perhaps specifically recruited by cue)
    % nonSpecificAvoidIndeces are just neurons that are more active before avoid VHEx 

    avoidWindow = [-1 0];
    avoidMeanZFR = anatData.allUnitsAvoidMeanZFR;
    reactMeanZFR = anatData.allUnitsReactMeanZFR;
    numUnits = size(avoidMeanZFR, 2);
    avoidSamps = [find(avoidWindow(1)>=timescaleSeg,1,'last') find(avoidWindow(2)>=timescaleSeg,1,'last')];
    nonPersistentCmap = [0 0 0];
    persistentCmap = [0 0.7 0.5];

    %first get means over repeats:
    avoidMean = squeeze(mean(avoidMeanZFR(:,:,avoidSamps(1):avoidSamps(2)),1));
    reactMean = squeeze(mean(reactMeanZFR(:,:,avoidSamps(1):avoidSamps(2)),1));
    diffMean = avoidMean - reactMean;
    diffMeanOverWindow = mean(diffMean,2);

    [~,sortIdx] = sort(diffMeanOverWindow, 'descend');
    avoidMeanSorted = avoidMean(sortIdx,:);
    reactMeanSorted = reactMean(sortIdx,:);
    diffMeanSorted = diffMean(sortIdx,:);

    specificAvoidIndeces = false(numUnits,1); %logical indeces
    nonSpecificAvoidIndeces = false(numUnits,1);
    nonSpecificAntiAvoidIndeces = false(numUnits,1);
    for i = 1:numUnits
        if (mean(diffMeanSorted(i,:)) > 0.1)
            nonSpecificAvoidIndeces(i) = true;
            if (mean(reactMeanSorted(i,:)) < 0.1)
                specificAvoidIndeces(i) = true;
            end
        elseif (mean(diffMeanSorted(i,:)) < -0.1)
            nonSpecificAntiAvoidIndeces(i) = true;
        end
    end
    percentSpecificAvoidNeurons = length(find(specificAvoidIndeces)) / numUnits;
    percentNonSpecificAvoidNeurons = length(find(nonSpecificAvoidIndeces)) / numUnits;
    percentNonSpecificAntiAvoidNeurons = length(find(nonSpecificAntiAvoidIndeces)) / numUnits;
    
%     lastAvoidNeuron = find(nonSpecificAvoidIndeces,1,'last');
%     lastAntiAvoidNeuron = find(nonSpecificAntiAvoidIndeces,1,'first');
%     figure; tiledlayout(1,3, "TileSpacing","tight");
%     nexttile; hold on; imagesc(avoidMeanSorted); clim([-1 1]); xLim = get(gca, 'xlim'); plot(xLim, [lastAvoidNeuron lastAvoidNeuron],'k--', 'LineWidth',2); plot(xLim, [lastAntiAvoidNeuron lastAntiAvoidNeuron],'r--', 'LineWidth',2);
%     title(label); axis tight;
%     nexttile; hold on; imagesc(reactMeanSorted); clim([-1 1]); plot(xLim, [lastAvoidNeuron lastAvoidNeuron],'k--', 'LineWidth',2); plot(xLim, [lastAntiAvoidNeuron lastAntiAvoidNeuron],'r--', 'LineWidth',2); axis tight;
%     nexttile; hold on; imagesc(diffMeanSorted); clim([-1 1]); plot(xLim, [lastAvoidNeuron lastAvoidNeuron],'k--', 'LineWidth',2); plot(xLim, [lastAntiAvoidNeuron lastAntiAvoidNeuron],'r--', 'LineWidth',2);
%     title(['% avoid neurons = ' num2str(percentNonSpecificAvoidNeurons*100)]); axis tight;
%     keyboard

    %now compare mean 'tuning' to other relevant events:
    cueWindow = [-0.5 2];
    cueSamps = [find(cueWindow(1)>=timescaleSeg,1,'last') find(cueWindow(2)>=timescaleSeg,1,'last')];
    timescaleCuePlot = timescaleSeg(cueSamps(1):cueSamps(2));
    cueData = squeeze(mean(anatData.allUnitsAvoidAtCueMeanZFR(:,:,cueSamps(1):cueSamps(2)),1)); %mean over repeats & window
    cueDataSorted = cueData(sortIdx,:); %use same sort as above so that xIndeces apply
%     figure; imagesc(cueDataSorted); clim([-1 1]);

    coldAirWindow = [-1 1];
    coldSamps = [find(coldAirWindow(1)>=timescaleSeg,1,'last') find(coldAirWindow(2)>=timescaleSeg,1,'last')];
    timescaleColdPlot = timescaleSeg(coldSamps(1):coldSamps(2));
    coldData = anatData.allUnitsColdAirOnMeanZFR(:,coldSamps(1):coldSamps(2)); %no repeats for cold air on data
    coldDataSorted = coldData(sortIdx,:); %use same sort as above so that xIndeces apply

    vhexWindow = [-1 3];
    vhexSamps = [find(vhexWindow(1)>=timescaleSeg,1,'last') find(vhexWindow(2)>=timescaleSeg,1,'last')];
    timescaleVhexPlot = timescaleSeg(vhexSamps(1):vhexSamps(2));
    postReactData = squeeze(mean(anatData.allUnitsReactMeanZFR(:,:,vhexSamps(1):vhexSamps(2)),1)); %mean over repeats & window
    postReactDataSorted = postReactData(sortIdx,:); %use same sort as above so that xIndeces apply

    postAvoidData = squeeze(mean(anatData.allUnitsAvoidMeanZFR(:,:,vhexSamps(1):vhexSamps(2)),1)); %mean over repeats & window
    postAvoidDataSorted = postAvoidData(sortIdx,:); %use same sort as above so that xIndeces apply

    figure;
    tiledlayout(1,4,'TileSpacing', 'tight');
    nexttile; hold on;
    [mean1, sem1] = grpstats(cueDataSorted(nonSpecificAvoidIndeces,:),[],{'mean' 'sem'}); %mean across columns/timepoints
    [h1,~] = boundedline(timescaleCuePlot, mean1, sem1, 'cmap', persistentCmap, 'alpha','y', 'transparency', 0.5, 'LineWidth', 3);
    [mean2, sem2] = grpstats(cueDataSorted(nonSpecificAntiAvoidIndeces,:),[],{'mean' 'sem'}); 
    [h2,~] = boundedline(timescaleCuePlot, mean2, sem2, 'cmap', nonPersistentCmap, 'alpha','y', 'transparency', 0.5, 'LineWidth', 3);
    title([label ' @ avoid cue']);
    zeroXDottedLine;
    makepretty;
    hold off;

    nexttile; hold on;
    [mean1, sem1] = grpstats(coldDataSorted(nonSpecificAvoidIndeces,:),[],{'mean' 'sem'}); %mean across columns/timepoints
    [h1,~] = boundedline(timescaleColdPlot, mean1, sem1, 'cmap', persistentCmap, 'alpha','y', 'transparency', 0.5, 'LineWidth', 3);
    [mean2, sem2] = grpstats(coldDataSorted(nonSpecificAntiAvoidIndeces,:),[],{'mean' 'sem'}); 
    [h2,~] = boundedline(timescaleColdPlot, mean2, sem2, 'cmap', nonPersistentCmap, 'alpha','y', 'transparency', 0.5, 'LineWidth', 3);
    title('@ cold air');
    zeroXDottedLine;
    makepretty;
    hold off;

    nexttile; hold on;
    [mean1, sem1] = grpstats(postReactDataSorted(nonSpecificAvoidIndeces,:),[],{'mean' 'sem'}); %mean across columns/timepoints
    [h1,~] = boundedline(timescaleVhexPlot, mean1, sem1, 'cmap', persistentCmap, 'alpha','y', 'transparency', 0.5, 'LineWidth', 3);
    [mean2, sem2] = grpstats(postReactDataSorted(nonSpecificAntiAvoidIndeces,:),[],{'mean' 'sem'}); 
    [h2,~] = boundedline(timescaleVhexPlot, mean2, sem2, 'cmap', nonPersistentCmap, 'alpha','y', 'transparency', 0.5, 'LineWidth', 3);
    title('@ react VHEx');
    zeroXDottedLine;
    makepretty;
    hold off;

    nexttile; hold on;
    [mean1, sem1] = grpstats(postAvoidDataSorted(nonSpecificAvoidIndeces,:),[],{'mean' 'sem'}); %mean across columns/timepoints
    [h1,~] = boundedline(timescaleVhexPlot, mean1, sem1, 'cmap', persistentCmap, 'alpha','y', 'transparency', 0.5, 'LineWidth', 3);
    [mean2, sem2] = grpstats(postAvoidDataSorted(nonSpecificAntiAvoidIndeces,:),[],{'mean' 'sem'}); 
    [h2,~] = boundedline(timescaleVhexPlot, mean2, sem2, 'cmap', nonPersistentCmap, 'alpha','y', 'transparency', 0.5, 'LineWidth', 3);
    title('@ avoid VHEx');
    zeroXDottedLine;
    makepretty;
    hold off;


end

function [specificAvoidIndeces, nonSpecificAvoidIndeces] = findAvoidAtCueNeurons(anatData, timescaleSeg, timescaleSegAll, label, thisCmap)
    %find "avoid neurons" that are more active before avoid cue
    %inputs are rpt x units x timeSeg 
    % specificAvoidIndeces are neurons that are more active before avoid VHEx AND relatively inactive before react (perhaps specifically recruited by cue)
    % nonSpecificAvoidIndeces are just neurons that are more active before avoid VHEx 

    avoidWindow = [0.0 0.5]; %after cue, use window where persistent activity separates and also include motor
    avoidSamps = [find(avoidWindow(1)>=timescaleSeg,1,'last') find(avoidWindow(2)>=timescaleSeg,1,'last')];
    nonPersistentCmap = [0 0 0];
    persistentCmap = [0 0.7 0.5];
    blueCmap = [0 0.4470 0.7410];
    orangeCmap = [0.8500 0.3250 0.0980];	
    medGreyCmap = [0.4 0.4 0.4];

    removeNull = 0;
    if removeNull 
%         avoidData = anatData.allUnitsTrialsAvoidNullMean; %option to remove null everywhere
%         reactData = anatData.allUnitsTrialsReactNullMean;
        avoidData = anatData.allUnitsTrialsAvoidMean; %option not to remove cold air null at VHEx here, since we only care about calculation of avoid/peristent neurons at cue
        reactData = anatData.allUnitsTrialsReactMean;
        avoidAtCueData = anatData.allUnitsTrialsAvoidAtCueNullMean;
        reactAtCueData = anatData.allUnitsTrialsReactAtCueNullMean;
    else
        avoidData = anatData.allUnitsTrialsAvoidMean;
        reactData = anatData.allUnitsTrialsReactMean;
        avoidAtCueData = anatData.allUnitsTrialsAvoidAtCueMean;
        reactAtCueData = anatData.allUnitsTrialsReactAtCueMean;
    end
    
    numMiceProbes = size(avoidData,2);
    numRepeats = 10;
    allMiceAvoidData = [];
    allMiceReactData = [];
    allMiceAvoidAtCueData = []; 
    allMiceReactAtCueData = [];
    %first get means over repeats & mice & combine into single neurons x time mean ZFR matrix:
    for i = 1:numMiceProbes
        thisMouseAvoidData = avoidData{i}; %unpack to get per-mouse cells that include all sample repeats
        thisMouseReactData = reactData{i};
        thisMouseAvoidAtCueData = avoidAtCueData{i}; 
        thisMouseReactAtCueData = reactAtCueData{i};

        numUnitsThisMouse = size(thisMouseAvoidData{1},1);
        numTimepoints = size(thisMouseAvoidData{1},2);
        thisMouseAvoidConcat = zeros(numRepeats,numUnitsThisMouse,numTimepoints);
        thisMouseReactConcat = zeros(numRepeats,numUnitsThisMouse,numTimepoints);
        thisMouseAvoidAtCueConcat = zeros(numRepeats,numUnitsThisMouse,numTimepoints);
        thisMouseReactAtCueConcat = zeros(numRepeats,numUnitsThisMouse,numTimepoints);

        for rpt = 1:numRepeats
            thisMouseAvoidConcat(rpt,:,:) = thisMouseAvoidData{rpt};
            thisMouseReactConcat(rpt,:,:) = thisMouseReactData{rpt};
            thisMouseAvoidAtCueConcat(rpt,:,:) = thisMouseAvoidAtCueData{rpt};
            thisMouseReactAtCueConcat(rpt,:,:) = thisMouseReactAtCueData{rpt};
        end
        %mean over repeats & concatenate over mice:
        allMiceAvoidData = [allMiceAvoidData; squeeze(mean(thisMouseAvoidConcat,1))];
        allMiceReactData = [allMiceReactData; squeeze(mean(thisMouseReactConcat,1))];
        allMiceAvoidAtCueData = [allMiceAvoidAtCueData; squeeze(mean(thisMouseAvoidAtCueConcat,1))]; 
        allMiceReactAtCueData = [allMiceReactAtCueData; squeeze(mean(thisMouseReactAtCueConcat,1))];
    end
    
    avoidAtCueMean = allMiceAvoidAtCueData(:,avoidSamps(1):avoidSamps(2));
    reactAtCueMean = allMiceReactAtCueData(:,avoidSamps(1):avoidSamps(2));
    diffAtCueMean = avoidAtCueMean - reactAtCueMean;
    diffAtCueMeanOverWindow = mean(diffAtCueMean,2); %this is 'tuning' to avoid @ cue

    [~,sortIdx] = sort(diffAtCueMeanOverWindow, 'descend'); %sort by avoid-react within window
    avoidMeanSorted = avoidAtCueMean(sortIdx,:);
    reactMeanSorted = reactAtCueMean(sortIdx,:);
    diffMeanSorted = diffAtCueMean(sortIdx,:);

    numUnits = size(avoidAtCueMean,1);
    specificAvoidIndeces = false(numUnits,1); %logical indeces; here "specific" are undefined so leave false
    nonSpecificAvoidIndeces = false(numUnits,1);
    nonSpecificAntiAvoidIndeces = false(numUnits,1);
    for i = 1:numUnits
        if (mean(diffMeanSorted(i,:)) > 0.1)
            nonSpecificAvoidIndeces(i) = true;
        elseif (mean(diffMeanSorted(i,:)) < -0.1)
            nonSpecificAntiAvoidIndeces(i) = true;
        end
    end
%     percentSpecificAvoidNeurons = length(find(specificAvoidIndeces)) / numUnits;
%     percentNonSpecificAvoidNeurons = length(find(nonSpecificAvoidIndeces)) / numUnits;
%     percentNonSpecificAntiAvoidNeurons = length(find(nonSpecificAntiAvoidIndeces)) / numUnits;
    
%     lastAvoidNeuron = find(nonSpecificAvoidIndeces,1,'last');
%     lastAntiAvoidNeuron = find(nonSpecificAntiAvoidIndeces,1,'first');
%     figure; tiledlayout(1,3, "TileSpacing","tight");
%     nexttile; hold on; imagesc(avoidMeanSorted); clim([-1 1]); xLim = get(gca, 'xlim'); plot(xLim, [lastAvoidNeuron lastAvoidNeuron],'k--', 'LineWidth',2); plot(xLim, [lastAntiAvoidNeuron lastAntiAvoidNeuron],'r--', 'LineWidth',2);
%     title(label); axis tight;
%     nexttile; hold on; imagesc(reactMeanSorted); clim([-1 1]); plot(xLim, [lastAvoidNeuron lastAvoidNeuron],'k--', 'LineWidth',2); plot(xLim, [lastAntiAvoidNeuron lastAntiAvoidNeuron],'r--', 'LineWidth',2); axis tight;
%     nexttile; hold on; imagesc(diffMeanSorted); clim([-1 1]); plot(xLim, [lastAvoidNeuron lastAvoidNeuron],'k--', 'LineWidth',2); plot(xLim, [lastAntiAvoidNeuron lastAntiAvoidNeuron],'r--', 'LineWidth',2);
%     title(['% avoid neurons = ' num2str(percentNonSpecificAvoidNeurons*100)]); axis tight;
%     keyboard

    %now compare mean 'tuning' to other relevant events:
    cueWindow = [-0.5 2];
    cueSamps = [find(cueWindow(1)>=timescaleSeg,1,'last') find(cueWindow(2)>=timescaleSeg,1,'last')];
    timescaleCuePlot = timescaleSeg(cueSamps(1):cueSamps(2));
    cueData = allMiceAvoidAtCueData(:,cueSamps(1):cueSamps(2)); %mean over repeats & window
    cueDataSorted = cueData(sortIdx,:); %use same sort as above so that xIndeces apply
%     figure; imagesc(cueDataSorted); clim([-1 1]);

    coldAirWindow = [-1 1];
    coldSamps = [find(coldAirWindow(1)>=timescaleSegAll,1,'last') find(coldAirWindow(2)>=timescaleSegAll,1,'last')];
    timescaleColdPlot = timescaleSegAll(coldSamps(1):coldSamps(2));
    coldData = anatData.allUnitsColdAirOnMeanZFR(:,coldSamps(1):coldSamps(2)); %no repeats for cold air on data
    coldDataSorted = coldData(sortIdx,:); %use same sort as above so that xIndeces apply

    vhexWindow = [-2 2];
    vhexSamps = [find(vhexWindow(1)>=timescaleSeg,1,'last') find(vhexWindow(2)>=timescaleSeg,1,'last')];
    timescaleVhexPlot = timescaleSeg(vhexSamps(1):vhexSamps(2));
    postReactData = allMiceReactData(:,vhexSamps(1):vhexSamps(2)); %mean over repeats & window
    postReactDataSorted = postReactData(sortIdx,:); %use same sort as above so that xIndeces apply

    postAvoidData = allMiceAvoidData(:,vhexSamps(1):vhexSamps(2)); %mean over repeats & window
    postAvoidDataSorted = postAvoidData(sortIdx,:); %use same sort as above so that xIndeces apply

    % mean activity plots across population %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     yLimAll = [-0.3 0.7];
%     f = figure;
%     tiledlayout(1,3,'TileSpacing', 'tight');
%     nexttile; hold on;
%     [mean1, sem1] = grpstats(cueDataSorted(nonSpecificAvoidIndeces,:),[],{'mean' 'sem'}); %mean across columns/timepoints
%     [h1,~] = boundedline(timescaleCuePlot, mean1, sem1, 'cmap', persistentCmap, 'alpha','y', 'transparency', 0.5, 'LineWidth', 3);
%     [mean2, sem2] = grpstats(cueDataSorted(nonSpecificAntiAvoidIndeces,:),[],{'mean' 'sem'}); 
%     [h2,~] = boundedline(timescaleCuePlot, mean2, sem2, 'cmap', nonPersistentCmap, 'alpha','y', 'transparency', 0.5, 'LineWidth', 3);
%     title([label ' @ avoid cue']);
%     zeroXDottedLine;
%     ylim(yLimAll)
%     makepretty;
%     hold off;
% 
% %     nexttile; hold on;
% %     [mean1, sem1] = grpstats(coldDataSorted(nonSpecificAvoidIndeces,:),[],{'mean' 'sem'}); %mean across columns/timepoints
% %     [h1,~] = boundedline(timescaleColdPlot, mean1, sem1, 'cmap', persistentCmap, 'alpha','y', 'transparency', 0.5, 'LineWidth', 3);
% %     [mean2, sem2] = grpstats(coldDataSorted(nonSpecificAntiAvoidIndeces,:),[],{'mean' 'sem'}); 
% %     [h2,~] = boundedline(timescaleColdPlot, mean2, sem2, 'cmap', nonPersistentCmap, 'alpha','y', 'transparency', 0.5, 'LineWidth', 3);
% %     title('@ cold air');
% %     zeroXDottedLine;
% %     ylim(yLimAll)
% %     makepretty;
% %     hold off;
% 
%     nexttile; hold on;
%     [mean1, sem1] = grpstats(postReactDataSorted(nonSpecificAvoidIndeces,:),[],{'mean' 'sem'}); %mean across columns/timepoints
%     [h1,~] = boundedline(timescaleVhexPlot, mean1, sem1, 'cmap', persistentCmap, 'alpha','y', 'transparency', 0.5, 'LineWidth', 3);
%     [mean2, sem2] = grpstats(postReactDataSorted(nonSpecificAntiAvoidIndeces,:),[],{'mean' 'sem'}); 
%     [h2,~] = boundedline(timescaleVhexPlot, mean2, sem2, 'cmap', nonPersistentCmap, 'alpha','y', 'transparency', 0.5, 'LineWidth', 3);
%     title('@ react VHEx');
%     zeroXDottedLine;
%     ylim(yLimAll)
%     makepretty;
%     hold off;
% 
%     nexttile; hold on;
%     [mean1, sem1] = grpstats(postAvoidDataSorted(nonSpecificAvoidIndeces,:),[],{'mean' 'sem'}); %mean across columns/timepoints
%     [h1,~] = boundedline(timescaleVhexPlot, mean1, sem1, 'cmap', persistentCmap, 'alpha','y', 'transparency', 0.5, 'LineWidth', 3);
%     [mean2, sem2] = grpstats(postAvoidDataSorted(nonSpecificAntiAvoidIndeces,:),[],{'mean' 'sem'}); 
%     [h2,~] = boundedline(timescaleVhexPlot, mean2, sem2, 'cmap', nonPersistentCmap, 'alpha','y', 'transparency', 0.5, 'LineWidth', 3);
%     title('@ avoid VHEx');
%     zeroXDottedLine;
%     ylim(yLimAll)
%     makepretty;
%     hold off;
%     f.Position(3:4) = [1500 300]; %set figure size appropriate for plots

    %scatterplots for 'tuning' comparisons %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    coldAirWindow = [0 0.1]; 
    coldSamps = [find(coldAirWindow(1)>=timescaleSegAll,1,'last') find(coldAirWindow(2)>=timescaleSegAll,1,'last')];
    coldData = anatData.allUnitsColdAirOnMeanZFR(:,coldSamps(1):coldSamps(2)); %no repeats for cold air on data
    itiData = anatData.itiMeanZFR(:,coldSamps(1):coldSamps(2));
    coldDataDiff = coldData - itiData; %subtract from mid-ITI data (of same window length) to match subtraction procedure for avoid-react
    coldDataMean = mean(coldDataDiff,2); 

    vhexWindow = [-0.5 0];
    vhexSamps = [find(vhexWindow(1)>=timescaleSeg,1,'last') find(vhexWindow(2)>=timescaleSeg,1,'last')];
    itiSamps = [find(vhexWindow(1)>=timescaleSegAll,1,'last') find(vhexWindow(2)>=timescaleSegAll,1,'last')];
    reactData = allMiceReactData(:,vhexSamps(1):vhexSamps(2)); %mean over repeats & window
    itiData = anatData.itiMeanZFR(:,itiSamps(1):itiSamps(2));
    reactDataDiff = reactData - itiData; 
    reactVhexMean = mean(reactDataDiff,2);

    avoidData = allMiceAvoidData(:,vhexSamps(1):vhexSamps(2)); %mean over repeats & window
    avoidDataDiff = avoidData - itiData; 
    avoidVhexMean = mean(avoidDataDiff,2);

    axisLim = [-3 3 -3 3];
    f2 = figure;
    tiledlayout(1,3,'TileSpacing', 'tight');
    nexttile; hold on;
    scatter(diffAtCueMeanOverWindow, reactVhexMean, 60, thisCmap,'o')
    xlabel('avoid-react @ cue')
    ylabel('react-ITI @ VHEx')
    axis(axisLim)
    zeroXDottedLine;
    zeroYDottedLine;
    makepretty;
    title(label)

    nexttile; hold on;
    scatter(diffAtCueMeanOverWindow, avoidVhexMean, 60, thisCmap,'o')
    xlabel('avoid-react @ cue')
    ylabel('avoid-ITI @ VHEx')
    axis(axisLim)
    zeroXDottedLine;
    zeroYDottedLine;
    makepretty;

    nexttile; hold on;
    scatter(diffAtCueMeanOverWindow, coldDataMean, 60, thisCmap,'o')
    xlabel('avoid-react @ cue')
    ylabel('coldAir-ITI @ VHEx')
    axis(axisLim)
    zeroXDottedLine;
    zeroYDottedLine;
    makepretty;
    hold off;

    f2.Position(3:4) = [1200 300]; %set figure size appropriate for plots

    %pie charts for 'tuning' comparisons %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %first find [avoid-react @ cue] tuned neurons above threshold
    avoidCueIndeces = false(numUnits,1); %logical indeces; here "specific" are undefined so leave false
    alsoReactVhexIndeces = false(numUnits,1);
    alsoAvoidVhexIndeces = false(numUnits,1);
    for i = 1:numUnits
        if (diffAtCueMeanOverWindow(i) > 0)
            avoidCueIndeces(i) = true;
            %now for these neurons, see if additionally tuned to other conditions above threshold:
            if (reactVhexMean(i) > 0)
                alsoReactVhexIndeces(i) = true;
            end
            if (avoidVhexMean(i) > 0)
                alsoAvoidVhexIndeces(i) = true;
            end
        end
    end

    percentAvoidCueTuned = length(find(avoidCueIndeces)) / numUnits;
    percentNotAvoidCueTuned = 1 - percentAvoidCueTuned;
    percentAlsoTunedReact = length(find(alsoReactVhexIndeces)) / numUnits;
    percentNotAlsoTunedReact = percentAvoidCueTuned - percentAlsoTunedReact;
    percentAlsoTunedAvoid = length(find(alsoAvoidVhexIndeces)) / numUnits;
    percentNotAlsoTunedAvoid = percentAvoidCueTuned - percentAlsoTunedAvoid;

    f3 = figure;
    t = tiledlayout(1,2,TileSpacing="loose");
    nexttile; hold on;
    explode = [0 1 1];
    p = pie([percentNotAvoidCueTuned percentAlsoTunedReact percentNotAlsoTunedReact],explode);
    title(label)
    pText = findobj(p,'Type','text');
    percentValues = get(pText,'String'); 
    txt = {'not tuned avoid cue: ';'also tuned react VHEx: ';'not also tuned react VHEx: '}; 
    combinedtxt = strcat(txt,percentValues); 
    pText(1).String = combinedtxt(1);
    pText(2).String = combinedtxt(2);
    pText(3).String = combinedtxt(3);
    pp = findobj(p,'Type','patch');
    pp(1).FaceColor = [0 0 0];
    pp(2).FaceColor = blueCmap;
    pp(3).FaceColor = medGreyCmap;
    axis off

    nexttile; hold on;
    explode = [0 1 1];
    p = pie([percentNotAvoidCueTuned percentAlsoTunedAvoid percentNotAlsoTunedAvoid],explode);
    pText = findobj(p,'Type','text');
    percentValues = get(pText,'String'); 
    txt = {' ';'also tuned avoid VHEx: ';'not also tuned avoid VHEx: '}; 
    combinedtxt = strcat(txt,percentValues); 
    pText(1).String = combinedtxt(1);
    pText(2).String = combinedtxt(2);
    pText(3).String = combinedtxt(3);
    pp = findobj(p,'Type','patch');
    pp(1).FaceColor = [0 0 0];
    pp(2).FaceColor = orangeCmap;
    pp(3).FaceColor = medGreyCmap;
    axis off

    set(gcf,'color','w');
    f3.Position(3:4) = [1200 300]; %set figure size appropriate for plots

end

function plotEphysMeanZscoreFREvent(allUnitsMeanZFR, timescaleSeg, titleStr) 
%plot output from above
    figure;
    hold on;
    imagesc(timescaleSeg, 1:1:size(allUnitsMeanZFR,1), allUnitsMeanZFR)
    set(gca, 'YDir', 'normal');
    axis tight;
    x = [0 0]; %dottend line at event = time 0
    yLims = get(gca, 'YLim');
    y = [yLims(1) yLims(2)];
    plot(x, y, 'k--', 'LineWidth', 0.1);

    colormap(blueOrange)
%     colormap(parula)
    clim([-1 1]);
    cb = colorbar;
%     cb.Position = cb.Position + [0.09 0 0 0];
    cb.Label.String = ['z-score ' '\Delta' ' firing rate (Hz)'];
    ylabel('unit #');
    xlabel('time (sec)');
    title(['PSTH ' titleStr]);
    hold off;
    makepretty;
    set(gcf,'color','w'); %set figure background white
end

function plotEphysZscoreFRConcatEvents(allUnitsZFR, numEvents, timescaleSeg, sortIdx, secBefore, secAfter) 
%plot concatenated trials across entire session
    smoothWindow = 5; % how many PSTH bins in window for Gaussian smoothing below; ex. 20ms bins x 5 bins = 100ms
    if sortIdx ~= 0
        dataToPlot = smoothdata(allUnitsZFR(sortIdx,:),2,'gaussian',smoothWindow); %smooth firing rates
    else
        dataToPlot = smoothdata(allUnitsZFR,2,'gaussian',smoothWindow); %default depth order
    end
    
    fullTime = (secBefore + secAfter)*numEvents;
    timeStep = timescaleSeg(2) - timescaleSeg(1);
    timescaleFull = 0:timeStep:fullTime;
    sampsPerSeg = length(timescaleSeg);
    numSampAfter = round(secAfter/timeStep);
    numSampBefore = round(secBefore/timeStep);
    zeroIdx = [numSampBefore:sampsPerSeg:sampsPerSeg*numEvents];
    startIdx = [1:sampsPerSeg:sampsPerSeg*numEvents];
    endIdx = [sampsPerSeg:sampsPerSeg:sampsPerSeg*numEvents];

    f = figure;
    hold on;
    imagesc(timescaleFull, 1:1:size(dataToPlot,1), dataToPlot)
    set(gca, 'YDir', 'reverse');
    axis tight;

    %for dotted lines at extensions & solid lines between segments:
    jumpIndecesDottedLinesX = nan(1,numEvents*3);
    jumpIndecesDottedLinesX(1:3:end) = timescaleFull(zeroIdx);
    jumpIndecesDottedLinesX(2:3:end) = timescaleFull(zeroIdx);
    jumpIndecesDottedLinesY = nan(1,numEvents*3);
    jumpIndecesDottedLinesY(1:3:end) = 0;
    jumpIndecesDottedLinesY(2:3:end) = ones(1,numEvents).*size(dataToPlot,1);
    plot(jumpIndecesDottedLinesX,jumpIndecesDottedLinesY, 'k--', 'LineWidth', 1);
    segmentLinesX = nan(1,numEvents*3);
    segmentLinesX(1:3:end) = timescaleFull(startIdx);
    segmentLinesX(2:3:end) = timescaleFull(startIdx);
    segmentLinesY = nan(1,numEvents*3);
    segmentLinesY(1:3:end) = 0;
    segmentLinesY(2:3:end) = ones(1,numEvents).*size(dataToPlot,1);
    plot(segmentLinesX,segmentLinesY, 'k-', 'LineWidth', 2);
    segmentLinesX(1:3:end) = timescaleFull(endIdx);
    segmentLinesX(2:3:end) = timescaleFull(endIdx);
    segmentLinesY(1:3:end) = 0;
    segmentLinesY(2:3:end) = ones(1,numEvents).*size(dataToPlot,1);
    plot(segmentLinesX,segmentLinesY, 'k-', 'LineWidth', 2);

    colormap(blueOrange)
    clim([-1 1]);
    cb = colorbar;
    cb.Label.String = ['z-score ' '\Delta' ' firing rate (Hz)'];
    ylabel('unit #');
    xlabel('time (sec)');
    hold off;
    f.Position(3:4) = [2000 500]; %set figure size appropriate for plots
    movegui(f,'center')
    makepretty;
    set(gcf,'color','w'); %set figure background white

end

function plotEphysRasterConcatEvents(shankData, tsData, idxData, eventTimes, sortIdx, secBefore, secAfter) 
%plot raster of concatenated trials across entire session w/ trial type marks
    numEvents = length(eventTimes);    
    trialSubset = [1:15]; %optionally choose a smaller (~10-20) subset of trials to see spiking details more clearly across different trial types
    %trialSubset = [1:numEvents];
    timescale = tsData.timescale; %behavioral timescale
    segmentTime = secBefore + secAfter;
    fullTime = segmentTime*numEvents;
    timeStep = timescale(2) - timescale(1);
    numSampAfter = round(secAfter/timeStep);
    numSampBefore = round(secBefore/timeStep);
    segmentTimescale = -secBefore:timeStep:(secAfter);
    stCell = shankData.stCell;
    stSampRate = shankData.sampRate;
    stTimeStep = 1/stSampRate;
    stTimescale = 0:(1/stSampRate):shankData.totalRecordingTime; %probe timescale

    avoidIdx = idxData.jumpAntIndeces;
    % reactIdx = ~idxData.jumpAntIndeces; %can assume react if not opto or avoid
    optoIdx = idxData.jumpOptoIndeces;
    blueCmap = [0 0.4470 0.7410];
    orangeCmap = [0.8500 0.3250 0.0980];

    if sortIdx ~= 0 
        rasterIdx = sortIdx; %sets row order for plotting below
    else
        rasterIdx = 1:1:length(stCell);
    end

    f = figure; tiledlayout(2,1,"TileSpacing","tight");
    currentStartTime = 0; %x-axis plot time
    for i = trialSubset
    % for each event, gather all spike times within segment/window around
        currentEventTime = eventTimes(i);
        currentEndTime = currentStartTime + segmentTime;
    
        nexttile(1); hold on;
        currentUnitIdx = 1;
        for k = rasterIdx %now go through each row,unit, rasterize & plot at each event
            rasterX = stCell{k}; %spike times
            firstSpike = find(rasterX>(currentEventTime-secBefore),1,'first');
            lastSpike = find(rasterX<(currentEventTime+secAfter),1,'last');
            rasterXSeg = rasterX(firstSpike:lastSpike) - currentEventTime + currentStartTime + secBefore; %get segment & shift times to current segment on plot
            %classic slash raster:
            rasterXslash = nan(1,length(rasterXSeg)*3); %use 3rd NaN column so plot is discontinuous
            rasterXslash(1:3:end) = rasterXSeg; %every first & second point is the spike time, but third point is NaN
            rasterXslash(2:3:end) = rasterXSeg;
            rasterY = nan(1,length(rasterXSeg)*3);
            minValues = ones(1,length(rasterXSeg)).*currentUnitIdx; %for each unit, the raster ticks are 1 a.u.
            rasterY(1:3:end) = minValues;
            rasterY(2:3:end) = minValues+1;
            plot(rasterXslash,rasterY,'Color',[0 0 0 0.5], 'Marker', 'none', 'LineWidth', 0.5); %4 element color vector for transparency
            %plot(rasterXslash,rasterY,'Color',[0 0 0], 'Marker', 'none', 'LineWidth', 0.1); %no transparency option for easier figure import
            currentUnitIdx = currentUnitIdx + 1;
        end
    
        % for current segment plot dotted line at VHEx colored by trial type & solid lines separator after segment:
        if optoIdx(i)
            dottedCmap = [0 0 0]; %black for opto trial
        elseif avoidIdx(i)
            dottedCmap = orangeCmap;
        else
            dottedCmap = blueCmap;
        end

        dottedLineX(1) = currentStartTime+secBefore;
        dottedLineX(2) = currentStartTime+secBefore;
        dottedLineY(1) = 0;
        dottedLineY(2) = size(stCell,2);
        plot(dottedLineX,dottedLineY, 'LineStyle','--', 'LineWidth', 3, 'Color', dottedCmap);

        segmentLineX(1) = currentStartTime+segmentTime;
        segmentLineX(2) = currentStartTime+segmentTime;
        plot(segmentLineX,dottedLineY, 'g-', 'LineWidth', 2);



        % now plot platform distance segments underneath
        nexttile(2); hold on;
        currentEventIdx = find(timescale>currentEventTime,1,'first');
        platSeg = tsData.jumpDistance(currentEventIdx-numSampBefore:currentEventIdx+numSampAfter-1) - median(tsData.jumpDistance);
        platTimescale = segmentTimescale+currentStartTime+secBefore;
        plot(platTimescale,platSeg, 'm-', 'LineWidth', 2);
        dottedLineY(2) = 5; %~max platform distance
        plot(dottedLineX,dottedLineY, 'LineStyle','--', 'LineWidth', 3, 'Color', dottedCmap);
        plot(segmentLineX,dottedLineY, 'g-', 'LineWidth', 2);
        
        currentStartTime = currentEndTime + timeStep;
    end

    % finish figure
    nexttile(1);
    axis tight
    ylabel('neurons');
    set(gca, 'XTick', []);
    hold off;
    makepretty;
    nexttile(2);
    axis tight
    %pbaspect([timestep/stTimeStep, 1, 1]); %set x axes to same aspect ratio
    ylabel({'platform'; 'distance (cm)'});
    set(gca, 'XTick', []);
    xlabel('time (sec)');
    hold off;
    makepretty;
    f.Position(3:4) = [2000 500]; %set figure size appropriate for plots
    movegui(f,'center')
    set(gcf,'color','w'); %set figure background white

end

function plotEphysMeanFREvent(allUnitsMeanFR, timescaleSeg, titleStr) 
%plot output from above, NOT zscored
    figure;
    hold on;
    imagesc(timescaleSeg, 1:1:size(allUnitsMeanFR,1), allUnitsMeanFR)
    set(gca, 'YDir', 'normal');
    axis tight;
    zeroXDottedLine;

    colormap(blueOrange)
    clim([-0.5 0.5]);
    cb = colorbar;
    cb.Label.String = ['z-score FR'];
    ylabel('unit #');
    xlabel('time (sec)');
    title(['PSTH ' titleStr]);
    hold off;
    makepretty;
    set(gcf,'color','w'); %set figure background white
end

function plotEphysMeanZscoreFR2Events(allUnitsMeanZFR1, allUnitsMeanZFR2, timescaleSeg) 
% plot 2 diff mean PSTHs at 2 events, side by side
    figure;
    tiledlayout(1,2,'TileSpacing', 'tight');
    climsPsth = [-0.5 0.5];
    fontSz = 16;

    nexttile
    hold on;
    imagesc(timescaleSeg, 1:1:size(allUnitsMeanZFR1,1), allUnitsMeanZFR1)
    colormap(blueOrange)
    set(gca, 'YDir', 'reverse');
    axis tight;
    zeroXDottedLine;
    clim(climsPsth);
    ylabel('unit #', 'FontSize', fontSz);
    xlabel('time (sec)', 'FontSize', fontSz);
    title('@ cue', 'FontSize', fontSz);
    hold off;

    nexttile
    hold on;
    imagesc(timescaleSeg, 1:1:size(allUnitsMeanZFR2,1), allUnitsMeanZFR2)
    colormap(blueOrange)
    set(gca, 'YDir', 'reverse');
    axis tight;
    zeroXDottedLine;
    clim(climsPsth);
    set(gca, 'YTick', []);
    cb = colorbar;
    cb.Label.String = ['z-score ' '\Delta' ' firing rate (Hz)'];
    cb.Label.FontSize = fontSz;
    xlabel('time (sec)', 'FontSize', fontSz);
    title('@ VHEx', 'FontSize', fontSz);
    hold off;

    set(gcf,'color','w'); %set figure background white

end

function plotEphysMeanZscoreFR3Events(allUnitsMeanZFR1, allUnitsMeanZFR2, allUnitsMeanZFR3, timescaleSeg) 
% plot 2 diff mean PSTHs at 2 events, side by side
    figure;
    tiledlayout(1,3,'TileSpacing', 'tight');
    climsPsth = [-0.5 0.5];
    fontSz = 16;

    nexttile
    hold on;
    imagesc(timescaleSeg, 1:1:size(allUnitsMeanZFR1,1), allUnitsMeanZFR1)
    colormap(blueOrange)
    set(gca, 'YDir', 'reverse');
    axis tight;
    zeroXDottedLine;
    clim(climsPsth);
    ylabel('unit #', 'FontSize', fontSz);
    xlabel('time (sec)', 'FontSize', fontSz);
    title('@ cue', 'FontSize', fontSz);
    hold off;

    nexttile
    hold on;
    imagesc(timescaleSeg, 1:1:size(allUnitsMeanZFR2,1), allUnitsMeanZFR2)
    colormap(blueOrange)
    set(gca, 'YDir', 'reverse');
    axis tight;
    zeroXDottedLine;
    clim(climsPsth);
    set(gca, 'YTick', []);
    xlabel('time (sec)', 'FontSize', fontSz);
    title('@ VHEx', 'FontSize', fontSz);
    hold off;

    nexttile
    hold on;
    imagesc(timescaleSeg, 1:1:size(allUnitsMeanZFR3,1), allUnitsMeanZFR3)
    colormap(blueOrange)
    set(gca, 'YDir', 'reverse');
    axis tight;
    zeroXDottedLine;
    clim(climsPsth);
    set(gca, 'YTick', []);
    cb = colorbar;
    cb.Label.String = ['z-score ' '\Delta' ' firing rate (Hz)'];
    cb.Label.FontSize = fontSz;
    xlabel('time (sec)', 'FontSize', fontSz);
    title('@ ITI', 'FontSize', fontSz);
    hold off;

    set(gcf,'color','w'); %set figure background white

end

function plotEphysOptoTestTrials(optoOnMeanZFR_dist1, optoOnMeanZFR_dist2, optoOnMeanZFR_dist3, optoOnMeanZFR_dist4, optoOnMeanZFR_dist5, optoOnMeanZFR_dist6, timescaleSeg)
% plot mean PSTHs for optoTest at different probe-opto distances, side by side
    figure;
    tiledlayout(1,6,'TileSpacing', 'tight');
    climsPsth = [-1 1];
    fontSz = 16;

    nexttile
    hold on;
    imagesc(timescaleSeg, 1:1:size(optoOnMeanZFR_dist1,1), optoOnMeanZFR_dist1)
    colormap(blueOrange)
    set(gca, 'YDir', 'normal');
    axis tight;
    zeroXDottedLine;
    clim(climsPsth);
    ylabel('unit #', 'FontSize', fontSz);
    xlabel('time (sec)', 'FontSize', fontSz);
    title('0.0 mm', 'FontSize', fontSz);
    hold off;

    nexttile
    hold on;
    imagesc(timescaleSeg, 1:1:size(optoOnMeanZFR_dist2,1), optoOnMeanZFR_dist2)
    colormap(blueOrange)
    set(gca, 'YDir', 'normal');
    axis tight;
    zeroXDottedLine;
    clim(climsPsth);
    xlabel('time (sec)', 'FontSize', fontSz);
    set(gca, 'YTick', []);
    title('0.5 mm', 'FontSize', fontSz);
    hold off;

    nexttile
    hold on;
    imagesc(timescaleSeg, 1:1:size(optoOnMeanZFR_dist3,1), optoOnMeanZFR_dist3)
    colormap(blueOrange)
    set(gca, 'YDir', 'normal');
    axis tight;
    zeroXDottedLine;
    clim(climsPsth);
    xlabel('time (sec)', 'FontSize', fontSz);
    set(gca, 'YTick', []);
    title('0.75 mm', 'FontSize', fontSz);
    hold off;

    nexttile
    hold on;
    imagesc(timescaleSeg, 1:1:size(optoOnMeanZFR_dist4,1), optoOnMeanZFR_dist4)
    colormap(blueOrange)
    set(gca, 'YDir', 'normal');
    axis tight;
    zeroXDottedLine;
    clim(climsPsth);
    xlabel('time (sec)', 'FontSize', fontSz);
    set(gca, 'YTick', []);
    title('1.0 mm', 'FontSize', fontSz);
    hold off;

    nexttile
    hold on;
    imagesc(timescaleSeg, 1:1:size(optoOnMeanZFR_dist5,1), optoOnMeanZFR_dist5)
    colormap(blueOrange)
    set(gca, 'YDir', 'normal');
    axis tight;
    zeroXDottedLine;
    clim(climsPsth);
    xlabel('time (sec)', 'FontSize', fontSz);
    set(gca, 'YTick', []);
    title('2.0 mm', 'FontSize', fontSz);
    hold off;

    nexttile
    hold on;
    imagesc(timescaleSeg, 1:1:size(optoOnMeanZFR_dist6,1), optoOnMeanZFR_dist6)
    colormap(blueOrange)
    set(gca, 'YDir', 'normal');
    axis tight;
    zeroXDottedLine;
    clim(climsPsth);
    cb = colorbar;
    cb.Label.String = ['z-score ' '\Delta' ' firing rate (Hz)'];
    cb.Label.FontSize = fontSz;
    xlabel('time (sec)', 'FontSize', fontSz);
    set(gca, 'YTick', []);
    title('3.0 mm', 'FontSize', fontSz);
    hold off;

    set(gcf,'color','w'); %set figure background white

end

function plotEphysOptoTestTrialsRaster(tsData, stCell, optoTimes, sortIdx)
% optoTest at different probe-opto distances, side by side: raster w/
% opto overlay
    figure;
    tiledlayout(1,6,'TileSpacing', 'tight');
    fontSz = 16;
    distances = [0.0 0.5 0.75 1.0 2.0 3.0];
    optoDuration = 2.0;
    optoRampDown = 0.5;

    for i = 1:length(distances)

        if sortIdx ~= 0 
            rasterIdx = sortIdx;
        else
            rasterIdx = 1:1:length(stCell); %default sort
        end
    
        optoTimeCurr = optoTimes((i-1)*5 + 3); %there were 5x repeats for each distance, so just pick Xth repeat as example
        secBefore = 1;
        secAfter = 4;
        sampTime = tsData.timescale(2) - tsData.timescale(1);
        sampRate = 1/sampTime;
        numSampAfter = round(sampRate*secAfter);
        numSampBefore = round(sampRate*secBefore);
        timescaleSeg = (-numSampBefore*sampTime):sampTime:(numSampAfter*sampTime - sampTime);
    
        nexttile
        hold on;

        currentIdx = 1;
        for k = rasterIdx %now go through each row,unit, rasterize & plot
            %classic slash raster option:
            rasterX = stCell{k}; %spike times
            firstSpike = find(rasterX>(optoTimeCurr-secBefore),1,'first');
            lastSpike = find(rasterX<(optoTimeCurr+secAfter),1,'last');
            rasterXSeg = rasterX(firstSpike:lastSpike) - optoTimeCurr; %subtract to center 0 at opto on
    
            rasterXslash = nan(1,length(rasterXSeg)*3); %use 3rd NaN column so plot is discontinuous
            rasterXslash(1:3:end) = rasterXSeg; %every first & second point is the spike time, but third point is NaN
            rasterXslash(2:3:end) = rasterXSeg;
            rasterY = nan(1,length(rasterXSeg)*3);
            minValues = ones(1,length(rasterXSeg)).*currentIdx; %for each unit, the raster ticks are 1 a.u.
            rasterY(1:3:end) = minValues;
            rasterY(2:3:end) = minValues+1;
            plot(rasterXslash,rasterY,'Color',[0 0 0 0.5], 'Marker', 'none', 'LineWidth', 0.01); %4 element color vector for transparency
            currentIdx = currentIdx + 1;
        end

    
        %colored rectangle patch when opto on, separate offset patch for ramp
        %down to apply gradient color elsewhere (or just leave off since matlab
        %groups patches on SVG export, & use graphics program to create)
        y = [1 length(stCell)+1 length(stCell)+1 1];
        x = [0 0 optoDuration optoDuration];
        patch(x, y, [0 0.5 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
    %     x = [optoDuration+0.01 optoDuration+0.01 optoDuration+optoRampDown optoDuration+optoRampDown];
    %     patch(x, y, [0 0.5 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
    
        hold off;
        axis tight;
        ylabel('cluster/unit #', 'FontSize', fontSz);
        xlabel('time (sec)', 'FontSize', fontSz);
        title( num2str(distances(i)), ' mm', 'FontSize', fontSz); 
        makepretty;
    end

    set(gcf,'color','w'); %set figure background white
end

function plotEphysOptoTestPercentages(optoOnMeanPercentFRsAll)
% optoTest at different probe-opto distances, side by side: raster w/
% opto overlay
    figure;
    hold on;
    fontSz = 16;
    %distances = {'0.0', '0.5', '0.75', '1.0', '2.0', '3.0'};
    %medGreyCmap = [0.4 0.4 0.4];
    [percentsMean, percentsSem] = grpstats(optoOnMeanPercentFRsAll,[],{'mean' 'sem'}); %mean across columns/timepoints
    %plotSpread(optoOnMeanPercentFRsAll, 'xNames', distances, 'distributionColors', {medGreyCmap medGreyCmap medGreyCmap medGreyCmap medGreyCmap medGreyCmap}); %works but a bit too busy
    boundedline([1:6], percentsMean, percentsSem, 'cmap', [0 0 0])
    hold off;
    axis tight;
    ylabel('% baseline firing rate', 'FontSize', fontSz);
    xlabel('distance from optical fiber to probe (mm)', 'FontSize', fontSz);
    makepretty;
    set(gcf,'color','w'); %set figure background white

    %nonparametric stats:

end

function plotBasicNeuronTuningGrid(allCueMeanZFR, allColdAirZFR, allVhexMeanZFR, avoidMeanZFR, reactMeanZFR, allUnitsAnatDepths, timescaleSeg, plotWindow, itiWindow) 
%plot ant-react order next to other events in grid layout
    plotSamps = [find(plotWindow(1)>=timescaleSeg,1,'last') find(plotWindow(2)>=timescaleSeg,1,'last')];
    itiSamps = [find(itiWindow(1)>=timescaleSeg,1,'last') find(itiWindow(2)>=timescaleSeg,1,'last')];
    timescalePlot = timescaleSeg(plotSamps(1):plotSamps(2));
    midPt = find(0>=timescaleSeg,1,'last'); % for sorting
    
    zplotTitles = {'@ cue', '@ cold air', '@ all VHEx', '@ avoid only', '@ react only'};
    depthSort = 0;
    % use depth sort or rastermap sort for each plot
    if depthSort
        [~,sortIdx] = sort(allUnitsAnatDepths, 'descend'); %'depths' are distances from probe tip, so smaller numbers are deepest in brain
        zplot1 = allCueMeanZFR(sortIdx,plotSamps(1):plotSamps(2));
        zplot2 = allColdAirZFR(sortIdx,plotSamps(1):plotSamps(2));
        zplot3 = allVhexMeanZFR(sortIdx,plotSamps(1):plotSamps(2));
        zplot4 = avoidMeanZFR(sortIdx,plotSamps(1):plotSamps(2));
        zplot5 = reactMeanZFR(sortIdx,plotSamps(1):plotSamps(2));
    else % else can use beginning - end diff for each, or try rastermap
        [~,sortIdx3] = sort( mean(allVhexMeanZFR(:,itiSamps(1):itiSamps(2)),2) - mean(allVhexMeanZFR(:,plotSamps(1):midPt),2), 'descend');

        [~,sortIdx1] = sort( mean(allCueMeanZFR(:,plotSamps(1):midPt),2) - mean(allCueMeanZFR(:,midPt+1:plotSamps(2)),2), 'descend');
        zplot1 = allCueMeanZFR(sortIdx3,plotSamps(1):plotSamps(2));
        [~,sortIdx2] = sort( mean(allColdAirZFR(:,plotSamps(1):midPt),2) - mean(allColdAirZFR(:,midPt+1:plotSamps(2)),2), 'descend');
        zplot2 = allColdAirZFR(sortIdx3,plotSamps(1):plotSamps(2));

        [~,sortIdx3] = sort( mean(allVhexMeanZFR(:,itiSamps(1):itiSamps(2)),2) - mean(allVhexMeanZFR(:,plotSamps(1):midPt),2), 'descend');
        zplot3 = allVhexMeanZFR(sortIdx3,plotSamps(1):plotSamps(2));
%         [~,sortIdx4] = sort(mean(avoidMeanZFR(:,1:midPt),2)-mean(avoidMeanZFR(:,midPt+1:end),2), 'descend');
        zplot4 = avoidMeanZFR(sortIdx3,plotSamps(1):plotSamps(2));
%         [~,sortIdx5] = sort(mean(reactMeanZFR(:,1:midPt),2)-mean(reactMeanZFR(:,midPt+1:end),2), 'descend');
        zplot5 = reactMeanZFR(sortIdx3,plotSamps(1):plotSamps(2));
    end
    

    figure;
    tiledlayout(1,5,'TileSpacing', 'tight');
    climsPsth = [-0.5 0.5];
    fontSz = 16;

    nexttile
    hold on;
    imagesc(timescalePlot, 1:1:size(zplot1,1), zplot1)
    colormap(blueOrange)
    set(gca, 'YDir', 'reverse');
    axis tight;
    zeroXDottedLine;
    clim(climsPsth);
    ylabel('unit #', 'FontSize', fontSz);
    xlabel('time (sec)', 'FontSize', fontSz);
    title([zplotTitles{1}], 'FontSize', fontSz);
    hold off;

    nexttile
    hold on;
    imagesc(timescalePlot, 1:1:size(zplot2,1), zplot2)
    colormap(blueOrange)
    set(gca, 'YDir', 'reverse');
    axis tight;
    zeroXDottedLine;
    clim(climsPsth);
    xlabel('time (sec)', 'FontSize', fontSz);
    set(gca, 'YTick', []);
    title([zplotTitles{2}], 'FontSize', fontSz);
    hold off;

    nexttile
    hold on;
    imagesc(timescalePlot, 1:1:size(zplot3,1), zplot3)
    colormap(blueOrange)
    set(gca, 'YDir', 'reverse');
    axis tight;
    zeroXDottedLine;
    clim(climsPsth);
    xlabel('time (sec)', 'FontSize', fontSz);
    set(gca, 'YTick', []);
    title([zplotTitles{3}], 'FontSize', fontSz);
    hold off;

    nexttile
    hold on;
    imagesc(timescalePlot, 1:1:size(zplot4,1), zplot4)
    colormap(blueOrange)
    set(gca, 'YDir', 'reverse');
    axis tight;
    zeroXDottedLine;
    clim(climsPsth);
    xlabel('time (sec)', 'FontSize', fontSz);
    set(gca, 'YTick', []);
    title([zplotTitles{4}], 'FontSize', fontSz);
    hold off;

    nexttile
    hold on;
    imagesc(timescalePlot, 1:1:size(zplot5,1), zplot5)
    colormap(blueOrange)
    set(gca, 'YDir', 'reverse');
    axis tight;
    zeroXDottedLine;
    clim(climsPsth);
    cb = colorbar;
    cb.Label.String = ['z-score ' '\Delta' ' firing rate (Hz)'];
    cb.Label.FontSize = fontSz;
    xlabel('time (sec)', 'FontSize', fontSz);
    set(gca, 'YTick', []);
    title([zplotTitles{5}], 'FontSize', fontSz);
    hold off;

    set(gcf,'color','w'); %set figure background white

end

function plotBasicNeuronTuningGridPerProbe(allProbes, timescaleSeg, plotWindow, itiWindow) 
%plot ant-react order next to other events in grid layout
    plotSamps = [find(plotWindow(1)>=timescaleSeg,1,'last') find(plotWindow(2)>=timescaleSeg,1,'last')];
    itiSamps = [find(itiWindow(1)>=timescaleSeg,1,'last') find(itiWindow(2)>=timescaleSeg,1,'last')];
    timescalePlot = timescaleSeg(plotSamps(1):plotSamps(2));
    midPt = find(0>=timescaleSeg,1,'last'); % for sorting

    figure;
    numProbes = 3;
    unitsPerProbe = allProbes.numProbeTypeUnits;
    yAxesAspectRatios = unitsPerProbe ./ max(unitsPerProbe);
    cumSumUnits = cumsum(allProbes.numProbeTypeUnits);
    probeStartUnits = cumSumUnits - unitsPerProbe + 1; %index into start of each probe
    numPlots = 4;
    tiledlayout(numProbes,numPlots,'TileSpacing', 'compact');
    climsPsth = [-0.5 0.5];
    fontSz = 16;
    zplotTitles = {'@ avoid only', '@ react only', '@ cue', '@ cold air',};
    depthSort = 0;
    currentPlotOffset = 0;
    for probeIdx = [2 1 3] % control order of probes plotted here: go rostral to caudal

        %first break up units by probe
        allUnitsAnatDepths = allProbes.stDepths(probeStartUnits(probeIdx):cumSumUnits(probeIdx));
        allCueMeanZFR = allProbes.allUnitsPlatformExtendMeanZFR(probeStartUnits(probeIdx):cumSumUnits(probeIdx),:);
%         allCueMeanZFR = allProbes.allUnitsPlatformCommandMeanZFR(probeStartUnits(probeIdx):cumSumUnits(probeIdx),:);
        allColdAirZFR = allProbes.allUnitsColdAirOnMeanZFR(probeStartUnits(probeIdx):cumSumUnits(probeIdx),:);
        allVhexMeanZFR = allProbes.allUnitsMeanZFR(probeStartUnits(probeIdx):cumSumUnits(probeIdx),:);
        avoidMeanZFR = allProbes.allUnitsAvoidMeanZFR(probeStartUnits(probeIdx):cumSumUnits(probeIdx),:);
        reactMeanZFR = allProbes.allUnitsReactMeanZFR(probeStartUnits(probeIdx):cumSumUnits(probeIdx),:);

        % use depth sort or rastermap sort for each plot
        if depthSort
            [~,sortIdx] = sort(allUnitsAnatDepths, 'descend'); %'depths' are distances from probe tip, so smaller numbers are deepest in brain
            zplot3 = allCueMeanZFR(sortIdx,plotSamps(1):plotSamps(2));
            zplot4 = allColdAirZFR(sortIdx,plotSamps(1):plotSamps(2));
%             zplot3 = allVhexMeanZFR(sortIdx,plotSamps(1):plotSamps(2));
            zplot2 = reactMeanZFR(sortIdx,plotSamps(1):plotSamps(2));
            zplot1 = avoidMeanZFR(sortIdx,plotSamps(1):plotSamps(2));
        else % else can use beginning - end diff for each, or try rastermap
            [~,sortIdx] = sort(mean(avoidMeanZFR(:,itiSamps(1):itiSamps(2)),2)-mean(avoidMeanZFR(:,round(plotSamps(1)):midPt),2), 'descend'); %sort all by avoid modulation
            %[~,sortIdx1] = sort( mean(allCueMeanZFR(:,plotSamps(1):midPt),2) - mean(allCueMeanZFR(:,midPt+1:plotSamps(2)),2), 'descend');
            zplot3 = allCueMeanZFR(sortIdx,plotSamps(1):plotSamps(2));
            %[~,sortIdx2] = sort( mean(allColdAirZFR(:,plotSamps(1):midPt),2) - mean(allColdAirZFR(:,midPt+1:plotSamps(2)),2), 'descend');
            zplot4 = allColdAirZFR(sortIdx,plotSamps(1):plotSamps(2));
%             zplot3 = allVhexMeanZFR(sortIdx,plotSamps(1):plotSamps(2));
            %[~,sortIdx4] = sort(mean(reactMeanZFR(:,1:midPt),2)-mean(reactMeanZFR(:,midPt+1:end),2), 'descend');
            zplot2 = reactMeanZFR(sortIdx,plotSamps(1):plotSamps(2));
            [~,sortIdx5] = sort(mean(avoidMeanZFR(:,itiSamps(1):itiSamps(2)),2)-mean(avoidMeanZFR(:,round(plotSamps(1)/2):midPt),2), 'descend');
            zplot1 = avoidMeanZFR(sortIdx,plotSamps(1):plotSamps(2));
        end
        
        nexttile(currentPlotOffset+1)
        hold on;
        imagesc(timescalePlot, 1:1:size(zplot1,1), zplot1)
        pbaspect([1, yAxesAspectRatios(probeIdx), 1]); %set all units to equal row height regardless of units per probe
        colormap(blueOrange)
        set(gca, 'YDir', 'reverse');
        axis tight;
        zeroXDottedLine;
        clim(climsPsth);
        ylabel('unit #', 'FontSize', fontSz);
        if probeIdx ~= 3
            set(gca, 'XTick', []);
        end
        if probeIdx == 2
            title([zplotTitles{1}], 'FontSize', fontSz);
        end
        hold off;
    
        nexttile(currentPlotOffset+2)
        hold on;
        imagesc(timescalePlot, 1:1:size(zplot2,1), zplot2)
        pbaspect([1, yAxesAspectRatios(probeIdx), 1]);
        colormap(blueOrange)
        set(gca, 'YDir', 'reverse');
        axis tight;
        zeroXDottedLine;
        clim(climsPsth);
        set(gca, 'YTick', []);
        if probeIdx ~= 3
            set(gca, 'XTick', []);
        end
        if probeIdx == 2
            title([zplotTitles{2}], 'FontSize', fontSz);
        end
        hold off;
    
        nexttile(currentPlotOffset+3)
        hold on;
        imagesc(timescalePlot, 1:1:size(zplot3,1), zplot3)
        pbaspect([1, yAxesAspectRatios(probeIdx), 1]);
        colormap(blueOrange)
        set(gca, 'YDir', 'reverse');
        axis tight;
        zeroXDottedLine;
        clim(climsPsth);
        if probeIdx == 3
            xlabel('time (sec)', 'FontSize', fontSz);
        end
        set(gca, 'YTick', []);
        if probeIdx ~= 3
            set(gca, 'XTick', []);
        end
        if probeIdx == 2
            title([zplotTitles{3}], 'FontSize', fontSz);
        end
        hold off;
    
        nexttile(currentPlotOffset+4)
        hold on;
        imagesc(timescalePlot, 1:1:size(zplot4,1), zplot4)
        pbaspect([1, yAxesAspectRatios(probeIdx), 1]);
        colormap(blueOrange)
        set(gca, 'YDir', 'reverse');
        axis tight;
        zeroXDottedLine;
        clim(climsPsth);
        cb = colorbar;
        cb.Label.String = ['z-score  firing rate (Hz)'];
        cb.Label.FontSize = fontSz;
        set(gca, 'YTick', []);
        if probeIdx ~= 3
            set(gca, 'XTick', []);
        end
        if probeIdx == 2
            title([zplotTitles{4}], 'FontSize', fontSz);
        end
        hold off;

    currentPlotOffset = currentPlotOffset + numPlots;
    end

    set(gcf,'color','w'); %set figure background white

end

function plotAvoidNeuronTuningGrid(diffMeanZFR, reactMeanZFR, platformMeanZFR, coldAirOnMeanZFR, allUnitsAnatHiers, timescaleSeg, plotWindow, sortIdx) 
%plot ant-react order next to other events in grid layout
    plotSamps = [find(plotWindow(1)>=timescaleSeg,1,'last') find(plotWindow(2)>=timescaleSeg,1,'last')];
    diffMeanZFRPlot = diffMeanZFR(sortIdx,plotSamps(1):plotSamps(2));
    reactMeanZFRPlot = reactMeanZFR(sortIdx,plotSamps(1):plotSamps(2));
    platformMeanZFRPlot = platformMeanZFR(sortIdx,plotSamps(1):plotSamps(2));
    coldAirOnMeanZFRPlot = coldAirOnMeanZFR(sortIdx,plotSamps(1):plotSamps(2));
    timescalePlot = timescaleSeg(plotSamps(1):plotSamps(2));
    
    [unitCmaps, unitCmapsIdx, anatCmapsPlot, anatNamesPlot] = anatHiersToCmaps(allUnitsAnatHiers); 
    unitCmaps = unitCmaps(sortIdx,:);
    unitCmapsIdx = unitCmapsIdx(sortIdx,:);

    figure;
    tiledlayout(1,4,'TileSpacing', 'tight');
    climsPsth = [-0.5 0.5];
    fontSz = 16;

    nexttile
    hold on;
    imagesc(timescalePlot, 1:1:size(diffMeanZFRPlot,1), diffMeanZFRPlot)
    colormap(blueOrange)
    set(gca, 'YDir', 'reverse');
    %plot small colored horizontal lines to denote anatomy on the left side of plot:
    lineLength = 0.1;
    numColors = size(anatCmapsPlot,1);
    totalWidth = lineLength*numColors;    
    for unitIdx = 1:size(unitCmaps,1)
        colorShift = (unitCmapsIdx(unitIdx)/numColors)*totalWidth;
        x1 = timescalePlot(1)-colorShift; %slightly shift each line depending on color index
        x2 = timescalePlot(1)-colorShift+lineLength;
        plot([x1 x2], [unitIdx unitIdx], 'LineWidth', 0.1, 'Color', unitCmaps(unitIdx,:))
    end
    axis tight;
    zeroXDottedLine;
    clim(climsPsth);
    ylabel('unit #', 'FontSize', fontSz);
    xlabel('time (sec)', 'FontSize', fontSz);
    title('avoid - react @ VHEx', 'FontSize', fontSz);
%     drawAnatColormapLegend(anatCmapsPlot, anatNamesPlot)
    hold off;

    nexttile
    hold on;
    imagesc(timescalePlot, 1:1:size(reactMeanZFRPlot,1), reactMeanZFRPlot)
    colormap(blueOrange)
    set(gca, 'YDir', 'reverse');
    axis tight;
    zeroXDottedLine;
    clim(climsPsth);
    xlabel('time (sec)', 'FontSize', fontSz);
    set(gca, 'YTick', []);
    title('react @ VHEx', 'FontSize', fontSz);
    hold off;

    nexttile
    hold on;
    imagesc(timescalePlot, 1:1:size(platformMeanZFRPlot,1), platformMeanZFRPlot)
    colormap(blueOrange)
    set(gca, 'YDir', 'reverse');
    axis tight;
    zeroXDottedLine;
    clim(climsPsth);
    xlabel('time (sec)', 'FontSize', fontSz);
    set(gca, 'YTick', []);
    title('@ platform motor', 'FontSize', fontSz);
    hold off;

    nexttile
    hold on;
    imagesc(timescalePlot, 1:1:size(coldAirOnMeanZFRPlot,1), coldAirOnMeanZFRPlot)
    colormap(blueOrange)
    set(gca, 'YDir', 'reverse');
    axis tight;
    zeroXDottedLine;
    clim(climsPsth);
    cb = colorbar;
    cb.Label.String = ['z-score ' '\Delta' ' firing rate (Hz)'];
    cb.Label.FontSize = fontSz;
    xlabel('time (sec)', 'FontSize', fontSz);
    set(gca, 'YTick', []);
    title('@ beep off / cold on', 'FontSize', fontSz);
    hold off;

    set(gcf,'color','w'); %set figure background white

end

function plotAvoidAtCueNeuronTuningGrid(diffAtCueMeanZFR, diffAtVhexMeanZFR, platformMeanZFR, coldAirOnMeanZFR, avoidMeanZFR, reactMeanZFR, allUnitsAnatHiers, timescaleSeg, plotWindow, sortIdx) 
%plot ant-react order next to other events in grid layout
    plotSamps = [find(plotWindow(1)>=timescaleSeg,1,'last') find(plotWindow(2)>=timescaleSeg,1,'last')];
    diffAtCueMeanZFRPlot = diffAtCueMeanZFR(sortIdx,plotSamps(1):plotSamps(2));
    diffAtVhexMeanZFRPlot = diffAtVhexMeanZFR(sortIdx,plotSamps(1):plotSamps(2));
    avoidMeanZFRPlot = avoidMeanZFR(sortIdx,plotSamps(1):plotSamps(2));
    reactMeanZFRPlot = reactMeanZFR(sortIdx,plotSamps(1):plotSamps(2));
    platformMeanZFRPlot = platformMeanZFR(sortIdx,plotSamps(1):plotSamps(2));
    coldAirOnMeanZFRPlot = coldAirOnMeanZFR(sortIdx,plotSamps(1):plotSamps(2));
    timescalePlot = timescaleSeg(plotSamps(1):plotSamps(2));
    
    [unitCmaps, unitCmapsIdx, anatCmapsPlot, anatNamesPlot] = anatHiersToCmaps(allUnitsAnatHiers); 
    unitCmaps = unitCmaps(sortIdx,:);
    unitCmapsIdx = unitCmapsIdx(sortIdx,:);

    figure;
%     tiledlayout(1,6,'TileSpacing', 'tight');
    climsPsth = [-0.5 0.5];
    fontSz = 16;

%     nexttile
    hold on;
    imagesc(timescalePlot, 1:1:size(diffAtCueMeanZFRPlot,1), diffAtCueMeanZFRPlot)
    colormap(blueOrange)
    set(gca, 'YDir', 'reverse');
    %plot small colored horizontal lines to denote anatomy on the left side of plot:
%     lineLength = 0.1;
%     numColors = size(anatCmapsPlot,1);
%     totalWidth = lineLength*numColors;    
%     for unitIdx = 1:size(unitCmaps,1)
%         colorShift = (unitCmapsIdx(unitIdx)/numColors)*totalWidth;
%         x1 = timescalePlot(1)-colorShift; %slightly shift each line depending on color index
%         x2 = timescalePlot(1)-colorShift+lineLength;
%         plot([x1 x2], [unitIdx unitIdx], 'LineWidth', 0.1, 'Color', unitCmaps(unitIdx,:))
%     end
    axis tight;
    zeroXDottedLine;
    clim(climsPsth);
    ylabel('unit #', 'FontSize', fontSz);
    xlabel('time (sec)', 'FontSize', fontSz);
    title('avoid - react @ cue', 'FontSize', fontSz);
%     drawAnatColormapLegend(anatCmapsPlot, anatNamesPlot)
    hold off;
% 
%     nexttile
%     hold on;
%     imagesc(timescalePlot, 1:1:size(diffAtVhexMeanZFRPlot,1), diffAtVhexMeanZFRPlot)
%     colormap(blueOrange)
%     set(gca, 'YDir', 'reverse');
%     axis tight;
%     zeroXDottedLine;
%     clim(climsPsth);
%     xlabel('time (sec)', 'FontSize', fontSz);
%     set(gca, 'YTick', []);
%     title('avoid - react @ VHEx', 'FontSize', fontSz);
%     hold off;
% 
%     nexttile
%     hold on;
%     imagesc(timescalePlot, 1:1:size(platformMeanZFRPlot,1), platformMeanZFRPlot)
%     colormap(blueOrange)
%     set(gca, 'YDir', 'reverse');
%     axis tight;
%     zeroXDottedLine;
%     clim(climsPsth);
%     xlabel('time (sec)', 'FontSize', fontSz);
%     set(gca, 'YTick', []);
%     title('@ platform cue', 'FontSize', fontSz);
%     hold off;
% 
%     nexttile
%     hold on;
%     imagesc(timescalePlot, 1:1:size(avoidMeanZFRPlot,1), avoidMeanZFRPlot)
%     colormap(blueOrange)
%     set(gca, 'YDir', 'reverse');
%     axis tight;
%     zeroXDottedLine;
%     clim(climsPsth);
%     xlabel('time (sec)', 'FontSize', fontSz);
%     set(gca, 'YTick', []);
%     title('avoid only', 'FontSize', fontSz);
%     hold off;
% 
%     nexttile
%     hold on;
%     imagesc(timescalePlot, 1:1:size(reactMeanZFRPlot,1), reactMeanZFRPlot)
%     colormap(blueOrange)
%     set(gca, 'YDir', 'reverse');
%     axis tight;
%     zeroXDottedLine;
%     clim(climsPsth);
%     xlabel('time (sec)', 'FontSize', fontSz);
%     set(gca, 'YTick', []);
%     title('react only', 'FontSize', fontSz);
%     hold off;
% 
%     nexttile
%     hold on;
%     imagesc(timescalePlot, 1:1:size(coldAirOnMeanZFRPlot,1), coldAirOnMeanZFRPlot)
%     colormap(blueOrange)
%     set(gca, 'YDir', 'reverse');
%     axis tight;
%     zeroXDottedLine;
%     clim(climsPsth);
%     cb = colorbar;
%     cb.Label.String = ['z-score ' '\Delta' ' firing rate (Hz)'];
%     cb.Label.FontSize = fontSz;
%     xlabel('time (sec)', 'FontSize', fontSz);
%     set(gca, 'YTick', []);
%     title('@ cold air on', 'FontSize', fontSz);
%     hold off;

    set(gcf,'color','w'); %set figure background white

end

function plotAvoidAtVhexNeuronTuningGrid(diffAtVhexMeanZFR, diffAtCueMeanZFR, platformMeanZFR, coldAirOnMeanZFR, avoidMeanZFR, reactMeanZFR, allUnitsAnatHiers, timescaleSeg, plotWindow, sortIdx) 
%plot ant-react order next to other events in grid layout
    plotSamps = [find(plotWindow(1)>=timescaleSeg,1,'last') find(plotWindow(2)>=timescaleSeg,1,'last')];
    diffAtCueMeanZFRPlot = diffAtCueMeanZFR(sortIdx,plotSamps(1):plotSamps(2));
    diffAtVhexMeanZFRPlot = diffAtVhexMeanZFR(sortIdx,plotSamps(1):plotSamps(2));
    avoidMeanZFRPlot = avoidMeanZFR(sortIdx,plotSamps(1):plotSamps(2));
    reactMeanZFRPlot = reactMeanZFR(sortIdx,plotSamps(1):plotSamps(2));
    platformMeanZFRPlot = platformMeanZFR(sortIdx,plotSamps(1):plotSamps(2));
    coldAirOnMeanZFRPlot = coldAirOnMeanZFR(sortIdx,plotSamps(1):plotSamps(2));
    timescalePlot = timescaleSeg(plotSamps(1):plotSamps(2));
    
    [unitCmaps, unitCmapsIdx, anatCmapsPlot, anatNamesPlot] = anatHiersToCmaps(allUnitsAnatHiers); 
    unitCmaps = unitCmaps(sortIdx,:);
    unitCmapsIdx = unitCmapsIdx(sortIdx,:);

    figure;
%     tiledlayout(1,6,'TileSpacing', 'tight');
    climsPsth = [-0.5 0.5];
    fontSz = 16;

%     nexttile
    hold on;
    imagesc(timescalePlot, 1:1:size(diffAtVhexMeanZFRPlot,1), diffAtVhexMeanZFRPlot)
    colormap(blueOrange)
    set(gca, 'YDir', 'reverse');
    %plot small colored horizontal lines to denote anatomy on the left side of plot:
%     lineLength = 0.1;
%     numColors = size(anatCmapsPlot,1);
%     totalWidth = lineLength*numColors;    
%     for unitIdx = 1:size(unitCmaps,1)
%         colorShift = (unitCmapsIdx(unitIdx)/numColors)*totalWidth;
%         x1 = timescalePlot(1)-colorShift; %slightly shift each line depending on color index
%         x2 = timescalePlot(1)-colorShift+lineLength;
%         plot([x1 x2], [unitIdx unitIdx], 'LineWidth', 0.1, 'Color', unitCmaps(unitIdx,:))
%     end
    axis tight;
    zeroXDottedLine;
    clim(climsPsth);
    ylabel('unit #', 'FontSize', fontSz);
    xlabel('time (sec)', 'FontSize', fontSz);
    title('avoid - react @ VHEx', 'FontSize', fontSz);
%     drawAnatColormapLegend(anatCmapsPlot, anatNamesPlot)
    hold off;

%     nexttile
%     hold on;
%     imagesc(timescalePlot, 1:1:size(diffAtCueMeanZFRPlot,1), diffAtCueMeanZFRPlot)
%     colormap(blueOrange)
%     set(gca, 'YDir', 'reverse');
%     axis tight;
%     zeroXDottedLine;
%     clim(climsPsth);
%     xlabel('time (sec)', 'FontSize', fontSz);
%     set(gca, 'YTick', []);
%     title('avoid - react @ cue', 'FontSize', fontSz);
%     hold off;
% 
%     nexttile
%     hold on;
%     imagesc(timescalePlot, 1:1:size(platformMeanZFRPlot,1), platformMeanZFRPlot)
%     colormap(blueOrange)
%     set(gca, 'YDir', 'reverse');
%     axis tight;
%     zeroXDottedLine;
%     clim(climsPsth);
%     xlabel('time (sec)', 'FontSize', fontSz);
%     set(gca, 'YTick', []);
%     title('@ platform cue', 'FontSize', fontSz);
%     hold off;
% 
%     nexttile
%     hold on;
%     imagesc(timescalePlot, 1:1:size(avoidMeanZFRPlot,1), avoidMeanZFRPlot)
%     colormap(blueOrange)
%     set(gca, 'YDir', 'reverse');
%     axis tight;
%     zeroXDottedLine;
%     clim(climsPsth);
%     xlabel('time (sec)', 'FontSize', fontSz);
%     set(gca, 'YTick', []);
%     title('avoid only', 'FontSize', fontSz);
%     hold off;
% 
%     nexttile
%     hold on;
%     imagesc(timescalePlot, 1:1:size(reactMeanZFRPlot,1), reactMeanZFRPlot)
%     colormap(blueOrange)
%     set(gca, 'YDir', 'reverse');
%     axis tight;
%     zeroXDottedLine;
%     clim(climsPsth);
%     xlabel('time (sec)', 'FontSize', fontSz);
%     set(gca, 'YTick', []);
%     title('react only', 'FontSize', fontSz);
%     hold off;
% 
%     nexttile
%     hold on;
%     imagesc(timescalePlot, 1:1:size(coldAirOnMeanZFRPlot,1), coldAirOnMeanZFRPlot)
%     colormap(blueOrange)
%     set(gca, 'YDir', 'reverse');
%     axis tight;
%     zeroXDottedLine;
%     clim(climsPsth);
%     cb = colorbar;
%     cb.Label.String = ['z-score ' '\Delta' ' firing rate (Hz)'];
%     cb.Label.FontSize = fontSz;
%     xlabel('time (sec)', 'FontSize', fontSz);
%     set(gca, 'YTick', []);
%     title('@ beep off / cold on', 'FontSize', fontSz);
%     hold off;

    set(gcf,'color','w'); %set figure background white

end

function plotAvoidVersusReactNeuronModulations(Probe0, Probe1, Probe2, timescaleSeg)
    %function plotAvoidVersusReactNeuronModulations(diffAtCueMeanZFR, diffAtVhexMeanZFR,diffAtItiMeanZFR, anatHiers, timescaleSeg)
    % function to plot up/down modulations for avoid versus react trials, separated by anatomy and whether at cue/extend/ITI
    allData.Probe0 = Probe0; allData.Probe1 = Probe1; allData.Probe2 = Probe2;
    cueWindow = [-1 0.1 0.1 1]; % for cue & ITI just compare second before and after
    beforeSampsCue = [find(cueWindow(1)>=timescaleSeg,1,'last') find(cueWindow(2)>=timescaleSeg,1,'last')];
    afterSampsCue = [find(cueWindow(3)>=timescaleSeg,1,'last') find(cueWindow(4)>=timescaleSeg,1,'last')];
    vhexWindow = [-6 -5 -1.5 0.1]; % for VHEx compare ITI to just before extension threshold
    beforeSampsVhex = [find(vhexWindow(1)>=timescaleSeg,1,'last') find(vhexWindow(2)>=timescaleSeg,1,'last')];
    afterSampsVhex = [find(vhexWindow(3)>=timescaleSeg,1,'last') find(vhexWindow(4)>=timescaleSeg,1,'last')];
    
    % for all neurons, compute indeces at cue/extend/ITI that describe whether mean activity is greater or less during avoid versus react trials
    probeToPlot = {'Probe0','Probe1','Probe0','Probe2'}; %matched pairs with below
    anatToPlot = {'CTX','CTX','BS','HB'};
    labels = {'M1','PFC','TH/HY','HB'}; 
    lineCmaps = linspecer(size(probeToPlot,2));
    
    for i = 1:size(probeToPlot,2) %main loop over probe structures
        numUnits = size(allData.(probeToPlot{1,i}).(anatToPlot{1,i}).allUnitsAvoidMeanZFR,1);
        cueMods.(probeToPlot{1,i}).(anatToPlot{1,i}) = zeros(numUnits,1);
        vhexMods.(probeToPlot{1,i}).(anatToPlot{1,i}) = zeros(numUnits,1);
        itiMods.(probeToPlot{1,i}).(anatToPlot{1,i}) = zeros(numUnits,1);
    
        diffAtCueMeanZFR = allData.(probeToPlot{1,i}).(anatToPlot{1,i}).allUnitsAvoidAtCueMeanZFR - allData.(probeToPlot{1,i}).(anatToPlot{1,i}).allUnitsReactAtCueMeanZFR;
        diffAtVhexMeanZFR = allData.(probeToPlot{1,i}).(anatToPlot{1,i}).allUnitsAvoidMeanZFR - allData.(probeToPlot{1,i}).(anatToPlot{1,i}).allUnitsReactMeanZFR;
        diffAtItiMeanZFR = allData.(probeToPlot{1,i}).(anatToPlot{1,i}).allUnitsAvoidAtItiMeanZFR - allData.(probeToPlot{1,i}).(anatToPlot{1,i}).allUnitsReactAtItiMeanZFR;
        for j = 1:numUnits
            cueMods.(probeToPlot{1,i}).(anatToPlot{1,i})(j) = mean(diffAtCueMeanZFR(j,afterSampsCue(1):afterSampsCue(2))) - mean(diffAtCueMeanZFR(j,beforeSampsCue(1):beforeSampsCue(2)));
            vhexMods.(probeToPlot{1,i}).(anatToPlot{1,i})(j) =  mean(diffAtVhexMeanZFR(j,afterSampsVhex(1):afterSampsVhex(2))) - mean(diffAtVhexMeanZFR(j,beforeSampsVhex(1):beforeSampsVhex(2)));
            itiMods.(probeToPlot{1,i}).(anatToPlot{1,i})(j) = mean(diffAtItiMeanZFR(j,afterSampsCue(1):afterSampsCue(2))) - mean(diffAtItiMeanZFR(j,beforeSampsCue(1):beforeSampsCue(2)));
        end
    end
    
    f = figure; hold on;
    
    data = {[cueMods.(probeToPlot{1,1}).(anatToPlot{1,1}) vhexMods.(probeToPlot{1,1}).(anatToPlot{1,1}) itiMods.(probeToPlot{1,1}).(anatToPlot{1,1})],...
            [cueMods.(probeToPlot{1,2}).(anatToPlot{1,2}) vhexMods.(probeToPlot{1,2}).(anatToPlot{1,2}) itiMods.(probeToPlot{1,2}).(anatToPlot{1,2})],...
            [cueMods.(probeToPlot{1,3}).(anatToPlot{1,3}) vhexMods.(probeToPlot{1,3}).(anatToPlot{1,3}) itiMods.(probeToPlot{1,3}).(anatToPlot{1,3})],...
            [cueMods.(probeToPlot{1,4}).(anatToPlot{1,4}) vhexMods.(probeToPlot{1,4}).(anatToPlot{1,4}) itiMods.(probeToPlot{1,4}).(anatToPlot{1,4})]};
    conditionLabels = {'@ cue','@ VHEx','@ mid-ITI'};
    
    % full violin plots with white embedded boxplots for 3x2 data
    daviolinplot(data,'xtlabels',conditionLabels,'violin','full', 'colors', lineCmaps,'box',2,'outliers',0); 
    
    zeroYDottedLine;
    drawAnatColormapLegend(lineCmaps, labels);
    ylabel('avoid minus react mean difference (z-score)')
    hold off;
    f.Position(3:4) = [1200 300]; %set figure size appropriate for plots
    set(gcf,'color','w'); %set figure background white
    makepretty

end

function plotOptoResponseGrid(probe0ctxFR, probe0nonctxFR, probe1ctxFR, probe1nonctxFR, probe2ctxFR, probe2nonctxFR, timescaleSeg) 
%plot opto responses for separate probes in grid layout; NOTE that this is
%in Hz not zcored units
    figure;
    tiledlayout(3,2,'TileSpacing', 'tight');
    climsPsth = [-1 1];
    fontSz = 16;

    %PFC first
    nexttile
    hold on;
    imagesc(timescaleSeg, 1:1:size(probe1ctxFR,1), probe1ctxFR)
    colormap(blueOrange)
    set(gca, 'YDir', 'reverse');
    axis tight;
    zeroXDottedLine;
    clim(climsPsth);
    ylabel('unit #', 'FontSize', fontSz);
    xlabel('time (sec)', 'FontSize', fontSz);
    title('probe 1 ctx', 'FontSize', fontSz);
    hold off;

    nexttile
    hold on;
    imagesc(timescaleSeg, 1:1:size(probe1nonctxFR,1), probe1nonctxFR)
    colormap(blueOrange)
    set(gca, 'YDir', 'reverse');
    axis tight;
    zeroXDottedLine;
    clim(climsPsth);
    ylabel('unit #', 'FontSize', fontSz);
    xlabel('time (sec)', 'FontSize', fontSz);
    title('probe 1 non-ctx', 'FontSize', fontSz);
    hold off;
    
    %then M1
    nexttile
    hold on;
    imagesc(timescaleSeg, 1:1:size(probe0ctxFR,1), probe0ctxFR)
    colormap(blueOrange)
    set(gca, 'YDir', 'reverse');
    axis tight;
    zeroXDottedLine;
    clim(climsPsth);
    ylabel('unit #', 'FontSize', fontSz);
    xlabel('time (sec)', 'FontSize', fontSz);
    title('probe 0 ctx', 'FontSize', fontSz);
    hold off;

    nexttile
    hold on;
    imagesc(timescaleSeg, 1:1:size(probe0nonctxFR,1), probe0nonctxFR)
    colormap(blueOrange)
    set(gca, 'YDir', 'reverse');
    axis tight;
    zeroXDottedLine;
    clim(climsPsth);
    ylabel('unit #', 'FontSize', fontSz);
    xlabel('time (sec)', 'FontSize', fontSz);
    title('probe 0 non-ctx', 'FontSize', fontSz);
    hold off;

    %then V1
    nexttile
    hold on;
    imagesc(timescaleSeg, 1:1:size(probe2ctxFR,1), probe2ctxFR)
    colormap(blueOrange)
    set(gca, 'YDir', 'reverse');
    axis tight;
    zeroXDottedLine;
    clim(climsPsth);
    ylabel('unit #', 'FontSize', fontSz);
    xlabel('time (sec)', 'FontSize', fontSz);
    title('probe 2 ctx', 'FontSize', fontSz);
    hold off;

    nexttile
    hold on;
    imagesc(timescaleSeg, 1:1:size(probe2nonctxFR,1), probe2nonctxFR)
    colormap(blueOrange)
    set(gca, 'YDir', 'reverse');
    axis tight;
    zeroXDottedLine;
    clim(climsPsth);
    ylabel('unit #', 'FontSize', fontSz);
    xlabel('time (sec)', 'FontSize', fontSz);
    title('probe 2 non-ctx', 'FontSize', fontSz);
    hold off;

    set(gcf,'color','w'); %set figure background white

end

function plotEphysFR(tsData, idxData, shankData, sortIdx)
    % plot change in firing rate from session mean; pass sortIdx = 0 if not sorting
    if sortIdx ~= 0
        dataToPlot = shankData.allLocalFR(:,sortIdx);
    else
        dataToPlot = shankData.allLocalFR; %default depth order
    end

    figure;
%     subplot(2,1,1);
%     plot(tsData.timescale, tsData.optoVoltage.*(18/5), 'r-', 'LineWidth', 2)
%     ylabel('optostimulation (norm.)')
%     axis tight;
%     makepretty;
%     
%     subplot(2,1,2);
    hold on;
    imagesc(shankData.psthTimes, 1:1:size(dataToPlot,2), shankData.dataToPlot')

    set(gca, 'YDir', 'normal');
    colormap(parula)
    clim([0 30]);
    cb = colorbar;
    %cb.Position = cb.Position + [0.09 0 0 0];
    cb.Label.String = 'deltaFR (Hz)';
    ylabel('unit #');
    xlabel('time (sec)')

    %for dotted lines at extensions:
    jumpIndecesDottedLinesX = nan(1,length(idxData.jumpIndeces)*3);
    jumpIndecesDottedLinesX(1:3:end) = tsData.timescale(idxData.jumpIndeces);
    jumpIndecesDottedLinesX(2:3:end) = tsData.timescale(idxData.jumpIndeces);
    jumpIndecesDottedLinesY = nan(1,length(idxData.jumpIndeces)*3);
    jumpIndecesDottedLinesY(1:3:end) = 0;
    jumpIndecesDottedLinesY(2:3:end) = ones(1,length(idxData.jumpIndeces)).*size(shankData.allLocalFR,2);
    plot(jumpIndecesDottedLinesX,jumpIndecesDottedLinesY,'k--', 'LineWidth', 2);

    %for dotted lines at opto:
%     optoIndecesDottedLinesX = nan(1,length(idxData.optoOnIndeces)*3);
%     optoIndecesDottedLinesX(1:3:end) = tsData.timescale(idxData.optoOnIndeces);
%     optoIndecesDottedLinesX(2:3:end) = tsData.timescale(idxData.optoOnIndeces);
%     optoIndecesDottedLinesY = nan(1,length(idxData.optoOnIndeces)*3);
%     optoIndecesDottedLinesY(1:3:end) = 0;
%     optoIndecesDottedLinesY(2:3:end) = ones(1,length(idxData.optoOnIndeces)).*size(stCell,2);
%     plot(optoIndecesDottedLinesX,optoIndecesDottedLinesY,'r--', 'LineWidth', 2);    
    
    
    hold off;
    axis tight;
    makepretty;
    set(gcf,'color','w'); %set figure background white
end

function plotEphysFRzScore(tsData, idxData, shankData, sortIdx)
    % plot z-scored firing rate changes across session;
    if sortIdx ~= 0
        dataToPlot = shankData.allDeltaFRzScore(:,sortIdx);
    else
        dataToPlot = shankData.allDeltaFRzScore; %default depth order
    end


    figure;
%     subplot(2,1,1);
%     plot(tsData.timescale, tsData.optoVoltage, 'r-', 'LineWidth', 2)
%     ylabel('optostimulation (a.u.)')
%     axis tight;
%     makepretty;
%     
%     subplot(2,1,2);
    hold on;
    
    imagesc(shankData.psthTimes, 1:1:size(dataToPlot,2), dataToPlot')
    set(gca, 'YDir', 'normal');
    colormap(blueOrange)
    clim([-2 2]);
    cb = colorbar;
    %cb.Position = cb.Position + [0.09 0 0 0];
    cb.Label.String = 'z-score deltaFR';
    ylabel('unit #');
    xlabel('time (sec)')

    
    %for dotted lines at extensions:
    jumpIndecesDottedLinesX = nan(1,length(idxData.jumpIndeces)*3);
    jumpIndecesDottedLinesX(1:3:end) = tsData.timescale(idxData.jumpIndeces);
    jumpIndecesDottedLinesX(2:3:end) = tsData.timescale(idxData.jumpIndeces);
    jumpIndecesDottedLinesY = nan(1,length(idxData.jumpIndeces)*3);
    jumpIndecesDottedLinesY(1:3:end) = 0;
    jumpIndecesDottedLinesY(2:3:end) = ones(1,length(idxData.jumpIndeces)).*size(dataToPlot,2);
    plot(jumpIndecesDottedLinesX,jumpIndecesDottedLinesY, 'g--', 'LineWidth', 0.1);

    %for dotted lines at opto:
    optoIndecesDottedLinesX = nan(1,length(idxData.optoOnIndeces)*3);
    optoIndecesDottedLinesX(1:3:end) = tsData.timescale(idxData.optoOnIndeces);
    optoIndecesDottedLinesX(2:3:end) = tsData.timescale(idxData.optoOnIndeces);
    optoIndecesDottedLinesY = nan(1,length(idxData.optoOnIndeces)*3);
    optoIndecesDottedLinesY(1:3:end) = 0;
    optoIndecesDottedLinesY(2:3:end) = ones(1,length(idxData.optoOnIndeces)).*size(dataToPlot,2);
    plot(optoIndecesDottedLinesX,optoIndecesDottedLinesY,'r--', 'LineWidth', 0.5);    
    
    
    hold off;
    axis tight;
    makepretty;
    set(gcf,'color','w'); %set figure background white
end

function plotAvoidVsReactPcTrajectoriesAllMice(allUnitsMeanZFR, allUnitsAvoidMeanZFR, allUnitsReactMeanZFR, allUnitsAnatHiers, numGoodUnitsTotal, numGoodUnitsPerProbeType, timescaleSeg, pcaWindow, plotWindow, titleStr)
% input is all N x T matrices that concatenate trial averages (over avoid / react trials) across mice/sessions
% also input numGoodUnits matrix that keeps track of # neuron splits between mice/sessions, and time windows for PCA calculation and plotting
    numProbes = length(numGoodUnitsTotal);
    pcaSamps = [find(pcaWindow(1)>=timescaleSeg,1,'last') find(pcaWindow(2)>=timescaleSeg,1,'last')];
    plotSamps = [find(plotWindow(1)>=timescaleSeg,1,'last') find(plotWindow(2)>=timescaleSeg,1,'last')];
    timescalePlot = timescaleSeg(plotSamps(1):plotSamps(2));
    fontSz = 20;
    % empirically determined number of PCs to use if aligning PCs in shared subspace (using varExplained below), to keep constant across different # neurons in recordings
    numPcsToUse = 10;
    usePcContrast = 0; % whether we perform contrastive PCA on data1 vs data2, or vanilla PCA on allData 
    alignPCs = 1; % whether to align all mice/sessions into common subspace for computing PCs; otherwise PCs are computed independently for each mouse/session (default for simplicity)
    smoothingSamples = [10 0]; %ex. smooth window 25 is 25*20ms PSTH bins = 500ms; for causal half-kernel use [X 0]
    plot3d = 1;
    plotSinglePCs = 0;
    pcaIdx1 = 1; %make this variable to easy to change which PCs to view if necessary
    pcaIdx2 = 2;
    pcaIdx3 = 3;

    %plotting constants
    numColorPts = 256; %standard colormap pts
    R = linspace(0.3,0.9,numColorPts)'; %go from half-black to full color
    G = linspace(0.1,0.5,numColorPts)';
    B = zeros(numColorPts,1);
    cmap1 = [R G B];
    R = zeros(numColorPts,1);
    G = linspace(0.0,0.4,numColorPts)';
    B = linspace(0.3,0.9,numColorPts)';
    cmap2 = [R G B];
    cmapIdx = numColorPts - round(numColorPts/4); %for single colors for individual PC plots
    psub1 = [];
    psub2 = [];
    
    [pcSeg, varExplained] = getEphysPCs(allUnitsMeanZFR, allUnitsAvoidMeanZFR, allUnitsReactMeanZFR, numGoodUnitsTotal, numProbes, numPcsToUse, pcaSamps, usePcContrast, alignPCs); 

    % optional plot to look at variance explained & distribution of PC weights
%     hVarFig = figure; 
%     plotCap = 15; %cap the num PCs to plot to see detail
%     varMin = 40; %min % to plot
%     hold on;
%     if alignPCs
%         plot(cumsum(varExplained), 'k-', 'LineWidth', 3);
%     else
%         for i = 1:numProbes
%             plot(cumsum(varExplained{i}), 'r--', 'LineWidth', 0.5);
%             varExplainedAllProbes(i,:) = varExplained{i}(1:plotCap); %collect for mean
%         end
%         varExplainedMeanAcrossProbes = mean(varExplainedAllProbes,1);
%         plot(cumsum(varExplainedMeanAcrossProbes), 'r-', 'LineWidth', 3);
%     end
%     plot([numPcsToUse numPcsToUse], [varMin 100], 'b--');
%     hold off;
%     axis([0 plotCap varMin 100])
%     xlabel('PCs')
%     ylabel('% variance explained')
%     title(titleStr)
%     set(gcf,'color','w'); %set figure background white

    % --------------------------------------------------------------------------
%     % plot PC weights to help interpret PCs, with neurons colored by anatomy and separators drawn across probes 0/1/2
%     [unitCmaps, unitCmapsIdx, anatCmapsPlot, anatNamesPlot] = anatHiersToCmaps(allUnitsAnatHiers); 
%     pcMat = []; % empty matrix to combine PCs from cells for diff probes
%     if alignPCs
%         pcMat = cell2mat(pcSeg');
%     else
%         for i = 1:numProbes
%             pcMat = [pcMat; pcSeg{i}(:,1:numPcsToUse)];
%         end
%     end
%     hPcLoadFig = figure; 
%     tiledlayout(3,1,'TileSpacing', 'compact');
%     nexttile; hold on;
%     for i = 1:size(pcMat,1)
%         plot(i, pcMat(i,pcaIdx1), 'Marker', 'o', 'MarkerSize', 8, 'MarkerFaceColor', unitCmaps(i,:), 'MarkerEdgeColor', 'k') % plot each unit/point with separate color, takes a while but fine
%     end
%     if ~isempty(numGoodUnitsPerProbeType)
%         drawProbeSeparatorsPlot(numGoodUnitsPerProbeType);
%     end
%     ylabel(['PC' num2str(pcaIdx1) ' coeff.'])
%     hold off;
%     title(titleStr)
%     nexttile; hold on;
%     for i = 1:size(pcMat,1)
%         plot(i, pcMat(i,pcaIdx2), 'Marker', 'o', 'MarkerSize', 8, 'MarkerFaceColor', unitCmaps(i,:), 'MarkerEdgeColor', 'k') % plot each point with separate color, takes a while but fine
%     end
%     if ~isempty(numGoodUnitsPerProbeType)
%         drawProbeSeparatorsPlot(numGoodUnitsPerProbeType);
%     end
%     ylabel(['PC' num2str(pcaIdx2) ' coeff.'])
%     hold off;
%     nexttile; hold on;
%     for i = 1:size(pcMat,1)
%         plot(i, pcMat(i,pcaIdx3), 'Marker', 'o', 'MarkerSize', 8, 'MarkerFaceColor', unitCmaps(i,:), 'MarkerEdgeColor', 'k') % plot each point with separate color, takes a while but fine
%     end
%     if ~isempty(numGoodUnitsPerProbeType)
%         drawProbeSeparatorsPlot(numGoodUnitsPerProbeType);
%     end
%     ylabel(['PC' num2str(pcaIdx3) ' coeff.'])
%     xlabel('unit #')
%     drawAnatColormapLegend(anatCmapsPlot, anatNamesPlot)
%     hold off;
%     hPcLoadFig.Position(3:4) = [2000 1000]; %set figure size appropriate for plots
%     movegui(hPcLoadFig,'center')
%     set(gcf,'color','w'); %set figure background white  
    % --------------------------------------------------------------------------

    currentEndUnit = 0;
    if plotSinglePCs
        pc1OnlyFigure = figure; 
        tiledlayout(1,3,'TileSpacing', 'compact');
    end
    
    for i = 1:numProbes
        currentStartUnit = currentEndUnit + 1;
        currentEndUnit = currentEndUnit + numGoodUnitsTotal(i);
        if numGoodUnitsTotal(i) < numPcsToUse
            continue; % go to next mouse if this particular recording/region has less dimensions than we intend to use
        end

        % project mean avoid/ react activity for each mouse onto its respective PC matrix to get trajectories for each mouse & overlay
        avoidMeanProj = pcSeg{i}'*allUnitsAvoidMeanZFR(currentStartUnit:currentEndUnit,plotSamps(1):plotSamps(2)); %apply PCs to the mean responses across trials
        reactMeanProj = pcSeg{i}'*allUnitsReactMeanZFR(currentStartUnit:currentEndUnit,plotSamps(1):plotSamps(2));
        avoidMeanProj = smoothdata(avoidMeanProj, 2, 'gaussian', smoothingSamples); %ex. smooth window 25 is 25*20ms PSTH bins = 500ms
        reactMeanProj = smoothdata(reactMeanProj, 2, 'gaussian', smoothingSamples);
    
        colorsIdx = floor(linspace(1, numColorPts, size(avoidMeanProj,2))); %evenly spaced colors that show time dimension
    
        if plot3d
            figure; hold on; %uncomment this and axes labels etc. below for new PCA fig every time
            euclDistance = zeros(size(avoidMeanProj,2)-1,1);
            zeroSample = find(timescalePlot>0,1,'first'); %for plotting midpoint
    
            for dataPoint = 1:size(avoidMeanProj,2)-1
                if dataPoint == size(avoidMeanProj,2)-1
                    psub1 = plot3(avoidMeanProj(pcaIdx1,dataPoint:dataPoint+1),avoidMeanProj(pcaIdx2,dataPoint:dataPoint+1),avoidMeanProj(pcaIdx3,dataPoint:dataPoint+1),'Color',cmap1(colorsIdx(dataPoint),:), 'LineWidth', 3, 'DisplayName','avoid'); %subset for legend display name
                    psub2 = plot3(reactMeanProj(pcaIdx1,dataPoint:dataPoint+1),reactMeanProj(pcaIdx2,dataPoint:dataPoint+1),reactMeanProj(pcaIdx3,dataPoint:dataPoint+1),'Color',cmap2(colorsIdx(dataPoint),:),  'LineWidth', 3, 'DisplayName','react'); %subset for legend display name
                else
                    plot3(avoidMeanProj(pcaIdx1,dataPoint:dataPoint+1),avoidMeanProj(pcaIdx2,dataPoint:dataPoint+1),avoidMeanProj(pcaIdx3,dataPoint:dataPoint+1),'Color',cmap1(colorsIdx(dataPoint),:), 'LineWidth', 3); 
                    plot3(reactMeanProj(pcaIdx1,dataPoint:dataPoint+1),reactMeanProj(pcaIdx2,dataPoint:dataPoint+1),reactMeanProj(pcaIdx3,dataPoint:dataPoint+1),'Color',cmap2(colorsIdx(dataPoint),:),  'LineWidth', 3); 
                    if dataPoint == zeroSample %plot marker at t=0 extension threshold
                        plot3(avoidMeanProj(1,dataPoint),avoidMeanProj(2,dataPoint),avoidMeanProj(3,dataPoint), 'Marker', 'o', 'MarkerFaceColor','m', 'MarkerSize', 9, 'MarkerEdgeColor','none'); 
                        plot3(reactMeanProj(1,dataPoint),reactMeanProj(2,dataPoint),reactMeanProj(3,dataPoint),'Marker', 'o','MarkerFaceColor','m', 'MarkerSize', 9, 'MarkerEdgeColor','none'); 
                    end
                    if dataPoint == 1 %plot marker beginning
                        plot3(avoidMeanProj(1,dataPoint),avoidMeanProj(2,dataPoint),avoidMeanProj(3,dataPoint), 'Marker', 'o', 'MarkerFaceColor','g', 'MarkerSize', 9, 'MarkerEdgeColor','none'); 
                        plot3(reactMeanProj(1,dataPoint),reactMeanProj(2,dataPoint),reactMeanProj(3,dataPoint),'Marker', 'o','MarkerFaceColor','g', 'MarkerSize', 9, 'MarkerEdgeColor','none'); 
                    end
                end
            end % end dataPoint loop

            xlabel('PC1 (a.u.)')
            ylabel('PC2 (a.u.)')
            zlabel('PC3 (a.u.)')
            if ~isempty(psub1)
                legend([psub1 psub2 psub4 psub3], 'Box', 'off');
            end
            hold off;
            view(45,30);
            axis tight;
            title(titleStr)
            set(gcf,'color','w'); %set figure background white
        end %end if plot3d
        
        if plotSinglePCs
            % separate single PC trajectories plots - check for temporal structure in higher PCs
            figure(pc1OnlyFigure); 
            nexttile(1); hold on;
            plot(timescalePlot, avoidMeanProj(pcaIdx1,:),'Color',cmap1(cmapIdx,:), 'LineWidth', 0.1, 'DisplayName','avoid'); 
            plot(timescalePlot, reactMeanProj(pcaIdx1,:),'Color',cmap2(cmapIdx,:), 'LineWidth', 0.1, 'DisplayName','react'); 
            pc1ProjsAvoid(i,:) = avoidMeanProj(pcaIdx1,:); % save for bounded line / mean plot below
            pc1ProjsReact(i,:) = reactMeanProj(pcaIdx1,:);
            nexttile(2); hold on;
            plot(timescalePlot, avoidMeanProj(pcaIdx2,:),'Color',cmap1(cmapIdx,:), 'LineWidth', 0.1, 'DisplayName','avoid'); 
            plot(timescalePlot, reactMeanProj(pcaIdx2,:),'Color',cmap2(cmapIdx,:), 'LineWidth', 0.1, 'DisplayName','react'); 
            pc2ProjsAvoid(i,:) = avoidMeanProj(pcaIdx2,:); % save for bounded line / mean plot below
            pc2ProjsReact(i,:) = reactMeanProj(pcaIdx2,:);
            nexttile(3); hold on;
            plot(timescalePlot, avoidMeanProj(pcaIdx3,:),'Color',cmap1(cmapIdx,:), 'LineWidth', 0.1, 'DisplayName','avoid'); 
            plot(timescalePlot, reactMeanProj(pcaIdx3,:),'Color',cmap2(cmapIdx,:), 'LineWidth', 0.1, 'DisplayName','react'); 
            pc3ProjsAvoid(i,:) = avoidMeanProj(pcaIdx3,:); % save for bounded line / mean plot below
            pc3ProjsReact(i,:) = reactMeanProj(pcaIdx3,:);
        end
    end %end loop over mice
    
    if plotSinglePCs
        % finish single PC figures:
        figure(pc1OnlyFigure);
        nexttile(1); 
        [meanAvoidPC1, semAvoidPC1] = grpstats(pc1ProjsAvoid,[],{'mean' 'sem'});
        [meanReactPC1, semReactPC1] = grpstats(pc1ProjsReact,[],{'mean' 'sem'});
        plot(timescalePlot, meanAvoidPC1, 'Color',cmap1(cmapIdx,:), 'LineWidth', 3);
        plot(timescalePlot, meanReactPC1, 'Color',cmap2(cmapIdx,:), 'LineWidth', 3);
    %     boundedline(timescalePlot, meanAvoidPC1, semAvoidPC1, 'cmap', cmap1(cmapIdx,:));
    %     boundedline(timescalePlot, meanReactPC1, semReactPC1, 'cmap', cmap2(cmapIdx,:));
        zeroXDottedLine;
        hold off;
        xlabel('time (sec)', 'FontSize', fontSz)
        ylabel(['PC' num2str(pcaIdx1) ' (a.u.)'], 'FontSize', fontSz)
        axis tight;
    
        nexttile(2); 
        [meanAvoidPC2, semAvoidPC2] = grpstats(pc2ProjsAvoid,[],{'mean' 'sem'});
        [meanReactPC2, semReactPC2] = grpstats(pc2ProjsReact,[],{'mean' 'sem'});
        plot(timescalePlot, meanAvoidPC2, 'Color',cmap1(cmapIdx,:), 'LineWidth', 3);
        plot(timescalePlot, meanReactPC2, 'Color',cmap2(cmapIdx,:), 'LineWidth', 3);
    %     boundedline(timescalePlot, meanAvoidPC2, semAvoidPC2, 'cmap', cmap1(cmapIdx,:));
    %     boundedline(timescalePlot, meanReactPC2, semReactPC2, 'cmap', cmap2(cmapIdx,:));
        zeroXDottedLine;
        hold off;
        xlabel('time (sec)', 'FontSize', fontSz)
        ylabel(['PC' num2str(pcaIdx2) ' (a.u.)'], 'FontSize', fontSz)
        axis tight;
    
        nexttile(3); 
        [meanAvoidPC3, semAvoidPC3] = grpstats(pc3ProjsAvoid,[],{'mean' 'sem'});
        [meanReactPC3, semReactPC3] = grpstats(pc3ProjsReact,[],{'mean' 'sem'});
        plot(timescalePlot, meanAvoidPC3, 'Color',cmap1(cmapIdx,:), 'LineWidth', 3);
        plot(timescalePlot, meanReactPC3, 'Color',cmap2(cmapIdx,:), 'LineWidth', 3);
    %     boundedline(timescalePlot, meanAvoidPC3, semAvoidPC3, 'cmap', cmap1(cmapIdx,:));
    %     boundedline(timescalePlot, meanReactPC3, semReactPC3, 'cmap', cmap2(cmapIdx,:));
        zeroXDottedLine;
        hold off;
        xlabel('time (sec)', 'FontSize', fontSz)
        ylabel(['PC' num2str(pcaIdx3) ' (a.u.)'], 'FontSize', fontSz)
        axis tight;
        set(gcf,'color','w'); %set figure background white
        title(titleStr)
    end
end

function plotAvoidVsReactVsOptoPcTrajectoriesAllMice(allUnitsMeanZFR, allUnitsAvoidMeanZFR, allUnitsReactMeanZFR, allUnitsOptoMeanZFR, allUnitsAnatHiers, numGoodUnitsTotal, numGoodUnitsPerProbeType, timescaleSeg, pcaWindow, plotWindow, titleStr)
% input is all N x T matrices that concatenate trial averages (over avoid / react trials) across mice/sessions
% also input numGoodUnits matrix that keeps track of # neuron splits between mice/sessions, and time windows for PCA calculation and plotting
    numProbes = length(numGoodUnitsTotal);
    pcaSamps = [find(pcaWindow(1)>=timescaleSeg,1,'last') find(pcaWindow(2)>=timescaleSeg,1,'last')];
    plotSamps = [find(plotWindow(1)>=timescaleSeg,1,'last') find(plotWindow(2)>=timescaleSeg,1,'last')];
    timescalePlot = timescaleSeg(plotSamps(1):plotSamps(2));
    fontSz = 20;
    blueCmap = [0 0.4470 0.7410];
    orangeCmap = [0.8500 0.3250 0.0980];
    % empirically determined number of PCs to use if aligning PCs in shared subspace (using varExplained below), to keep constant across different # neurons in recordings
    numPcsToUse = 10;
    usePcContrast = 1; % whether we perform contrastive PCA on data1 vs data2, or vanilla PCA on allData 
    alignPCs = 1; % whether to align all mice/sessions into common subspace for computing PCs; otherwise PCs are computed independently for each mouse/session (default for simplicity)
    smoothingSamples = [25 0]; %ex. smooth window 25 is 25*20ms PSTH bins = 500ms; for causal half-kernel use [X 0]
    plot3d = 0;
    plotSinglePCs = 0;
    pcaIdx1 = 1; %make this variable to easy to change which PCs to view if necessary
    pcaIdx2 = 2;
    pcaIdx3 = 3;

    %plotting constants
    numColorPts = 256; %standard colormap pts
    R = linspace(0.3,0.9,numColorPts)'; %go from half-black to full color
    G = linspace(0.1,0.5,numColorPts)';
    B = zeros(numColorPts,1);
    cmap1 = [R G B];
    R = zeros(numColorPts,1);
    G = linspace(0.0,0.4,numColorPts)';
    B = linspace(0.3,0.9,numColorPts)';
    cmap2 = [R G B];
    R = ones(numColorPts,1)*0.1;
    G = ones(numColorPts,1)*0.1;
    B = ones(numColorPts,1)*0.1;
    cmap3 = [R G B];
    cmapIdx = numColorPts - round(numColorPts/4); %for single colors for individual PC plots
    psub1 = [];
    psub2 = [];

    [pcSeg, varExplained] = getEphysPCs(allUnitsMeanZFR, allUnitsAvoidMeanZFR, allUnitsReactMeanZFR, numGoodUnitsTotal, numProbes, numPcsToUse, pcaSamps, usePcContrast, alignPCs); 

    currentEndUnit = 0;
    if plotSinglePCs
        pc1OnlyFigure = figure; 
        tiledlayout(1,3,'TileSpacing', 'compact');
    end
    
    for i = 1:numProbes
        currentStartUnit = currentEndUnit + 1;
        currentEndUnit = currentEndUnit + numGoodUnitsTotal(i);
        if numGoodUnitsTotal(i) < numPcsToUse
            continue; % go to next mouse if this particular recording/region has less dimensions than we intend to use
        end

        % project mean avoid/ react activity for each mouse onto its respective PC matrix to get trajectories for each mouse & overlay
        avoidMeanProj = pcSeg{i}'*allUnitsAvoidMeanZFR(currentStartUnit:currentEndUnit,plotSamps(1):plotSamps(2)); %apply PCs to the mean responses across trials
        reactMeanProj = pcSeg{i}'*allUnitsReactMeanZFR(currentStartUnit:currentEndUnit,plotSamps(1):plotSamps(2));
        optoMeanProj = pcSeg{i}'*allUnitsOptoMeanZFR(currentStartUnit:currentEndUnit,plotSamps(1):plotSamps(2));
        avoidMeanProj = smoothdata(avoidMeanProj, 2, 'gaussian', smoothingSamples); %ex. smooth window 25 is 25*20ms PSTH bins = 500ms
        reactMeanProj = smoothdata(reactMeanProj, 2, 'gaussian', smoothingSamples);
        optoMeanProj = smoothdata(optoMeanProj, 2, 'gaussian', smoothingSamples);

        colorsIdx = floor(linspace(1, numColorPts, size(avoidMeanProj,2))); %evenly spaced colors that show time dimension
    
        if plot3d
            figure; hold on; %uncomment this and axes labels etc. below for new PCA fig every time
            euclDistanceOptoAvoid = zeros(size(avoidMeanProj,2)-1,1);
            zeroSample = find(timescalePlot>0,1,'first'); %for plotting midpoint
    
            for dataPoint = 1:size(avoidMeanProj,2)-1
                if dataPoint == size(avoidMeanProj,2)-1
                    psub1 = plot3(avoidMeanProj(pcaIdx1,dataPoint:dataPoint+1),avoidMeanProj(pcaIdx2,dataPoint:dataPoint+1),avoidMeanProj(pcaIdx3,dataPoint:dataPoint+1),'Color',cmap1(colorsIdx(dataPoint),:), 'LineWidth', 3, 'DisplayName','avoid'); %subset for legend display name
                    psub2 = plot3(reactMeanProj(pcaIdx1,dataPoint:dataPoint+1),reactMeanProj(pcaIdx2,dataPoint:dataPoint+1),reactMeanProj(pcaIdx3,dataPoint:dataPoint+1),'Color',cmap2(colorsIdx(dataPoint),:),  'LineWidth', 3, 'DisplayName','react'); %subset for legend display name
                    psub3 = plot3(optoMeanProj(pcaIdx1,dataPoint:dataPoint+1),optoMeanProj(pcaIdx2,dataPoint:dataPoint+1),optoMeanProj(pcaIdx3,dataPoint:dataPoint+1),'Color',[0 0 0],  'LineWidth', 3, 'DisplayName','opto'); %subset for legend display name
                else
                    plot3(avoidMeanProj(pcaIdx1,dataPoint:dataPoint+1),avoidMeanProj(pcaIdx2,dataPoint:dataPoint+1),avoidMeanProj(pcaIdx3,dataPoint:dataPoint+1),'Color',cmap1(colorsIdx(dataPoint),:), 'LineWidth', 3); 
                    plot3(reactMeanProj(pcaIdx1,dataPoint:dataPoint+1),reactMeanProj(pcaIdx2,dataPoint:dataPoint+1),reactMeanProj(pcaIdx3,dataPoint:dataPoint+1),'Color',cmap2(colorsIdx(dataPoint),:),  'LineWidth', 3); 
                    plot3(optoMeanProj(pcaIdx1,dataPoint:dataPoint+1),optoMeanProj(pcaIdx2,dataPoint:dataPoint+1),optoMeanProj(pcaIdx3,dataPoint:dataPoint+1),'Color',[0 0 0],  'LineWidth', 3);
                    if dataPoint == zeroSample %plot marker at t=0 extension threshold
                        plot3(avoidMeanProj(1,dataPoint),avoidMeanProj(2,dataPoint),avoidMeanProj(3,dataPoint), 'Marker', 'o', 'MarkerFaceColor','m', 'MarkerSize', 9, 'MarkerEdgeColor','none'); 
                        plot3(reactMeanProj(1,dataPoint),reactMeanProj(2,dataPoint),reactMeanProj(3,dataPoint),'Marker', 'o','MarkerFaceColor','m', 'MarkerSize', 9, 'MarkerEdgeColor','none'); 
                        plot3(optoMeanProj(1,dataPoint),optoMeanProj(2,dataPoint),optoMeanProj(3,dataPoint),'Marker', 'o','MarkerFaceColor','m', 'MarkerSize', 9, 'MarkerEdgeColor','none'); 
                    end
                    if dataPoint == 1 %plot marker beginning
                        plot3(avoidMeanProj(1,dataPoint),avoidMeanProj(2,dataPoint),avoidMeanProj(3,dataPoint), 'Marker', 'o', 'MarkerFaceColor','g', 'MarkerSize', 9, 'MarkerEdgeColor','none'); 
                        plot3(reactMeanProj(1,dataPoint),reactMeanProj(2,dataPoint),reactMeanProj(3,dataPoint),'Marker', 'o','MarkerFaceColor','g', 'MarkerSize', 9, 'MarkerEdgeColor','none'); 
                        plot3(optoMeanProj(1,dataPoint),optoMeanProj(2,dataPoint),optoMeanProj(3,dataPoint),'Marker', 'o','MarkerFaceColor','g', 'MarkerSize', 9, 'MarkerEdgeColor','none'); 
                    end
                end
            end % end dataPoint loop

            xlabel('PC1 (a.u.)')
            ylabel('PC2 (a.u.)')
            zlabel('PC3 (a.u.)')
            if ~isempty(psub1)
                legend([psub1 psub2 psub3], 'Box', 'off');
            end
            hold off;
            view(45,30);
            axis tight;
            title(titleStr)
            set(gcf,'color','w'); %set figure background white
        end %end if plot3d
        
        if plotSinglePCs
            % separate single PC trajectories plots - check for temporal structure in higher PCs
            figure(pc1OnlyFigure); 
            nexttile(1); hold on;
%             h(1) = plot(timescalePlot, avoidMeanProj(pcaIdx1,:),'Color',cmap1(cmapIdx,:), 'LineWidth', 0.1, 'DisplayName','avoid'); 
%             h(2) = plot(timescalePlot, reactMeanProj(pcaIdx1,:),'Color',cmap2(cmapIdx,:), 'LineWidth', 0.1, 'DisplayName','react'); 
%             h(3) = plot(timescalePlot, optoMeanProj(pcaIdx1,:),'Color',cmap3(cmapIdx,:), 'LineWidth', 0.1, 'DisplayName','opto');
            pc1ProjsAvoid(i,:) = avoidMeanProj(pcaIdx1,:); % save for bounded line / mean plot below
            pc1ProjsReact(i,:) = reactMeanProj(pcaIdx1,:);
            pc1ProjsOpto(i,:) = optoMeanProj(pcaIdx1,:);
            nexttile(2); hold on;
%             plot(timescalePlot, avoidMeanProj(pcaIdx2,:),'Color',cmap1(cmapIdx,:), 'LineWidth', 0.1, 'DisplayName','avoid'); 
%             plot(timescalePlot, reactMeanProj(pcaIdx2,:),'Color',cmap2(cmapIdx,:), 'LineWidth', 0.1, 'DisplayName','react'); 
%             plot(timescalePlot, optoMeanProj(pcaIdx2,:),'Color',cmap3(cmapIdx,:), 'LineWidth', 0.1, 'DisplayName','opto'); 
            pc2ProjsAvoid(i,:) = avoidMeanProj(pcaIdx2,:); % save for bounded line / mean plot below
            pc2ProjsReact(i,:) = reactMeanProj(pcaIdx2,:);
            pc2ProjsOpto(i,:) = optoMeanProj(pcaIdx2,:);
            nexttile(3); hold on;
%             plot(timescalePlot, avoidMeanProj(pcaIdx3,:),'Color',cmap1(cmapIdx,:), 'LineWidth', 0.1, 'DisplayName','avoid'); 
%             plot(timescalePlot, reactMeanProj(pcaIdx3,:),'Color',cmap2(cmapIdx,:), 'LineWidth', 0.1, 'DisplayName','react'); 
%             plot(timescalePlot, optoMeanProj(pcaIdx3,:),'Color',cmap3(cmapIdx,:), 'LineWidth', 0.1, 'DisplayName','opto'); 
            pc3ProjsAvoid(i,:) = avoidMeanProj(pcaIdx3,:); % save for bounded line / mean plot below
            pc3ProjsReact(i,:) = reactMeanProj(pcaIdx3,:);
            pc3ProjsOpto(i,:) = optoMeanProj(pcaIdx3,:);
        end
    end %end loop over mice
    
    if plotSinglePCs
        % finish single PC figures:
        figure(pc1OnlyFigure);
        nexttile(1); 
        [meanAvoidPC1, semAvoidPC1] = grpstats(pc1ProjsAvoid,[],{'mean' 'sem'});
        [meanReactPC1, semReactPC1] = grpstats(pc1ProjsReact,[],{'mean' 'sem'});
        [meanOptoPC1, semOptoPC1] = grpstats(pc1ProjsOpto,[],{'mean' 'sem'});
        h(1) = plot(timescalePlot, meanAvoidPC1, 'Color',cmap1(cmapIdx,:), 'LineWidth', 3, 'DisplayName','avoid');
        h(2) = plot(timescalePlot, meanReactPC1, 'Color',cmap2(cmapIdx,:), 'LineWidth', 3, 'DisplayName','react');
        h(3) = plot(timescalePlot, meanOptoPC1, 'Color',cmap3(cmapIdx,:), 'LineWidth', 3, 'DisplayName','opto');
    %     boundedline(timescalePlot, meanAvoidPC1, semAvoidPC1, 'cmap', cmap1(cmapIdx,:));
    %     boundedline(timescalePlot, meanReactPC1, semReactPC1, 'cmap', cmap2(cmapIdx,:));

        zeroXDottedLine;
        hold off;
        xlabel('time (sec)', 'FontSize', fontSz)
        ylabel(['PC' num2str(pcaIdx1) ' (a.u.)'], 'FontSize', fontSz)
        legend(h, 'Box', 'off');
        axis tight;
    
        nexttile(2); 
        [meanAvoidPC2, semAvoidPC2] = grpstats(pc2ProjsAvoid,[],{'mean' 'sem'});
        [meanReactPC2, semReactPC2] = grpstats(pc2ProjsReact,[],{'mean' 'sem'});
        [meanOptoPC2, semOptoPC2] = grpstats(pc2ProjsOpto,[],{'mean' 'sem'});
        plot(timescalePlot, meanAvoidPC2, 'Color',cmap1(cmapIdx,:), 'LineWidth', 3);
        plot(timescalePlot, meanReactPC2, 'Color',cmap2(cmapIdx,:), 'LineWidth', 3);
        plot(timescalePlot, meanOptoPC2, 'Color',cmap3(cmapIdx,:), 'LineWidth', 3);
    %     boundedline(timescalePlot, meanAvoidPC2, semAvoidPC2, 'cmap', cmap1(cmapIdx,:));
    %     boundedline(timescalePlot, meanReactPC2, semReactPC2, 'cmap', cmap2(cmapIdx,:));

        zeroXDottedLine;
        hold off;
        xlabel('time (sec)', 'FontSize', fontSz)
        ylabel(['PC' num2str(pcaIdx2) ' (a.u.)'], 'FontSize', fontSz)
        axis tight;
    
        nexttile(3); 
        [meanAvoidPC3, semAvoidPC3] = grpstats(pc3ProjsAvoid,[],{'mean' 'sem'});
        [meanReactPC3, semReactPC3] = grpstats(pc3ProjsReact,[],{'mean' 'sem'});
        [meanOptoPC3, semOptoPC3] = grpstats(pc3ProjsOpto,[],{'mean' 'sem'});
        plot(timescalePlot, meanAvoidPC3, 'Color',cmap1(cmapIdx,:), 'LineWidth', 3);
        plot(timescalePlot, meanReactPC3, 'Color',cmap2(cmapIdx,:), 'LineWidth', 3);
        plot(timescalePlot, meanOptoPC3, 'Color',cmap3(cmapIdx,:), 'LineWidth', 3);
    %     boundedline(timescalePlot, meanAvoidPC3, semAvoidPC3, 'cmap', cmap1(cmapIdx,:));
    %     boundedline(timescalePlot, meanReactPC3, semReactPC3, 'cmap', cmap2(cmapIdx,:));

        zeroXDottedLine;
        hold off;
        xlabel('time (sec)', 'FontSize', fontSz)
        ylabel(['PC' num2str(pcaIdx3) ' (a.u.)'], 'FontSize', fontSz)
        axis tight;
        set(gcf,'color','w'); %set figure background white
        title(titleStr)
    end

    savedOptoAvoidEuclDistance = [];
    savedOptoReactEuclDistance = [];
    for i = 1:numProbes
        if numGoodUnitsTotal(i) < numPcsToUse
            continue; % go to next mouse if this particular recording/region has less dimensions than we intend to use
        end
    
        euclDistanceOptoAvoid = zeros(size(avoidMeanProj,2)-1,1);
        euclDistanceOptoReact = zeros(size(avoidMeanProj,2)-1,1);
        for dataPoint = 1:size(avoidMeanProj,2)-1
            euclDistanceOptoAvoid(dataPoint) = norm(optoMeanProj(:,dataPoint) - avoidMeanProj(:,dataPoint));
            euclDistanceOptoReact(dataPoint) = norm(optoMeanProj(:,dataPoint) - reactMeanProj(:,dataPoint));
        end
    
        % calculate Euclidean distance between trajectories over time, normalized by sqrt(# neurons recorded), since Euclidean distance scales by the sqrt of the number of dimensions
%         plot(timescalePlot(1:end-1), euclDistance./sqrt(numGoodUnits(i)))
        savedOptoAvoidEuclDistance(i, :) = euclDistanceOptoAvoid./sqrt(numGoodUnitsTotal(i)); %save for bounded line below
        savedOptoReactEuclDistance(i, :) = euclDistanceOptoReact./sqrt(numGoodUnitsTotal(i));
    end %end loop over mice

    figure; hold on;
    [meanOptoAvoidEuclDistance, semOptoAvoidEuclDistance] = grpstats(savedOptoAvoidEuclDistance,[],{'mean' 'sem'}); %mean across columns/timepoints
    [h1,~] = boundedline(timescalePlot(1:end-1), meanOptoAvoidEuclDistance, semOptoAvoidEuclDistance, 'cmap', orangeCmap, 'alpha','y', 'transparency', 0.5);
    [meanOptoReactEuclDistance, semOptoReactEuclDistance] = grpstats(savedOptoReactEuclDistance,[],{'mean' 'sem'}); %mean across columns/timepoints
    [h2,~] = boundedline(timescalePlot(1:end-1), meanOptoReactEuclDistance, semOptoReactEuclDistance, 'cmap', blueCmap, 'alpha','y', 'transparency', 0.5);
    ylim([0 1.2])
    zeroXDottedLine;
    hold off;
    legend([h1, h2], {'opto - avoid', 'opto - react'}, 'Box', 'off');
    xlabel('time (sec)', 'FontSize', fontSz)
    ylabel('distance between PC trajectories (z-Hz/neuron)', 'FontSize', fontSz)
    set(gcf,'color','w'); %set figure background white
    title(titleStr);
    makepretty;

end

function plotAvoidVsReactVsOptoPcTrajectoriesAllMiceEqualTrials(allUnitsMeanZFR, allUnitsAvoidMeanZFR, allUnitsReactMeanZFR, allUnitsOptoMeanZFR,  allUnitsMeanZFR_shuf1, allUnitsMeanZFR_shuf2, allUnitsAnatHiers, numGoodUnitsTotal, numGoodUnitsPerProbeType, timescaleSeg, pcaWindow, plotWindow, titleStr)
% % input is all #repeats x N x T matrices that concatenate trial averages (over avoid / react /opto trials) across mice/sessions
% only 1 of the 2 shuffles is used, since we're looking at difference between opto and avoid/react/shuffle
% also input numGoodUnits matrix that keeps track of # neuron splits between mice/sessions, and time windows for PCA calculation and plotting

    numMiceProbes = length(numGoodUnitsTotal);
    pcaSamps = [find(pcaWindow(1)>=timescaleSeg,1,'last') find(pcaWindow(2)>=timescaleSeg,1,'last')];
    plotSamps = [find(plotWindow(1)>=timescaleSeg,1,'last') find(plotWindow(2)>=timescaleSeg,1,'last')];
    timescalePlot = timescaleSeg(plotSamps(1):plotSamps(2));

    usePcContrast = 1; % whether we perform contrastive PCA on data1 vs data2, or vanilla PCA on allData 
    numPcsToUse = 10;
    alignPCs = 1; % whether to align all mice/sessions into common subspace for computing PCs; otherwise PCs are computed independently for each mouse/session (default for simplicity)
    smoothingSamples = [25 0]; %ex. smooth window 25 is 25*20ms PSTH bins = 500ms; for causal half-kernel use [X 0]
    %plotting constants
    fontSz = 20;
    blueCmap = [0 0.4470 0.7410];
    orangeCmap = [0.8500 0.3250 0.0980];
    optoCmap = [0 0 0];
    plot3d = 1;
    pcaIdx1 = 1; %make this variable to easy to change which PCs to view if necessary
    pcaIdx2 = 2;
    pcaIdx3 = 3;

    savedOptoAvoidEuclDistanceRpts = [];
    savedOptoReactEuclDistanceRpts = [];
    savedOptoShuffleEuclDistanceRpts = [];

    numRepeats = 10; %should match number in other EqualTrials function
    for rpt = 1:numRepeats 
        % for each sampling repeat, take get PCs and distances for each mouse
        [pcSeg, varExplained] = getEphysPCs(allUnitsMeanZFR, squeeze(allUnitsAvoidMeanZFR(rpt,:,:)), squeeze(allUnitsReactMeanZFR(rpt,:,:)), numGoodUnitsTotal, numMiceProbes, numPcsToUse, pcaSamps, usePcContrast, alignPCs); 
    
        currentEndUnit = 0;      
        for i = 1:numMiceProbes
            currentStartUnit = currentEndUnit + 1;
            currentEndUnit = currentEndUnit + numGoodUnitsTotal(i);
            if numGoodUnitsTotal(i) < numPcsToUse
                continue; % go to next mouse if this particular recording/region has less dimensions than we intend to use
            end

            % project mean avoid/ react activity for each mouse onto its respective PC matrix to get trajectories for each mouse & overlay
            avoidMeanProj = pcSeg{i}'*squeeze(allUnitsAvoidMeanZFR(rpt,currentStartUnit:currentEndUnit,plotSamps(1):plotSamps(2))); %apply PCs to the mean responses across trials
            reactMeanProj = pcSeg{i}'*squeeze(allUnitsReactMeanZFR(rpt,currentStartUnit:currentEndUnit,plotSamps(1):plotSamps(2)));
            optoMeanProj = pcSeg{i}'*squeeze(allUnitsOptoMeanZFR(rpt,currentStartUnit:currentEndUnit,plotSamps(1):plotSamps(2)));
            shuf1MeanProj = pcSeg{i}'*squeeze(allUnitsMeanZFR_shuf1(rpt,currentStartUnit:currentEndUnit,plotSamps(1):plotSamps(2))); %option to apply same PCs as above to shuffled means
            avoidMeanProj = smoothdata(avoidMeanProj, 2, 'gaussian', smoothingSamples); %ex. smooth window 25 is 25*20ms PSTH bins = 500ms
            reactMeanProj = smoothdata(reactMeanProj, 2, 'gaussian', smoothingSamples);
            optoMeanProj = smoothdata(optoMeanProj, 2, 'gaussian', smoothingSamples);
            shuf1MeanProj = smoothdata(shuf1MeanProj, 2, 'gaussian', smoothingSamples); 

            euclDistanceOptoAvoid = zeros(size(avoidMeanProj,2)-1,1);
            euclDistanceOptoReact = zeros(size(avoidMeanProj,2)-1,1);
            euclDistanceOptoShuffle = zeros(size(avoidMeanProj,2)-1,1);
            for dataPoint = 1:size(avoidMeanProj,2)-1
                euclDistanceOptoAvoid(dataPoint) = norm(optoMeanProj(:,dataPoint) - avoidMeanProj(:,dataPoint));
                euclDistanceOptoReact(dataPoint) = norm(optoMeanProj(:,dataPoint) - reactMeanProj(:,dataPoint));
                euclDistanceOptoShuffle(dataPoint) = norm(optoMeanProj(:,dataPoint) - shuf1MeanProj(:,dataPoint));
            end
       
            if plot3d
                if rpt == 1 %just do 1 repeat to avoid excessive figure generation
                    figure; hold on; %plot example trajectory fig every mouse
                    euclDistanceOptoAvoid = zeros(size(avoidMeanProj,2)-1,1);
                    zeroSample = find(timescalePlot>0,1,'first'); %for plotting midpoint
            
                    for dataPoint = 1:size(avoidMeanProj,2)-1
                        if dataPoint == size(avoidMeanProj,2)-1
                            psub1 = plot3(avoidMeanProj(pcaIdx1,dataPoint:dataPoint+1),avoidMeanProj(pcaIdx2,dataPoint:dataPoint+1),avoidMeanProj(pcaIdx3,dataPoint:dataPoint+1),'Color',orangeCmap, 'LineWidth', 3, 'DisplayName','avoid'); %subset for legend display name
                            psub2 = plot3(reactMeanProj(pcaIdx1,dataPoint:dataPoint+1),reactMeanProj(pcaIdx2,dataPoint:dataPoint+1),reactMeanProj(pcaIdx3,dataPoint:dataPoint+1),'Color',blueCmap,  'LineWidth', 3, 'DisplayName','react'); %subset for legend display name
                            psub3 = plot3(optoMeanProj(pcaIdx1,dataPoint:dataPoint+1),optoMeanProj(pcaIdx2,dataPoint:dataPoint+1),optoMeanProj(pcaIdx3,dataPoint:dataPoint+1),'Color',optoCmap,  'LineWidth', 3, 'DisplayName','opto'); %subset for legend display name
                        else
                            plot3(avoidMeanProj(pcaIdx1,dataPoint:dataPoint+1),avoidMeanProj(pcaIdx2,dataPoint:dataPoint+1),avoidMeanProj(pcaIdx3,dataPoint:dataPoint+1),'Color',orangeCmap, 'LineWidth', 3); 
                            plot3(reactMeanProj(pcaIdx1,dataPoint:dataPoint+1),reactMeanProj(pcaIdx2,dataPoint:dataPoint+1),reactMeanProj(pcaIdx3,dataPoint:dataPoint+1),'Color',blueCmap,  'LineWidth', 3); 
                            plot3(optoMeanProj(pcaIdx1,dataPoint:dataPoint+1),optoMeanProj(pcaIdx2,dataPoint:dataPoint+1),optoMeanProj(pcaIdx3,dataPoint:dataPoint+1),'Color',optoCmap,  'LineWidth', 3);
                            if dataPoint == zeroSample %plot marker at t=0 extension threshold
                                plot3(avoidMeanProj(1,dataPoint),avoidMeanProj(2,dataPoint),avoidMeanProj(3,dataPoint), 'Marker', 'o', 'MarkerFaceColor','m', 'MarkerSize', 9, 'MarkerEdgeColor','none'); 
                                plot3(reactMeanProj(1,dataPoint),reactMeanProj(2,dataPoint),reactMeanProj(3,dataPoint),'Marker', 'o','MarkerFaceColor','m', 'MarkerSize', 9, 'MarkerEdgeColor','none'); 
                                plot3(optoMeanProj(1,dataPoint),optoMeanProj(2,dataPoint),optoMeanProj(3,dataPoint),'Marker', 'o','MarkerFaceColor','m', 'MarkerSize', 9, 'MarkerEdgeColor','none'); 
                            end
                            if dataPoint == 1 %plot marker beginning
                                plot3(avoidMeanProj(1,dataPoint),avoidMeanProj(2,dataPoint),avoidMeanProj(3,dataPoint), 'Marker', 'o', 'MarkerFaceColor','g', 'MarkerSize', 9, 'MarkerEdgeColor','none'); 
                                plot3(reactMeanProj(1,dataPoint),reactMeanProj(2,dataPoint),reactMeanProj(3,dataPoint),'Marker', 'o','MarkerFaceColor','g', 'MarkerSize', 9, 'MarkerEdgeColor','none'); 
                                plot3(optoMeanProj(1,dataPoint),optoMeanProj(2,dataPoint),optoMeanProj(3,dataPoint),'Marker', 'o','MarkerFaceColor','g', 'MarkerSize', 9, 'MarkerEdgeColor','none'); 
                            end
                        end
                    end % end dataPoint loop
        
                    xlabel('PC1 (a.u.)')
                    ylabel('PC2 (a.u.)')
                    zlabel('PC3 (a.u.)')
                    if ~isempty(psub1)
                        legend([psub1 psub2 psub3], 'Box', 'off');
                    end
                    hold off;
                    view(45,30);
                    axis tight;
                    title(titleStr)
                    set(gcf,'color','w'); %set figure background white
                end
            end %end if plot3d


            savedOptoAvoidEuclDistanceRpts(rpt, i, :) = euclDistanceOptoAvoid./sqrt(numGoodUnitsTotal(i)); %save for bounded line below
            savedOptoReactEuclDistanceRpts(rpt, i, :) = euclDistanceOptoReact./sqrt(numGoodUnitsTotal(i));
            savedOptoShuffleEuclDistanceRpts(rpt, i, :) = euclDistanceOptoShuffle./sqrt(numGoodUnitsTotal(i));

        end %end loop over mice
    end %end sampling repeats

    % now take mean across repeats for each mouse:
    savedOptoAvoidEuclDistance = squeeze(mean(savedOptoAvoidEuclDistanceRpts,1));
    savedOptoReactEuclDistance = squeeze(mean(savedOptoReactEuclDistanceRpts,1));
    savedOptoShuffleEuclDistance = squeeze(mean(savedOptoShuffleEuclDistanceRpts,1));

    figure; hold on;
%     % optionally plot individual mice to see variability:
%     for i = 1:numMiceProbes
%         plot(timescalePlot(1:end-1), savedOptoAvoidEuclDistance(i,:), 'Color', orangeCmap, 'LineWidth', 0.1)
%         plot(timescalePlot(1:end-1), savedOptoReactEuclDistance(i,:), 'Color', blueCmap, 'LineWidth', 0.1)
% %         plot(timescalePlot(1:end-1), savedOptoShuffleEuclDistance(i,:), 'Color', [0 0 0], 'LineWidth', 0.1)
%     end
    [meanOptoAvoidEuclDistance, semOptoAvoidEuclDistance] = grpstats(savedOptoAvoidEuclDistance,[],{'mean' 'sem'}); %mean across columns/timepoints
    [h1,~] = boundedline(timescalePlot(1:end-1), meanOptoAvoidEuclDistance, semOptoAvoidEuclDistance, 'cmap', orangeCmap, 'alpha','y', 'transparency', 0.5);
    [meanOptoReactEuclDistance, semOptoReactEuclDistance] = grpstats(savedOptoReactEuclDistance,[],{'mean' 'sem'}); %mean across columns/timepoints
    [h2,~] = boundedline(timescalePlot(1:end-1), meanOptoReactEuclDistance, semOptoReactEuclDistance, 'cmap', blueCmap, 'alpha','y', 'transparency', 0.5);
%     [meanOptoShuffleEuclDistance, semOptoShuffleEuclDistance] = grpstats(savedOptoShuffleEuclDistance,[],{'mean' 'sem'}); %mean across columns/timepoints
%     [h3,~] = boundedline(timescalePlot(1:end-1), meanOptoShuffleEuclDistance, semOptoShuffleEuclDistance, 'cmap', [0 0 0], 'alpha','y', 'transparency', 0.5);
    ylim([0 0.7])
    xlim([timescalePlot(1)+0.5 timescalePlot(end)]) %adjust for smoothing edge effects
    zeroXDottedLine;
    hold off;
    legend([h1 h2], {'opto - avoid', 'opto - react'}, 'Box', 'off');
%     legend([h1 h2 h3], {'opto - avoid', 'opto - react', 'opto - shuffle'}, 'Box', 'off'); %shuffle doesn't add much here because always in between avoid and react
    xlabel('time (sec)', 'FontSize', fontSz)
    ylabel('distance between PC trajectories (z-Hz/neuron)', 'FontSize', fontSz)
    set(gcf,'color','w'); %set figure background white
    title(titleStr);
    makepretty;

end

function plotAvoidVsReactPcDistancesAllMice(allUnitsMeanZFR, allUnitsAvoidMeanZFR, allUnitsReactMeanZFR, allUnitsMeanZFR_shuf1, allUnitsMeanZFR_shuf2, numGoodUnits, timescaleSeg, pcaWindow, plotWindow, titleStr)
% input is all N x T matrices that concatenate trial averages (over avoid / react trials) across mice/sessions
% also input numGoodUnits matrix that keeps track of # neuron splits between mice/sessions, and time windows for PCA calculation and plotting
    numMiceProbes = length(numGoodUnits);
    pcaSamps = [find(pcaWindow(1)>=timescaleSeg,1,'last') find(pcaWindow(2)>=timescaleSeg,1,'last')];
    plotSamps = [find(plotWindow(1)>=timescaleSeg,1,'last') find(plotWindow(2)>=timescaleSeg,1,'last')];
    timescalePlot = timescaleSeg(plotSamps(1):plotSamps(2));

    usePcContrast = 1; % whether we perform contrastive PCA on data1 vs data2, or vanilla PCA on allData 
    numPcsToUse = 10;
    alignPCs = 1; % whether to align all mice/sessions into common subspace for computing PCs; otherwise PCs are computed independently for each mouse/session (default for simplicity)
    smoothingSamples = [25 0]; %ex. smooth window 25 is 25*20ms PSTH bins = 500ms; for causal half-kernel use [X 0]
    %plotting constants
    fontSz = 16;

    [pcSeg, varExplained] = getEphysPCs(allUnitsMeanZFR, allUnitsAvoidMeanZFR, allUnitsReactMeanZFR, numGoodUnits, numMiceProbes, numPcsToUse, pcaSamps, usePcContrast, alignPCs); 

    savedEuclDistance = [];
    currentEndUnit = 0;
    euclFig = figure; hold on;
    for i = 1:numMiceProbes %plot all probes separately
        currentStartUnit = currentEndUnit + 1;
        currentEndUnit = currentEndUnit + numGoodUnits(i);
        if numGoodUnits(i) < numPcsToUse
            continue; % go to next mouse if this particular recording/region has less dimensions than we intend to use
        end

        % project mean avoid/ react activity for each mouse onto its respective PC matrix to get trajectories for each mouse & overlay
        avoidMeanProj = pcSeg{i}'*allUnitsAvoidMeanZFR(currentStartUnit:currentEndUnit,plotSamps(1):plotSamps(2)); %apply PCs to the mean responses across trials
        reactMeanProj = pcSeg{i}'*allUnitsReactMeanZFR(currentStartUnit:currentEndUnit,plotSamps(1):plotSamps(2));
        avoidMeanProj = smoothdata(avoidMeanProj, 2, 'gaussian', smoothingSamples); %ex. smooth window 25 is 25*20ms PSTH bins = 500ms
        reactMeanProj = smoothdata(reactMeanProj, 2, 'gaussian', smoothingSamples);
    
        euclDistance = zeros(size(avoidMeanProj,2)-1,1);
        for dataPoint = 1:size(avoidMeanProj,2)-1
            euclDistance(dataPoint) = norm(avoidMeanProj(:,dataPoint) - reactMeanProj(:,dataPoint));
        end
    
        % calculate Euclidean distance between trajectories over time, normalized by sqrt(# neurons recorded), since Euclidean distance scales by the sqrt of the number of dimensions
        % plot(timescalePlot(1:end-1), euclDistance./sqrt(numGoodUnits(i)), 'Color', [0 0 1], 'LineWidth', 0.1)
        savedEuclDistance(i, :) = euclDistance./sqrt(numGoodUnits(i)); %save for bounded line below
    end %end loop over mice
   
    if ~isempty(savedEuclDistance)
        [meanEuclDistance, semEuclDistance] = grpstats(savedEuclDistance,[],{'mean' 'sem'}); %mean across columns/timepoints
        [hl1,~] = boundedline(timescalePlot(1:end-1), meanEuclDistance, semEuclDistance, 'cmap', [0 0 1], 'alpha','y', 'transparency', 0.5, 'LineWidth', 3);
    end


%%%%%%%%%%%%%%%% NOW PLOT CONTROL / SHUFFLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use same PC coefficients as above, but just project onto shuffled trial averages

% optionally used the shuffled trial means to create their own PCs and respective distances, but this is NOT a fair comparison:
%[pcSegShuf, varExplainedShuf] = getEphysPCs(allUnitsMeanZFR, allUnitsMeanZFR_shuf1, allUnitsMeanZFR_shuf2, numGoodUnits, numMiceProbes, numPcsToUse, pcaSamps, usePcContrast, alignPCs);

    savedEuclDistanceShuf = [];
    currentEndUnit = 0;

    for i = 1:numMiceProbes
        currentStartUnit = currentEndUnit + 1;
        currentEndUnit = currentEndUnit + numGoodUnits(i);
        if numGoodUnits(i) < numPcsToUse
            continue; % go to next mouse if this particular recording/region has less dimensions than we intend to use
        end

        % project mean shuffled activity for each mouse onto its respective PC matrix to get trajectories for each mouse & overlay
        shuf1MeanProj = pcSeg{i}'*allUnitsMeanZFR_shuf1(currentStartUnit:currentEndUnit,plotSamps(1):plotSamps(2)); %option to apply same PCs as above to shuffled means
        shuf2MeanProj = pcSeg{i}'*allUnitsMeanZFR_shuf2(currentStartUnit:currentEndUnit,plotSamps(1):plotSamps(2));
        %shuf1MeanProj = pcSegShuf{i}'*allUnitsMeanZFR_shuf1(currentStartUnit:currentEndUnit,plotSamps(1):plotSamps(2)); %option to apply new shuffled PCs to shuffled means
        %shuf2MeanProj = pcSegShuf{i}'*allUnitsMeanZFR_shuf2(currentStartUnit:currentEndUnit,plotSamps(1):plotSamps(2));
        shuf1MeanProj = smoothdata(shuf1MeanProj, 2, 'gaussian', smoothingSamples); %ex. smooth window 25 is 25*20ms PSTH bins = 500ms
        shuf2MeanProj = smoothdata(shuf2MeanProj, 2, 'gaussian', smoothingSamples);
    
        euclDistance = zeros(size(shuf1MeanProj,2)-1,1);
        for dataPoint = 1:size(shuf1MeanProj,2)-1
            euclDistance(dataPoint) = norm(shuf1MeanProj(:,dataPoint) - shuf2MeanProj(:,dataPoint));
        end
    
        % calculate Euclidean distance between trajectories over time, normalized by sqrt(# neurons recorded), since Euclidean distance scales by the sqrt of the number of dimensions
%         plot(timescalePlot(1:end-1), euclDistance./sqrt(numGoodUnits(i)))
        savedEuclDistanceShuf(i, :) = euclDistance./sqrt(numGoodUnits(i)); %save for bounded line below
    end %end loop over mice


    if ~isempty(savedEuclDistanceShuf)
        [meanEuclDistance, semEuclDistance] = grpstats(savedEuclDistanceShuf,[],{'mean' 'sem'}); %mean across columns/timepoints
        [hl2,~] = boundedline(timescalePlot(1:end-1), meanEuclDistance, semEuclDistance, 'cmap', [0 0 0], 'alpha','y', 'transparency', 0.5, 'LineWidth', 3);
        ylim([0 0.5])
        zeroXDottedLine;
        hold off;
        legend([hl1, hl2], {'avoid - react', 'shuffle'}, 'Box', 'off');
        xlabel('time (sec)', 'FontSize', fontSz)
        ylabel('distance between PC trajectories (z-Hz/neuron)', 'FontSize', fontSz)
        set(gcf,'color','w'); %set figure background white
        title(titleStr);
        makepretty;
    end

end

function [shufPcDist] = plotAvoidVsReactPcDistancesAllMiceEqualTrials(allUnitsMeanZFR, allUnitsAvoidMeanZFR, allUnitsReactMeanZFR, allUnitsMeanZFR_shuf1, allUnitsMeanZFR_shuf2, numGoodUnits, timescaleSeg, pcaWindow, plotWindow, titleStr, cmapAnat)
% input is all #repeats x N x T matrices that concatenate trial averages (over avoid / react trials) across mice/sessions
% also input numGoodUnits matrix that keeps track of # neuron splits between mice/sessions, and time windows for PCA calculation and plotting
% return shuffled PC distance control for later averaging across regions

    numMiceProbes = length(numGoodUnits);
    pcaSamps = [find(pcaWindow(1)>=timescaleSeg,1,'last') find(pcaWindow(2)>=timescaleSeg,1,'last')];
    plotSamps = [find(plotWindow(1)>=timescaleSeg,1,'last') find(plotWindow(2)>=timescaleSeg,1,'last')];
    timescalePlot = timescaleSeg(plotSamps(1):plotSamps(2));

    usePcContrast = 1; % whether we perform contrastive PCA on data1 vs data2, or vanilla PCA on allData 
    numPcsToUse = 10;
    alignPCs = 1; % whether to align all mice/sessions into common subspace for computing PCs; otherwise PCs are computed independently for each mouse/session (default for simplicity)
    smoothingSamples = [25 0]; %ex. smooth window 25 is 25*20ms PSTH bins = 500ms; for causal half-kernel use [X 0]
    %plotting constants
    fontSz = 16;
     
    savedEuclDistanceRpts = [];
    savedEuclDistanceShufRpts = [];
    numRepeats = 10; %should match number in other EqualTrials functions
    for rpt = 1:numRepeats 
        % for each sampling repeat, take get PCs and distances for each mouse
        [pcSeg, varExplained] = getEphysPCs(allUnitsMeanZFR, squeeze(allUnitsAvoidMeanZFR(rpt,:,:)), squeeze(allUnitsReactMeanZFR(rpt,:,:)), numGoodUnits, numMiceProbes, numPcsToUse, pcaSamps, usePcContrast, alignPCs); 
   
        currentEndUnit = 0;
        for i = 1:numMiceProbes %plot all probes separately
            currentStartUnit = currentEndUnit + 1;
            currentEndUnit = currentEndUnit + numGoodUnits(i);
            if numGoodUnits(i) < numPcsToUse
                continue; % go to next mouse if this particular recording/region has less dimensions than we intend to use
            end
            % REAL DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % project mean avoid/ react activity for each mouse onto its respective PC matrix to get trajectories for each mouse & overlay
            avoidMeanProj = pcSeg{i}'*squeeze(allUnitsAvoidMeanZFR(rpt,currentStartUnit:currentEndUnit,plotSamps(1):plotSamps(2))); %apply PCs to the mean responses across trials
            reactMeanProj = pcSeg{i}'*squeeze(allUnitsReactMeanZFR(rpt,currentStartUnit:currentEndUnit,plotSamps(1):plotSamps(2)));
            avoidMeanProj = smoothdata(avoidMeanProj, 2, 'gaussian', smoothingSamples); %ex. smooth window 25 is 25*20ms PSTH bins = 500ms
            reactMeanProj = smoothdata(reactMeanProj, 2, 'gaussian', smoothingSamples);
            euclDistance = zeros(size(avoidMeanProj,2)-1,1);
            for dataPoint = 1:size(avoidMeanProj,2)-1
                euclDistance(dataPoint) = norm(avoidMeanProj(:,dataPoint) - reactMeanProj(:,dataPoint));
            end
            % calculate Euclidean distance between trajectories over time, normalized by sqrt(# neurons recorded), since Euclidean distance scales by the sqrt of the number of dimensions
            savedEuclDistanceRpts(rpt, i, :) = euclDistance./sqrt(numGoodUnits(i)); %save for bounded line below
       
            % SHUFFLED DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % project mean shuffled activity for each mouse onto its respective PC matrix to get trajectories for each mouse & overlay
            shuf1MeanProj = pcSeg{i}'*squeeze(allUnitsMeanZFR_shuf1(rpt,currentStartUnit:currentEndUnit,plotSamps(1):plotSamps(2))); %option to apply same PCs as above to shuffled means
            shuf2MeanProj = pcSeg{i}'*squeeze(allUnitsMeanZFR_shuf2(rpt,currentStartUnit:currentEndUnit,plotSamps(1):plotSamps(2)));
            shuf1MeanProj = smoothdata(shuf1MeanProj, 2, 'gaussian', smoothingSamples); %ex. smooth window 25 is 25*20ms PSTH bins = 500ms
            shuf2MeanProj = smoothdata(shuf2MeanProj, 2, 'gaussian', smoothingSamples);
            euclDistanceShuf = zeros(size(shuf1MeanProj,2)-1,1);
            for dataPoint = 1:size(shuf1MeanProj,2)-1
                euclDistanceShuf(dataPoint) = norm(shuf1MeanProj(:,dataPoint) - shuf2MeanProj(:,dataPoint));
            end
            savedEuclDistanceShufRpts(rpt, i, :) = euclDistanceShuf./sqrt(numGoodUnits(i)); %save for bounded line below

        end %end loop over mice
    end %end sampling repeats

    % now take mean across repeats for each mouse:
    savedEuclDistance = squeeze(mean(savedEuclDistanceRpts,1));
    savedEuclDistanceShuf = squeeze(mean(savedEuclDistanceShufRpts,1));
    shufPcDist = savedEuclDistanceShuf;

    figure; hold on;
%     % optionally plot individual mice to see variability:
%     for i = 1:numMiceProbes
%         plot(timescalePlot(1:end-1), savedEuclDistance(i,:), 'Color', cmapAnat, 'LineWidth', 0.1)
%         plot(timescalePlot(1:end-1), savedEuclDistanceShuf(i,:), 'Color', [0 0 0], 'LineWidth', 0.1)
%     end
    [meanEuclDistance, semEuclDistance] = grpstats(savedEuclDistance,[],{'mean' 'sem'}); %mean across columns/timepoints
    [hl1,~] = boundedline(timescalePlot(1:end-1), meanEuclDistance, semEuclDistance, 'cmap', cmapAnat, 'alpha','y', 'transparency', 0.5, 'LineWidth', 3);
    [meanEuclDistance, semEuclDistance] = grpstats(savedEuclDistanceShuf,[],{'mean' 'sem'}); %mean across columns/timepoints
    [hl2,~] = boundedline(timescalePlot(1:end-1), meanEuclDistance, semEuclDistance, 'cmap', [0 0 0], 'alpha','y', 'transparency', 0.5, 'LineWidth', 3);
    ylim([0 0.7])
    xlim([timescalePlot(1)+0.5 timescalePlot(end)]) %adjust for smoothing edge effects
    zeroXDottedLine;
    hold off;
    legend([hl1, hl2], {'avoid - react', 'shuffle'}, 'Box', 'off');
    xlabel('time (sec)', 'FontSize', fontSz)
    ylabel('distance between PC trajectories (z-Hz/neuron)', 'FontSize', fontSz)
    title(titleStr);
    makepretty;

end

function [shufPcDist] = plotAvoidVsReactPcDistancesAllMiceEqualTrialsNullRemoved(allUnitsMeanZFR, allUnitsAvoidMeanZFR, allUnitsReactMeanZFR, allUnitsNullMeanZFR, allUnitsMeanZFR_shuf1, allUnitsMeanZFR_shuf2, numGoodUnits, timescaleSeg, pcaWindow, plotWindow, titleStr, cmapAnat)
% input is all #repeats x N x T matrices that concatenate trial averages (over avoid / react trials) across mice/sessions
% also input numGoodUnits matrix that keeps track of # neuron splits between mice/sessions, and time windows for PCA calculation and plotting
% return shuffled PC distance control for later averaging across regions

    numMiceProbes = length(numGoodUnits);
    pcaSamps = [find(pcaWindow(1)>=timescaleSeg,1,'last') find(pcaWindow(2)>=timescaleSeg,1,'last')];
    plotSamps = [find(plotWindow(1)>=timescaleSeg,1,'last') find(plotWindow(2)>=timescaleSeg,1,'last')];
    timescalePlot = timescaleSeg(plotSamps(1):plotSamps(2));
%     nullspaceWindow = [-0.5 0.2]; %for cold air
%     nullspaceWindow = [0 1]; %for cue
    nullspaceWindow = [-1.5 0]; %for VHEx
    nullSamps = [find(nullspaceWindow(1)>=timescaleSeg,1,'last') find(nullspaceWindow(2)>=timescaleSeg,1,'last')];

    usePcContrast = 1; % whether we perform contrastive PCA on data1 vs data2, or vanilla PCA on allData 
    numPcsToUse = 10;
    smoothingSamples = [25 0]; %ex. smooth window 25 is 25*20ms PSTH bins = 500ms; for causal half-kernel use [X 0]
    %plotting constants
    fontSz = 16;
     
    savedEuclDistanceRpts = [];
    savedEuclDistanceShufRpts = [];
    numRepeats = 10; %should match number in other EqualTrials functions
    for rpt = 1:numRepeats 
        % for each sampling repeat, take get PCs and distances for each mouse, after removing covariance patterns from nullspace (e.g. cold air on)
        [pcSeg, varExplained] = getEphysPCsWithNullspace(allUnitsMeanZFR, squeeze(allUnitsAvoidMeanZFR(rpt,:,:)), squeeze(allUnitsReactMeanZFR(rpt,:,:)), allUnitsNullMeanZFR, numGoodUnits, numMiceProbes, numPcsToUse, pcaSamps, nullSamps, usePcContrast); %

        currentEndUnit = 0;
        for i = 1:numMiceProbes %plot all probes separately
            currentStartUnit = currentEndUnit + 1;
            currentEndUnit = currentEndUnit + numGoodUnits(i);
            if numGoodUnits(i) < numPcsToUse
                continue; % go to next mouse if this particular recording/region has less dimensions than we intend to use
            end
            % REAL DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % project mean avoid/ react activity for each mouse onto its respective PC matrix to get trajectories for each mouse & overlay
            avoidMeanProj = pcSeg{i}'*squeeze(allUnitsAvoidMeanZFR(rpt,currentStartUnit:currentEndUnit,plotSamps(1):plotSamps(2))); %apply PCs to the mean responses across trials
            reactMeanProj = pcSeg{i}'*squeeze(allUnitsReactMeanZFR(rpt,currentStartUnit:currentEndUnit,plotSamps(1):plotSamps(2)));
            avoidMeanProj = smoothdata(avoidMeanProj, 2, 'gaussian', smoothingSamples); %ex. smooth window 25 is 25*20ms PSTH bins = 500ms
            reactMeanProj = smoothdata(reactMeanProj, 2, 'gaussian', smoothingSamples);
            euclDistance = zeros(size(avoidMeanProj,2)-1,1);
            for dataPoint = 1:size(avoidMeanProj,2)-1
                euclDistance(dataPoint) = norm(avoidMeanProj(:,dataPoint) - reactMeanProj(:,dataPoint));
            end
            % calculate Euclidean distance between trajectories over time, normalized by sqrt(# neurons recorded), since Euclidean distance scales by the sqrt of the number of dimensions
            savedEuclDistanceRpts(rpt, i, :) = euclDistance./sqrt(numGoodUnits(i)); %save for bounded line below
       
            % SHUFFLED DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % project mean shuffled activity for each mouse onto its respective PC matrix to get trajectories for each mouse & overlay
            shuf1MeanProj = pcSeg{i}'*squeeze(allUnitsMeanZFR_shuf1(rpt,currentStartUnit:currentEndUnit,plotSamps(1):plotSamps(2))); %option to apply same PCs as above to shuffled means
            shuf2MeanProj = pcSeg{i}'*squeeze(allUnitsMeanZFR_shuf2(rpt,currentStartUnit:currentEndUnit,plotSamps(1):plotSamps(2)));
            shuf1MeanProj = smoothdata(shuf1MeanProj, 2, 'gaussian', smoothingSamples); %ex. smooth window 25 is 25*20ms PSTH bins = 500ms
            shuf2MeanProj = smoothdata(shuf2MeanProj, 2, 'gaussian', smoothingSamples);
            euclDistanceShuf = zeros(size(shuf1MeanProj,2)-1,1);
            for dataPoint = 1:size(shuf1MeanProj,2)-1
                euclDistanceShuf(dataPoint) = norm(shuf1MeanProj(:,dataPoint) - shuf2MeanProj(:,dataPoint));
            end
            savedEuclDistanceShufRpts(rpt, i, :) = euclDistanceShuf./sqrt(numGoodUnits(i)); %save for bounded line below

        end %end loop over mice
    end %end sampling repeats

    % now take mean across repeats for each mouse:
    savedEuclDistance = squeeze(mean(savedEuclDistanceRpts,1));
    savedEuclDistanceShuf = squeeze(mean(savedEuclDistanceShufRpts,1));
    shufPcDist = savedEuclDistanceShuf;

    figure; hold on;
%     % optionally plot individual mice to see variability:
%     for i = 1:numMiceProbes
%         plot(timescalePlot(1:end-1), savedEuclDistance(i,:), 'Color', cmapAnat, 'LineWidth', 0.1)
%         plot(timescalePlot(1:end-1), savedEuclDistanceShuf(i,:), 'Color', [0 0 0], 'LineWidth', 0.1)
%     end
    [meanEuclDistance, semEuclDistance] = grpstats(savedEuclDistance,[],{'mean' 'sem'}); %mean across columns/timepoints
    [hl1,~] = boundedline(timescalePlot(1:end-1), meanEuclDistance, semEuclDistance, 'cmap', cmapAnat, 'alpha','y', 'transparency', 0.5, 'LineWidth', 3);
    [meanEuclDistance, semEuclDistance] = grpstats(savedEuclDistanceShuf,[],{'mean' 'sem'}); %mean across columns/timepoints
    [hl2,~] = boundedline(timescalePlot(1:end-1), meanEuclDistance, semEuclDistance, 'cmap', [0 0 0], 'alpha','y', 'transparency', 0.5, 'LineWidth', 3);
    ylim([0 0.7])
    xlim([timescalePlot(1)+0.5 timescalePlot(end)]) %adjust for smoothing edge effects
    zeroXDottedLine;
    hold off;
    legend([hl1, hl2], {'avoid - react', 'shuffle'}, 'Box', 'off');
    xlabel('time (sec)', 'FontSize', fontSz)
    ylabel('distance between PC trajectories (z-Hz/neuron)', 'FontSize', fontSz)
    title(titleStr);
    makepretty;

end

function plotAvoidVsReactPcDistancesFromSubspace(allDataSubspace, removeNull)
% input is subspace that already contains trials projected into PCA space, so just take mean across trials and then calculate Euclidean distance

    plotWindow = [-2 1]; % [-1 2] for cue, [-2 1] for vhex
    % make structure to iterate brain regions in a loop to plot all at once
    params(1).probeName = 'Probe1'; params(1).anatName = 'CTX'; params(1).labelName = 'PFC'; params(1).atCue = true; params(1).cmapIdx = 2;
    params(2).probeName = 'Probe1'; params(2).anatName = 'CTX'; params(2).labelName = 'PFC'; params(2).atCue = false; params(2).cmapIdx = 2;
    params(3).probeName = 'Probe0'; params(3).anatName = 'CTX'; params(3).labelName = 'M1'; params(3).atCue = true; params(3).cmapIdx = 1;
    params(4).probeName = 'Probe0'; params(4).anatName = 'CTX'; params(4).labelName = 'M1'; params(4).atCue = false; params(4).cmapIdx = 1;
    params(5).probeName = 'Probe0'; params(5).anatName = 'BS'; params(5).labelName = 'TH/HY'; params(5).atCue = true; params(5).cmapIdx = 3;
    params(6).probeName = 'Probe0'; params(6).anatName = 'BS'; params(6).labelName = 'TH/HY'; params(6).atCue = false; params(6).cmapIdx = 3;
    params(7).probeName = 'Probe2'; params(7).anatName = 'HB'; params(7).labelName = 'HB'; params(7).atCue = true; params(7).cmapIdx = 4;
    params(8).probeName = 'Probe2'; params(8).anatName = 'HB'; params(8).labelName = 'HB'; params(8).atCue = false; params(8).cmapIdx = 4;

    f = figure; tiledlayout(4,2,"TileSpacing","compact");
    for p = 1:length(params)
    
        atCueOrVhex = params(p).atCue;
        probeName = params(p).probeName;
        anatName = params(p).anatName;
        labelName = params(p).labelName;
        cmapIdx = params(p).cmapIdx;
    
        if atCueOrVhex %true for cue timing
            if removeNull
                condAvoid = 'allUnitsTrialsAvoidAtCueNullMean'; %#ok<*NASGU> 
                condReact = 'allUnitsTrialsReactAtCueNullMean'; 
                condShuf1 = 'allUnitsTrialsShufAtCue1NullMean';
                condShuf2 = 'allUnitsTrialsShufAtCue2NullMean';
                titleStr = [labelName '(-null) at cue'];
            else
                condAvoid = 'allUnitsTrialsAvoidAtCueMean'; 
                condReact = 'allUnitsTrialsReactAtCueMean'; 
                condShuf1 = 'allUnitsTrialsShufAtCue1Mean';
                condShuf2 = 'allUnitsTrialsShufAtCue2Mean';
                titleStr = [labelName ' at cue'];
            end
        else %else @ VHEx timing
            if removeNull
                condAvoid = 'allUnitsTrialsAvoidNullMean';
                condReact = 'allUnitsTrialsReactNullMean';
                condShuf1 = 'allUnitsTrialsShuf1NullMean';
                condShuf2 = 'allUnitsTrialsShuf2NullMean';
                titleStr = [labelName '(-null) at VHEx'];
            else
                condAvoid = 'allUnitsTrialsAvoidMean';
                condReact = 'allUnitsTrialsReactMean';
                condShuf1 = 'allUnitsTrialsShuf1Mean';
                condShuf2 = 'allUnitsTrialsShuf2Mean';
                titleStr = [labelName ' at VHEx'];
            end
        end
    
        numMiceProbes = size(allDataSubspace.(probeName).(anatName).(condAvoid),2);
    %     numDimensions = size(allDataSubspace.(probeName).(anatName).(condAvoid){1}{1},1); %subspace dimensions for normalization below
        plotSamps = [find(plotWindow(1)>=allDataSubspace.timescaleSeg,1,'last') find(plotWindow(2)>=allDataSubspace.timescaleSeg,1,'last')];
        timescalePlot = allDataSubspace.timescaleSeg(plotSamps(1):plotSamps(2));
    
        %other constants
        fontSz = 16;
        cmapLine = allDataSubspace.lineCmaps(cmapIdx,:);
        numRepeats = 10;
        smoothingSamples = [10 0]; %ex. smooth window 10 is 10*20ms PSTH bins = 200ms; for causal half-kernel use [X 0]
    
        savedEuclDistanceRpts = [];
        savedEuclDistanceRptsShuf = [];
    
        dimStart = 1; %option for using a subset of the subspace dimensions
        dimStop = 10; %option for using a subset of the subspace dimensions
        numDimensions = dimStop - dimStart + 1;
        % for each mouse/repeat, compute Euclidean distance from subspace trial means, then average over trial sampling repeats
        for i = 1:numMiceProbes 
%         for i = 1:5 % local opto option
%         for i = 6:10 % large scale opto option
            thisMouseAvoidData = allDataSubspace.(probeName).(anatName).(condAvoid){i}; %unpack to get per-mouse cells that include all sample repeats
            thisMouseReactData = allDataSubspace.(probeName).(anatName).(condReact){i};
            thisMouseShuf1Data = allDataSubspace.(probeName).(anatName).(condShuf1){i};
            thisMouseShuf2Data = allDataSubspace.(probeName).(anatName).(condShuf2){i};
            for rpt = 1:numRepeats
                thisRptAvoidData = thisMouseAvoidData{rpt}(dimStart:dimStop,plotSamps(1):plotSamps(2)); %unpack one more time to get each repeat cell - now this cell for each trial mean
                thisRptReactData = thisMouseReactData{rpt}(dimStart:dimStop,plotSamps(1):plotSamps(2));
                thisRptShuf1Data = thisMouseShuf1Data{rpt}(dimStart:dimStop,plotSamps(1):plotSamps(2));
                thisRptShuf2Data = thisMouseShuf2Data{rpt}(dimStart:dimStop,plotSamps(1):plotSamps(2));
                thisRptAvoidData = smoothdata(thisRptAvoidData, 2, 'gaussian', smoothingSamples); 
                thisRptReactData = smoothdata(thisRptReactData, 2, 'gaussian', smoothingSamples); 
                thisRptShuf1Data = smoothdata(thisRptShuf1Data, 2, 'gaussian', smoothingSamples); 
                thisRptShuf2Data = smoothdata(thisRptShuf2Data, 2, 'gaussian', smoothingSamples); 
                euclDistance = zeros(size(thisRptAvoidData,2)-1,1);
                euclDistanceShuf = zeros(size(thisRptAvoidData,2)-1,1);
    %             figure; imagesc(thisRptAvoidData); clim([-1 1]);
    %             figure; imagesc(thisRptReactData); clim([-1 1]);
                for dataPoint = 1:size(thisRptAvoidData,2)-1
                    euclDistance(dataPoint) = norm(thisRptAvoidData(:,dataPoint) - thisRptReactData(:,dataPoint));
                    euclDistanceShuf(dataPoint) = norm(thisRptShuf1Data(:,dataPoint) - thisRptShuf2Data(:,dataPoint));
                end
                % calculate Euclidean distance between trajectories over time (normalized by sqrt(# subspace dimensions) for equal comparison with or without null removed)
                savedEuclDistanceRpts(rpt, i, :) = euclDistance./sqrt(numDimensions); %save for bounded line below
                savedEuclDistanceRptsShuf(rpt, i, :) = euclDistanceShuf./sqrt(numDimensions); 
            end %end sampling repeats
        end %end numMice
    
        % now take mean across repeats for each mouse:
        savedEuclDistance = squeeze(mean(savedEuclDistanceRpts,1));
        savedEuclDistanceShuf = squeeze(mean(savedEuclDistanceRptsShuf,1));
    
       nexttile; hold on;
        % optionally plot individual mice to see variability:
%         for i = 1:numMiceProbes
% %             plot(timescalePlot(1:end-1), savedEuclDistance(i,:), 'Color', cmapLine, 'LineWidth', 0.1)
%             plot(timescalePlot(1:end-1), savedEuclDistanceShuf(i,:), 'Color', [0 0 0], 'LineWidth', 0.1)
%         end
        [meanEuclDistance, semEuclDistance] = grpstats(savedEuclDistance,[],{'mean' 'sem'}); %mean across columns/timepoints
        [h1(p),~] = boundedline(timescalePlot(1:end-1), meanEuclDistance, semEuclDistance, 'cmap', cmapLine, 'alpha','y', 'transparency', 0.5, 'LineWidth', 3);
        [meanEuclDistanceShuf, semEuclDistanceShuf] = grpstats(savedEuclDistanceShuf,[],{'mean' 'sem'}); %mean across columns/timepoints
        [h2,~] = boundedline(timescalePlot(1:end-1), meanEuclDistanceShuf, semEuclDistanceShuf, 'cmap', [0 0 0], 'alpha','y', 'transparency', 0.5, 'LineWidth', 3);
        ylim([0 1.0]) %2.0 for full subspace, 1.0 for nullRemoved
%         xlim([timescalePlot(1)+0.5 timescalePlot(end)]) %adjust for smoothing edge effects
        zeroXDottedLine;
        makepretty;
        hold off;
        if p==1
            title('at cue')
        elseif p==2
            title('at VHEx')
        elseif p==5
            ylabel('distance between PC trajectories', 'FontSize', fontSz)
        elseif p==8
            xlabel('time (sec)', 'FontSize', fontSz)
        end

    end %end params loop

    f.Position(3:4) = [1200 900]; %set figure size appropriate for plots
%     legend([h1(1), h1(3), h1(5), h1(7), h2], {params(1).labelName, params(3).labelName, params(5).labelName, params(7).labelName,'shuffle'}, 'Box', 'off', 'Location','northeast');

end

function plotAvoidVsReactPcTrajectoriesFromSubspace(allDataSubspace, removeNull)
% input is subspace that already contains trials projected into PCA space, so just take mean across trials and then calculate Euclidean distance

    plotWindow = [-2 2];
    % make structure to iterate brain regions in a loop to plot all at once
%     params(1).probeName = 'Probe1'; params(1).anatName = 'CTX'; params(1).labelName = 'PFC'; params(1).atCue = true; params(1).cmapIdx = 2;
%     params(2).probeName = 'Probe1'; params(2).anatName = 'CTX'; params(2).labelName = 'PFC'; params(2).atCue = false; params(2).cmapIdx = 2;
%     params(3).probeName = 'Probe0'; params(3).anatName = 'CTX'; params(3).labelName = 'M1'; params(3).atCue = true; params(3).cmapIdx = 1;
%     params(4).probeName = 'Probe0'; params(4).anatName = 'CTX'; params(4).labelName = 'M1'; params(4).atCue = false; params(4).cmapIdx = 1;
%     params(5).probeName = 'Probe0'; params(5).anatName = 'BS'; params(5).labelName = 'TH/HY'; params(5).atCue = true; params(5).cmapIdx = 3;
%     params(6).probeName = 'Probe0'; params(6).anatName = 'BS'; params(6).labelName = 'TH/HY'; params(6).atCue = false; params(6).cmapIdx = 3;
%     params(7).probeName = 'Probe2'; params(7).anatName = 'HB'; params(7).labelName = 'HB'; params(7).atCue = true; params(7).cmapIdx = 4;
%     params(8).probeName = 'Probe2'; params(8).anatName = 'HB'; params(8).labelName = 'HB'; params(8).atCue = false; params(8).cmapIdx = 4;

    params(1).probeName = 'Probe1'; params(1).anatName = 'CTX'; params(1).labelName = 'PFC'; params(1).atCue = false; params(1).cmapIdx = 2;
    params(2).probeName = 'Probe2'; params(2).anatName = 'HB'; params(2).labelName = 'HB'; params(2).atCue = false; params(2).cmapIdx = 4;

    plot3d = 1;
    pcaIdx1 = 1; %make this variable to easy to change which PCs to view if necessary
    pcaIdx2 = 2;
    pcaIdx3 = 3;

    %plotting constants
    numColorPts = 256; %standard colormap pts
    R = linspace(0.3,0.9,numColorPts)'; %go from half-black to full color
    G = linspace(0.1,0.5,numColorPts)';
    B = zeros(numColorPts,1);
    cmap1 = [R G B];
    R = zeros(numColorPts,1);
    G = linspace(0.0,0.4,numColorPts)';
    B = linspace(0.3,0.9,numColorPts)';
    cmap2 = [R G B];
    cmapIdx = numColorPts - round(numColorPts/4); %for single colors for individual PC plots
    psub1 = [];
    psub2 = [];

    for p = 1:length(params)
    
        % 1 plot for each region/condition
        if plot3d
            pc3dFigure(p) = figure; 
            tiledlayout('flow');
        end

        atCueOrVhex = params(p).atCue;
        probeName = params(p).probeName;
        anatName = params(p).anatName;
        labelName = params(p).labelName;
        cmapIdx = params(p).cmapIdx;
    
        if atCueOrVhex %true for cue timing
            if removeNull
                condAvoid = 'allUnitsTrialsAvoidAtCueNullMean'; %#ok<*NASGU> 
                condReact = 'allUnitsTrialsReactAtCueNullMean'; 
                condShuf1 = 'allUnitsTrialsShufAtCue1NullMean';
                condShuf2 = 'allUnitsTrialsShufAtCue2NullMean';
                titleStr = [labelName '(-null) at cue'];
            else
                condAvoid = 'allUnitsTrialsAvoidAtCueMean'; 
                condReact = 'allUnitsTrialsReactAtCueMean'; 
                condShuf1 = 'allUnitsTrialsShufAtCue1Mean';
                condShuf2 = 'allUnitsTrialsShufAtCue2Mean';
                titleStr = [labelName ' at cue'];
            end
        else %else @ VHEx timing
            if removeNull
                condAvoid = 'allUnitsTrialsAvoidNullMean';
                condReact = 'allUnitsTrialsReactNullMean';
                condShuf1 = 'allUnitsTrialsShuf1NullMean';
                condShuf2 = 'allUnitsTrialsShuf2NullMean';
                titleStr = [labelName '(-null) at VHEx'];
            else
                condAvoid = 'allUnitsTrialsAvoidMean';
                condReact = 'allUnitsTrialsReactMean';
                condShuf1 = 'allUnitsTrialsShuf1Mean';
                condShuf2 = 'allUnitsTrialsShuf2Mean';
                titleStr = [labelName ' at VHEx'];
            end
        end
    
        numMiceProbes = size(allDataSubspace.(probeName).(anatName).(condAvoid),2);
        plotSamps = [find(plotWindow(1)>=allDataSubspace.timescaleSeg,1,'last') find(plotWindow(2)>=allDataSubspace.timescaleSeg,1,'last')];
        timescalePlot = allDataSubspace.timescaleSeg(plotSamps(1):plotSamps(2));
    
        %other constants
        fontSz = 16;
        numRepeats = 1; %only use 1 for example plot
        smoothingSamples = [25 0]; %ex. smooth window 10 is 10*20ms PSTH bins = 200ms; for causal half-kernel use [X 0]
    
    
        dimStart = 1; %option for using a subset of the subspace dimensions
        dimStop = 10; %option for using a subset of the subspace dimensions
        numDimensions = dimStop - dimStart + 1;

        for i = 1:numMiceProbes 
            thisMouseAvoidData = allDataSubspace.(probeName).(anatName).(condAvoid){i}; %unpack to get per-mouse cells that include all sample repeats
            thisMouseReactData = allDataSubspace.(probeName).(anatName).(condReact){i};
            thisMouseShuf1Data = allDataSubspace.(probeName).(anatName).(condShuf1){i};
            thisMouseShuf2Data = allDataSubspace.(probeName).(anatName).(condShuf2){i};
            for rpt = 1:numRepeats
                thisRptAvoidData = thisMouseAvoidData{rpt}(dimStart:dimStop,plotSamps(1):plotSamps(2)); %unpack one more time to get each repeat cell - now this cell for each trial mean
                thisRptReactData = thisMouseReactData{rpt}(dimStart:dimStop,plotSamps(1):plotSamps(2));
                thisRptShuf1Data = thisMouseShuf1Data{rpt}(dimStart:dimStop,plotSamps(1):plotSamps(2));
                thisRptShuf2Data = thisMouseShuf2Data{rpt}(dimStart:dimStop,plotSamps(1):plotSamps(2));
                thisRptAvoidData = smoothdata(thisRptAvoidData, 2, 'gaussian', smoothingSamples); 
                thisRptReactData = smoothdata(thisRptReactData, 2, 'gaussian', smoothingSamples); 
                thisRptShuf1Data = smoothdata(thisRptShuf1Data, 2, 'gaussian', smoothingSamples); 
                thisRptShuf2Data = smoothdata(thisRptShuf2Data, 2, 'gaussian', smoothingSamples); 
 
                colorsIdx = floor(linspace(1, numColorPts, size(thisRptAvoidData,2))); %evenly spaced colors that show time dimension
    
                % loop for plotting PC trajectories in 3D
                if plot3d
                    figure(pc3dFigure(p))
                    nexttile; 
                    hold on;
                    zeroSample = find(timescalePlot>0,1,'first'); %for plotting midpoint
                    for dataPoint = 1:size(thisRptAvoidData,2)-1
            
                        %3D option
%                         if dataPoint == size(thisRptAvoidData,2)-1
%                             psub1 = plot3(thisRptAvoidData(pcaIdx1,dataPoint:dataPoint+1),thisRptAvoidData(pcaIdx2,dataPoint:dataPoint+1),thisRptAvoidData(pcaIdx3,dataPoint:dataPoint+1),'Color',cmap1(colorsIdx(dataPoint),:), 'LineWidth', 3, 'DisplayName','avoid'); %subset for legend display name
%                             psub2 = plot3(thisRptReactData(pcaIdx1,dataPoint:dataPoint+1),thisRptReactData(pcaIdx2,dataPoint:dataPoint+1),thisRptReactData(pcaIdx3,dataPoint:dataPoint+1),'Color',cmap2(colorsIdx(dataPoint),:),  'LineWidth', 3, 'DisplayName','react'); %subset for legend display name
%                         else
%                             plot3(thisRptAvoidData(pcaIdx1,dataPoint:dataPoint+1),thisRptAvoidData(pcaIdx2,dataPoint:dataPoint+1),thisRptAvoidData(pcaIdx3,dataPoint:dataPoint+1),'Color',cmap1(colorsIdx(dataPoint),:), 'LineWidth', 3); 
%                             plot3(thisRptReactData(pcaIdx1,dataPoint:dataPoint+1),thisRptReactData(pcaIdx2,dataPoint:dataPoint+1),thisRptReactData(pcaIdx3,dataPoint:dataPoint+1),'Color',cmap2(colorsIdx(dataPoint),:),  'LineWidth', 3); 
%                             if dataPoint == zeroSample %plot marker at t=0 extension threshold
%                                 plot3(thisRptAvoidData(1,dataPoint),thisRptAvoidData(2,dataPoint),thisRptAvoidData(3,dataPoint), 'Marker', 'o', 'MarkerFaceColor','m', 'MarkerSize', 9, 'MarkerEdgeColor','none'); 
%                                 plot3(thisRptReactData(1,dataPoint),thisRptReactData(2,dataPoint),thisRptReactData(3,dataPoint),'Marker', 'o','MarkerFaceColor','m', 'MarkerSize', 9, 'MarkerEdgeColor','none'); 
%                             end
%                             if dataPoint == 1 %plot marker beginning
%                                 plot3(thisRptAvoidData(1,dataPoint),thisRptAvoidData(2,dataPoint),thisRptAvoidData(3,dataPoint), 'Marker', 'o', 'MarkerFaceColor','g', 'MarkerSize', 9, 'MarkerEdgeColor','none'); 
%                                 plot3(thisRptReactData(1,dataPoint),thisRptReactData(2,dataPoint),thisRptReactData(3,dataPoint),'Marker', 'o','MarkerFaceColor','g', 'MarkerSize', 9, 'MarkerEdgeColor','none'); 
%                             end
%                         end
                        
                        % simpler 2d plot option
                        if dataPoint == size(thisRptAvoidData,2)-1
                            psub1 = plot(thisRptAvoidData(pcaIdx1,dataPoint:dataPoint+1),thisRptAvoidData(pcaIdx2,dataPoint:dataPoint+1),'Color',cmap1(colorsIdx(dataPoint),:), 'LineWidth', 3, 'DisplayName','avoid'); %subset for legend display name
                            psub2 = plot(thisRptReactData(pcaIdx1,dataPoint:dataPoint+1),thisRptReactData(pcaIdx2,dataPoint:dataPoint+1),'Color',cmap2(colorsIdx(dataPoint),:),  'LineWidth', 3, 'DisplayName','react'); %subset for legend display name
                        else
                            plot(thisRptAvoidData(pcaIdx1,dataPoint:dataPoint+1),thisRptAvoidData(pcaIdx2,dataPoint:dataPoint+1),'Color',cmap1(colorsIdx(dataPoint),:), 'LineWidth', 3); 
                            plot(thisRptReactData(pcaIdx1,dataPoint:dataPoint+1),thisRptReactData(pcaIdx2,dataPoint:dataPoint+1),'Color',cmap2(colorsIdx(dataPoint),:),  'LineWidth', 3); 
                            if dataPoint == zeroSample %plot marker at t=0 extension threshold
                                plot(thisRptAvoidData(1,dataPoint),thisRptAvoidData(2,dataPoint), 'Marker', 'o', 'MarkerFaceColor','m', 'MarkerSize', 9, 'MarkerEdgeColor','none'); 
                                plot(thisRptReactData(1,dataPoint),thisRptReactData(2,dataPoint),'Marker', 'o','MarkerFaceColor','m', 'MarkerSize', 9, 'MarkerEdgeColor','none'); 
                            end
                            if dataPoint == 1 %plot marker beginning
                                plot(thisRptAvoidData(1,dataPoint),thisRptAvoidData(2,dataPoint), 'Marker', 'o', 'MarkerFaceColor','g', 'MarkerSize', 9, 'MarkerEdgeColor','none'); 
                                plot(thisRptReactData(1,dataPoint),thisRptReactData(2,dataPoint),'Marker', 'o','MarkerFaceColor','g', 'MarkerSize', 9, 'MarkerEdgeColor','none'); 
                            end
                        end
        
        
                    end % end dataPoint loop
                    hold off;
        %             view(45,30);
                    axis tight;
                end % end if(plot3d)


            end %end sampling repeats
        end %end numMice

        pc3dFigure(p).Position(3:4) = [1200 900]; %set figure size appropriate for plots
        title(params(p).labelName)
    end %end params loop

    






end

function plotAvoidVsReactPcDistancesWithOptoFromSubspace(allDataSubspace, removeNull)
% input is subspace that already contains trials projected into PCA space, so just take mean across trials and then calculate Euclidean distance

    fontSz = 20;
    blueCmap = [0 0.4470 0.7410];
    orangeCmap = [0.8500 0.3250 0.0980];
    optoCmap = [0 0 0];
    plotWindow = [-1 2]; % [-1 2] for cue, [-2 1] for vhex
    % make structure to iterate brain regions in a loop to plot all at once
    params(1).probeName = 'Probe1'; params(1).anatName = 'CTX'; params(1).labelName = 'PFC'; params(1).atCue = true; params(1).cmapIdx = 2;
    params(2).probeName = 'Probe1'; params(2).anatName = 'CTX'; params(2).labelName = 'PFC'; params(2).atCue = false; params(2).cmapIdx = 2;
    params(3).probeName = 'Probe0'; params(3).anatName = 'CTX'; params(3).labelName = 'M1'; params(3).atCue = true; params(3).cmapIdx = 1;
    params(4).probeName = 'Probe0'; params(4).anatName = 'CTX'; params(4).labelName = 'M1'; params(4).atCue = false; params(4).cmapIdx = 1;
    params(5).probeName = 'Probe0'; params(5).anatName = 'BS'; params(5).labelName = 'TH/HY'; params(5).atCue = true; params(5).cmapIdx = 3;
    params(6).probeName = 'Probe0'; params(6).anatName = 'BS'; params(6).labelName = 'TH/HY'; params(6).atCue = false; params(6).cmapIdx = 3;
    params(7).probeName = 'Probe2'; params(7).anatName = 'HB'; params(7).labelName = 'HB'; params(7).atCue = true; params(7).cmapIdx = 4;
    params(8).probeName = 'Probe2'; params(8).anatName = 'HB'; params(8).labelName = 'HB'; params(8).atCue = false; params(8).cmapIdx = 4;

    f = figure; tiledlayout(4,2,"TileSpacing","compact");
    for p = 1:length(params)
    
        atCueOrVhex = params(p).atCue;
        probeName = params(p).probeName;
        anatName = params(p).anatName;
        labelName = params(p).labelName;
        cmapIdx = params(p).cmapIdx;
    
        if atCueOrVhex %true for cue timing
            if removeNull
                condAvoid = 'allUnitsTrialsAvoidAtCueNullMean'; %#ok<*NASGU> 
                condReact = 'allUnitsTrialsReactAtCueNullMean'; 
                condOpto = 'allUnitsTrialsOptoAtCueNullMean';
                condShuf1 = 'allUnitsTrialsShufAtCue1NullMean';
                titleStr = [labelName '(-null) at cue'];
            else
                condAvoid = 'allUnitsTrialsAvoidAtCueMean'; 
                condReact = 'allUnitsTrialsReactAtCueMean'; 
                condOpto = 'allUnitsTrialsOptoAtCueMean';
                condShuf1 = 'allUnitsTrialsShufAtCue1Mean';
                titleStr = [labelName ' at cue'];
            end
        else %else @ VHEx timing
            if removeNull
                condAvoid = 'allUnitsTrialsAvoidNullMean';
                condReact = 'allUnitsTrialsReactNullMean';
                condOpto = 'allUnitsTrialsOptoNullMean';
                condShuf1 = 'allUnitsTrialsShuf1NullMean';
                titleStr = [labelName '(-null) at VHEx'];
            else
                condAvoid = 'allUnitsTrialsAvoidMean';
                condReact = 'allUnitsTrialsReactMean';
                condOpto = 'allUnitsTrialsOptoMean';
                condShuf1 = 'allUnitsTrialsShuf1Mean';
                titleStr = [labelName ' at VHEx'];
            end
        end
    
        numMiceProbes = size(allDataSubspace.(probeName).(anatName).(condAvoid),2);
    %     numDimensions = size(allDataSubspace.(probeName).(anatName).(condAvoid){1}{1},1); %subspace dimensions for normalization below
        plotSamps = [find(plotWindow(1)>=allDataSubspace.timescaleSeg,1,'last') find(plotWindow(2)>=allDataSubspace.timescaleSeg,1,'last')];
        timescalePlot = allDataSubspace.timescaleSeg(plotSamps(1):plotSamps(2));
    
        %other constants
        fontSz = 16;
        cmapLine = allDataSubspace.lineCmaps(cmapIdx,:);
        numRepeats = 10;
        smoothingSamples = [10 0]; %ex. smooth window 25 is 25*20ms PSTH bins = 500ms; for causal half-kernel use [X 0]

        savedOptoAvoidEuclDistanceRpts = [];
        savedOptoReactEuclDistanceRpts = [];
        savedOptoShuffleEuclDistanceRpts = [];
    
        dimStart = 1; %option for using a subset of the subspace dimensions
        dimStop = 10; %option for using a subset of the subspace dimensions
        numDimensions = dimStop - dimStart + 1;
        % for each mouse/repeat, compute Euclidean distance from subspace trial means, then average over trial sampling repeats
        for i = 1:numMiceProbes
%         for i = 1:5 % local opto option
%         for i = 6:10 % large scale opto option 
            thisMouseAvoidData = allDataSubspace.(probeName).(anatName).(condAvoid){i}; %unpack to get per-mouse cells that include all sample repeats
            thisMouseReactData = allDataSubspace.(probeName).(anatName).(condReact){i};
            thisMouseOptoData = allDataSubspace.(probeName).(anatName).(condOpto){i};
            thisMouseShuf1Data = allDataSubspace.(probeName).(anatName).(condShuf1){i};
            for rpt = 1:numRepeats
                thisRptAvoidData = thisMouseAvoidData{rpt}(dimStart:dimStop,plotSamps(1):plotSamps(2)); %unpack one more time to get each repeat cell - now this cell for each trial mean
                thisRptReactData = thisMouseReactData{rpt}(dimStart:dimStop,plotSamps(1):plotSamps(2));
                thisRptOptoData = thisMouseOptoData{rpt}(dimStart:dimStop,plotSamps(1):plotSamps(2));
                thisRptShuf1Data = thisMouseShuf1Data{rpt}(dimStart:dimStop,plotSamps(1):plotSamps(2));

                thisRptAvoidData = smoothdata(thisRptAvoidData, 2, 'gaussian', smoothingSamples); 
                thisRptReactData = smoothdata(thisRptReactData, 2, 'gaussian', smoothingSamples); 
                thisRptOptoData = smoothdata(thisRptOptoData, 2, 'gaussian', smoothingSamples); 
                thisRptShuf1Data = smoothdata(thisRptShuf1Data, 2, 'gaussian', smoothingSamples); 

                euclDistanceOptoAvoid = zeros(size(thisRptAvoidData,2)-1,1);
                euclDistanceOptoReact = zeros(size(thisRptAvoidData,2)-1,1);
                euclDistanceOptoShuffle = zeros(size(thisRptAvoidData,2)-1,1);

                for dataPoint = 1:size(thisRptAvoidData,2)-1
                    euclDistanceOptoAvoid(dataPoint) = norm(thisRptOptoData(:,dataPoint) - thisRptAvoidData(:,dataPoint));
                    euclDistanceOptoReact(dataPoint) = norm(thisRptOptoData(:,dataPoint) - thisRptReactData(:,dataPoint));
                    euclDistanceOptoShuffle(dataPoint) = norm(thisRptOptoData(:,dataPoint) - thisRptShuf1Data(:,dataPoint));
                end

                % calculate Euclidean distance between trajectories over time (normalized by sqrt(# subspace dimensions) for equal comparison with or without null removed)
                savedOptoAvoidEuclDistanceRpts(rpt, i, :) = euclDistanceOptoAvoid./sqrt(numDimensions); %save for bounded line below
                savedOptoReactEuclDistanceRpts(rpt, i, :) = euclDistanceOptoReact./sqrt(numDimensions);
                savedOptoShuffleEuclDistanceRpts(rpt, i, :) = euclDistanceOptoShuffle./sqrt(numDimensions);
            end %end sampling repeats
        end %end numMice

        % now take mean across repeats for each mouse:
        savedOptoAvoidEuclDistance = squeeze(mean(savedOptoAvoidEuclDistanceRpts,1));
        savedOptoReactEuclDistance = squeeze(mean(savedOptoReactEuclDistanceRpts,1));
        savedOptoShuffleEuclDistance = squeeze(mean(savedOptoShuffleEuclDistanceRpts,1));
    
       nexttile; hold on;
    %     % optionally plot individual mice to see variability:
    %     for i = 1:numMiceProbes
    %         plot(timescalePlot(1:end-1), savedOptoAvoidEuclDistance(i,:), 'Color', orangeCmap, 'LineWidth', 0.1)
    %         plot(timescalePlot(1:end-1), savedOptoReactEuclDistance(i,:), 'Color', blueCmap, 'LineWidth', 0.1)
    % %         plot(timescalePlot(1:end-1), savedOptoShuffleEuclDistance(i,:), 'Color', [0 0 0], 'LineWidth', 0.1)
    %     end
        [meanOptoAvoidEuclDistance, semOptoAvoidEuclDistance] = grpstats(savedOptoAvoidEuclDistance,[],{'mean' 'sem'}); %mean across columns/timepoints
        [h1(p),~] = boundedline(timescalePlot(1:end-1), meanOptoAvoidEuclDistance, semOptoAvoidEuclDistance, 'cmap', orangeCmap, 'alpha','y', 'transparency', 0.5, 'LineWidth', 3);
        [meanOptoReactEuclDistance, semOptoReactEuclDistance] = grpstats(savedOptoReactEuclDistance,[],{'mean' 'sem'}); %mean across columns/timepoints
        [h2(p),~] = boundedline(timescalePlot(1:end-1), meanOptoReactEuclDistance, semOptoReactEuclDistance, 'cmap', blueCmap, 'alpha','y', 'transparency', 0.5, 'LineWidth', 3);
    %     [meanOptoShuffleEuclDistance, semOptoShuffleEuclDistance] = grpstats(savedOptoShuffleEuclDistance,[],{'mean' 'sem'}); %mean across columns/timepoints
    %     [h3(p),~] = boundedline(timescalePlot(1:end-1), meanOptoShuffleEuclDistance, semOptoShuffleEuclDistance, 'cmap', [0 0 0], 'alpha','y', 'transparency', 0.5);

        ylim([0 2.5])
%         xlim([timescalePlot(1)+0.5 timescalePlot(end)]) %adjust for smoothing edge effects
        zeroXDottedLine;
        title(titleStr);
        makepretty;
        hold off;
        if p==5
            ylabel('distance between PC trajectories', 'FontSize', fontSz)
        elseif p==8
            xlabel('time (sec)', 'FontSize', fontSz)
        end
    end %end params loop

    f.Position(3:4) = [1200 900]; %set figure size appropriate for plots
%     legend([h1(1) h2(1)], {'opto - avoid', 'opto - react'}, 'Box', 'off');
end

function plotAvoidVsReactPcDistancesFromStruct(allData, removeNull)
% input is subspace that already contains trials projected into PCA space, so just take mean across trials and then calculate Euclidean distance

    % make structure to iterate brain regions in a loop to plot all at once
    params(1).probeName = 'Probe1'; params(1).anatName = 'CTX'; params(1).labelName = 'PFC'; params(1).atCue = true; params(1).cmapIdx = 2;
    params(2).probeName = 'Probe1'; params(2).anatName = 'CTX'; params(2).labelName = 'PFC'; params(2).atCue = false; params(2).cmapIdx = 2;
    params(3).probeName = 'Probe0'; params(3).anatName = 'CTX'; params(3).labelName = 'M1'; params(3).atCue = true; params(3).cmapIdx = 1;
    params(4).probeName = 'Probe0'; params(4).anatName = 'CTX'; params(4).labelName = 'M1'; params(4).atCue = false; params(4).cmapIdx = 1;
    params(5).probeName = 'Probe0'; params(5).anatName = 'BS'; params(5).labelName = 'TH/HY'; params(5).atCue = true; params(5).cmapIdx = 3;
    params(6).probeName = 'Probe0'; params(6).anatName = 'BS'; params(6).labelName = 'TH/HY'; params(6).atCue = false; params(6).cmapIdx = 3;
    params(7).probeName = 'Probe2'; params(7).anatName = 'HB'; params(7).labelName = 'HB'; params(7).atCue = true; params(7).cmapIdx = 4;
    params(8).probeName = 'Probe2'; params(8).anatName = 'HB'; params(8).labelName = 'HB'; params(8).atCue = false; params(8).cmapIdx = 4;

    f = figure; tiledlayout(4,2,"TileSpacing","compact");
    for p = 1:length(params)
    
        atCueOrVhex = params(p).atCue;
        probeName = params(p).probeName;
        anatName = params(p).anatName;
        labelName = params(p).labelName;
        cmapIdx = params(p).cmapIdx;
    
        if atCueOrVhex %true for cue timing
            if removeNull
                condAvoid = 'allUnitsTrialsAvoidAtCueNullMean'; %#ok<*NASGU> 
                condReact = 'allUnitsTrialsReactAtCueNullMean'; 
                condShuf1 = 'allUnitsTrialsShufAtCue1NullMean';
                condShuf2 = 'allUnitsTrialsShufAtCue2NullMean';
                titleStr = [labelName '(-null) at cue'];
            else
                condAvoid = 'allUnitsTrialsAvoidAtCueMean'; 
                condReact = 'allUnitsTrialsReactAtCueMean'; 
                condShuf1 = 'allUnitsTrialsShufAtCue1Mean';
                condShuf2 = 'allUnitsTrialsShufAtCue2Mean';
                titleStr = [labelName ' at cue'];
            end
            pcaWindow = [-0.1 1.5]; %@ cue timing
            pcaCondAvoid = 'allUnitsTrialsAvoidAtCueMean';
            pcaCondReact = 'allUnitsTrialsReactAtCueMean';
        else %else @ VHEx timing
            if removeNull
                condAvoid = 'allUnitsTrialsAvoidNullMean';
                condReact = 'allUnitsTrialsReactNullMean';
                condShuf1 = 'allUnitsTrialsShuf1NullMean';
                condShuf2 = 'allUnitsTrialsShuf2NullMean';
                titleStr = [labelName '(-null) at VHEx'];
            else
                condAvoid = 'allUnitsTrialsAvoidMean';
                condReact = 'allUnitsTrialsReactMean';
                condShuf1 = 'allUnitsTrialsShuf1Mean';
                condShuf2 = 'allUnitsTrialsShuf2Mean';
                titleStr = [labelName ' at VHEx'];
            end
            pcaWindow = [-1.5 0.1]; %@VHEx timing
            pcaCondAvoid = 'allUnitsTrialsAvoidMean';
            pcaCondReact = 'allUnitsTrialsReactMean';
        end
    
        numMiceProbes = size(allData.(probeName).(anatName).(condAvoid),2);
    %     numDimensions = size(allDataSubspace.(probeName).(anatName).(condAvoid){1}{1},1); %subspace dimensions for normalization below
        plotWindow = [-2 2];
        plotSamps = [find(plotWindow(1)>=allData.timescaleSeg,1,'last') find(plotWindow(2)>=allData.timescaleSeg,1,'last')];
        pcaSamps = [find(pcaWindow(1)>=allData.timescaleSeg,1,'last') find(pcaWindow(2)>=allData.timescaleSeg,1,'last')];
        timescalePlot = allData.timescaleSeg(plotSamps(1):plotSamps(2));
        %other constants
        fontSz = 16;
        cmapLine = allData.lineCmaps(cmapIdx,:);
        numRepeats = 1;
        usePcContrast = false; % whether we perform contrastive PCA on data1 vs data2, or vanilla PCA on combined data
        smoothingSamples = [10 0]; %ex. smooth window 25 is 25*20ms PSTH bins = 500ms; for causal half-kernel use [X 0]
    
        savedEuclDistanceRpts = [];
        savedEuclDistanceRptsShuf = [];
        % for each mouse/repeat, compute Euclidean distance from subspace trial means, then average over trial sampling repeats
        for i = 1:numMiceProbes
%         for i = 1:5 % local opto option
%         for i = 6:10 % large scale opto option
            thisMouseAvoidData = allData.(probeName).(anatName).(condAvoid){i}; %unpack to get per-mouse cells that include all sample repeats
            thisMouseReactData = allData.(probeName).(anatName).(condReact){i};
            thisMouseShuf1Data = allData.(probeName).(anatName).(condShuf1){i};
            thisMouseShuf2Data = allData.(probeName).(anatName).(condShuf2){i};
            numUnits = size(thisMouseAvoidData{1},1);
            numPCsToUse = numUnits; %use numUnits here for full dimensionality in distance calculation
            
            for rpt = 1:numRepeats
                thisRptAvoidData = thisMouseAvoidData{rpt}; %unpack one more time to get each repeat cell - now this cell for each trial mean
                thisRptReactData = thisMouseReactData{rpt};
                thisRptShuf1Data = thisMouseShuf1Data{rpt};
                thisRptShuf2Data = thisMouseShuf2Data{rpt};

                % for each sampling repeat, get PCs and distances for each mouse
                % always use data without null removal for PC calculation, then only apply null data when projecting, if necessary
                data1 = allData.(probeName).(anatName).(pcaCondAvoid){i}{rpt};
                data2 = allData.(probeName).(anatName).(pcaCondReact){i}{rpt};
                if usePcContrast
                    % re-center data; already z-scored but need to re-center after making trial segment means
                    data1All = data1(:, pcaSamps(1):pcaSamps(2));
                    data1Cent = data1All - mean(data1All,2);
                    data2All = data2(:, pcaSamps(1):pcaSamps(2));
                    data2Cent = data2All - mean(data2All,2);
                    % contrastive PCA: subtract covariances of above matrices, then find PCs of that matrix to get linear combinations that maximize the differences in mean activity patterns between the cases
                    [PC,~,VarExplained] = pcaContrast(data1Cent, data2Cent);  %PCs in columns in decreasing VarExplained
                else
                    dataCombined = [data1 data2];
                    dataCent = dataCombined - mean(dataCombined,2);
%                     [PC,~,VarExplained] = pcaCov(dataCent);
                    % option to calculate on data instead of covariance:
                    [U,S,V] = svd(dataCent');
                    varExplained = 100*(diag(S).^2)./sum(diag(S).^2);
%                     if i==1
%                         fv(p) = figure; plot(cumsum(varExplained)); hold on; title(params(p).labelName);
%                     else
%                         figure(fv(p)); plot(cumsum(varExplained));
%                     end
                    PC = V;
                end

                avoidMeanProj = PC'*thisRptAvoidData(:,plotSamps(1):plotSamps(2)); %apply PCs to the mean responses across trials
                reactMeanProj = PC*thisRptReactData(:,plotSamps(1):plotSamps(2));
                avoidMeanProj = smoothdata(avoidMeanProj, 2, 'gaussian', smoothingSamples); %ex. smooth window 25 is 25*20ms PSTH bins = 500ms
                reactMeanProj = smoothdata(reactMeanProj, 2, 'gaussian', smoothingSamples);

                shuf1MeanProj = PC'*thisRptShuf1Data(:,plotSamps(1):plotSamps(2)); %apply same PCs as above to shuffled means
                shuf2MeanProj = PC'*thisRptShuf2Data(:,plotSamps(1):plotSamps(2));
                shuf1MeanProj = smoothdata(shuf1MeanProj, 2, 'gaussian', smoothingSamples); %ex. smooth window 25 is 25*20ms PSTH bins = 500ms
                shuf2MeanProj = smoothdata(shuf2MeanProj, 2, 'gaussian', smoothingSamples);

                euclDistance = zeros(size(avoidMeanProj,2)-1,1);
                euclDistanceShuf = zeros(size(avoidMeanProj,2)-1,1);
    %             figure; imagesc(thisRptAvoidData); clim([-1 1]);
    %             figure; imagesc(thisRptReactData); clim([-1 1]);
                for dataPoint = 1:size(avoidMeanProj,2)-1
                    euclDistance(dataPoint) = norm(avoidMeanProj(1:numPCsToUse,dataPoint) - reactMeanProj(1:numPCsToUse,dataPoint));
                    euclDistanceShuf(dataPoint) = norm(shuf1MeanProj(1:numPCsToUse,dataPoint) - shuf2MeanProj(1:numPCsToUse,dataPoint));
                end
                % calculate Euclidean distance between trajectories over time (normalized by sqrt(# subspace dimensions) for equal comparison with or without null removed)
                savedEuclDistanceRpts(rpt, i, :) = euclDistance./sqrt(numPCsToUse); %save for bounded line below
                savedEuclDistanceRptsShuf(rpt, i, :) = euclDistanceShuf./sqrt(numPCsToUse); 
            end %end sampling repeats
        end %end numMice
    
        % now take mean across repeats for each mouse:
        savedEuclDistance = squeeze(mean(savedEuclDistanceRpts,1));
        savedEuclDistanceShuf = squeeze(mean(savedEuclDistanceRptsShuf,1));
        
       figure(f) 
       nexttile; hold on;
        % optionally plot individual mice to see variability:
%         for i = 1:numMiceProbes
% %             plot(timescalePlot(1:end-1), savedEuclDistance(i,:), 'Color', cmapLine, 'LineWidth', 0.1)
%             plot(timescalePlot(1:end-1), savedEuclDistanceShuf(i,:), 'Color', [0 0 0], 'LineWidth', 0.1)
%         end
        [meanEuclDistance, semEuclDistance] = grpstats(savedEuclDistance,[],{'mean' 'sem'}); %mean across columns/timepoints
        [h1(p),~] = boundedline(timescalePlot(1:end-1), meanEuclDistance, semEuclDistance, 'cmap', cmapLine, 'alpha','y', 'transparency', 0.5, 'LineWidth', 3);
        [meanEuclDistanceShuf, semEuclDistanceShuf] = grpstats(savedEuclDistanceShuf,[],{'mean' 'sem'}); %mean across columns/timepoints
        [h2,~] = boundedline(timescalePlot(1:end-1), meanEuclDistanceShuf, semEuclDistanceShuf, 'cmap', [0 0 0], 'alpha','y', 'transparency', 0.5, 'LineWidth', 3);
        ylim([0 1.5])
%         xlim([timescalePlot(1)+0.5 timescalePlot(end)]) %adjust for smoothing edge effects
        zeroXDottedLine;
        makepretty;
        hold off;
        if p==1
            title('at cue')
        elseif p==2
            title('at VHEx')
        elseif p==5
            ylabel('distance between PC trajectories', 'FontSize', fontSz)
        elseif p==8
            xlabel('time (sec)', 'FontSize', fontSz)
        end

    end %end params loop

    f.Position(3:4) = [1200 900]; %set figure size appropriate for plots
%     legend([h1(1), h1(3), h1(5), h1(7), h2], {params(1).labelName, params(3).labelName, params(5).labelName, params(7).labelName,'shuffle'}, 'Box', 'off', 'Location','northeast');

end

function plotAvoidVsReactPcDistancesWithOptoFromStruct(allData, removeNull)
% input is subspace that already contains trials projected into PCA space, so just take mean across trials and then calculate Euclidean distance

    fontSz = 20;
    blueCmap = [0 0.4470 0.7410];
    orangeCmap = [0.8500 0.3250 0.0980];
    optoCmap = [0 0 0];
    plotWindow = [-2 5];
    % make structure to iterate brain regions in a loop to plot all at once
    params(1).probeName = 'Probe1'; params(1).anatName = 'CTX'; params(1).labelName = 'PFC'; params(1).atCue = true; params(1).cmapIdx = 2;
    params(2).probeName = 'Probe1'; params(2).anatName = 'CTX'; params(2).labelName = 'PFC'; params(2).atCue = false; params(2).cmapIdx = 2;
    params(3).probeName = 'Probe0'; params(3).anatName = 'CTX'; params(3).labelName = 'M1'; params(3).atCue = true; params(3).cmapIdx = 1;
    params(4).probeName = 'Probe0'; params(4).anatName = 'CTX'; params(4).labelName = 'M1'; params(4).atCue = false; params(4).cmapIdx = 1;
    params(5).probeName = 'Probe0'; params(5).anatName = 'BS'; params(5).labelName = 'TH/HY'; params(5).atCue = true; params(5).cmapIdx = 3;
    params(6).probeName = 'Probe0'; params(6).anatName = 'BS'; params(6).labelName = 'TH/HY'; params(6).atCue = false; params(6).cmapIdx = 3;
    params(7).probeName = 'Probe2'; params(7).anatName = 'HB'; params(7).labelName = 'HB'; params(7).atCue = true; params(7).cmapIdx = 4;
    params(8).probeName = 'Probe2'; params(8).anatName = 'HB'; params(8).labelName = 'HB'; params(8).atCue = false; params(8).cmapIdx = 4;

    f = figure; tiledlayout(4,2,"TileSpacing","compact");
    for p = 1:length(params)
    
        atCueOrVhex = params(p).atCue;
        probeName = params(p).probeName;
        anatName = params(p).anatName;
        labelName = params(p).labelName;
        cmapIdx = params(p).cmapIdx;
    
        if atCueOrVhex %true for cue timing
            if removeNull
                condAvoid = 'allUnitsTrialsAvoidAtCueNullMean'; %#ok<*NASGU> 
                condReact = 'allUnitsTrialsReactAtCueNullMean'; 
                condOpto = 'allUnitsTrialsOptoAtCueNullMean';
                condShuf1 = 'allUnitsTrialsShufAtCue1NullMean';
                titleStr = [labelName '(-null) at cue'];
            else
                condAvoid = 'allUnitsTrialsAvoidAtCueMean'; 
                condReact = 'allUnitsTrialsReactAtCueMean'; 
                condOpto = 'allUnitsTrialsOptoAtCueMean';
                condShuf1 = 'allUnitsTrialsShufAtCue1Mean';
                titleStr = [labelName ' at cue'];
            end
            pcaWindow = [-0.1 1.5]; %@ cue timing
            pcaCondAvoid = 'allUnitsTrialsAvoidAtCueMean';
            pcaCondReact = 'allUnitsTrialsReactAtCueMean';
        else %else @ VHEx timing
            if removeNull
                condAvoid = 'allUnitsTrialsAvoidNullMean';
                condReact = 'allUnitsTrialsReactNullMean';
                condOpto = 'allUnitsTrialsOptoNullMean';
                condShuf1 = 'allUnitsTrialsShuf1NullMean';
                titleStr = [labelName '(-null) at VHEx'];
            else
                condAvoid = 'allUnitsTrialsAvoidMean';
                condReact = 'allUnitsTrialsReactMean';
                condOpto = 'allUnitsTrialsOptoMean';
                condShuf1 = 'allUnitsTrialsShuf1Mean';
                titleStr = [labelName ' at VHEx'];
            end
            pcaWindow = [-1.5 0.1]; %@VHEx timing
            pcaCondAvoid = 'allUnitsTrialsAvoidMean';
            pcaCondReact = 'allUnitsTrialsReactMean';
        end
    
        numMiceProbes = size(allData.(probeName).(anatName).(condAvoid),2);
    %     numDimensions = size(allDataSubspace.(probeName).(anatName).(condAvoid){1}{1},1); %subspace dimensions for normalization below
        plotSamps = [find(plotWindow(1)>=allData.timescaleSeg,1,'last') find(plotWindow(2)>=allData.timescaleSeg,1,'last')];
        timescalePlot = allData.timescaleSeg(plotSamps(1):plotSamps(2));
        usePcContrast = false; % whether we perform contrastive PCA on data1 vs data2, or vanilla PCA on combined data
    
        %other constants
        fontSz = 16;
        cmapLine = allData.lineCmaps(cmapIdx,:);
        numRepeats = 10;
        smoothingSamples = [10 0]; %ex. smooth window 25 is 25*20ms PSTH bins = 500ms; for causal half-kernel use [X 0]

        savedOptoAvoidEuclDistanceRpts = [];
        savedOptoReactEuclDistanceRpts = [];
        savedOptoShuffleEuclDistanceRpts = [];
    
        % for each mouse/repeat, compute Euclidean distance from subspace trial means, then average over trial sampling repeats
%         for i = 1:numMiceProbes
        for i = 1:5 % local opto option
%         for i = 6:10 % large scale opto option
            thisMouseAvoidData = allData.(probeName).(anatName).(condAvoid){i}; %unpack to get per-mouse cells that include all sample repeats
            thisMouseReactData = allData.(probeName).(anatName).(condReact){i};
            thisMouseOptoData = allData.(probeName).(anatName).(condOpto){i};
            thisMouseShuf1Data = allData.(probeName).(anatName).(condShuf1){i};
            numUnits = size(thisMouseAvoidData{1},1);
            numPCsToUse = numUnits; %use numUnits here for full dimensionality in distance calculation

            for rpt = 1:numRepeats
                thisRptAvoidData = thisMouseAvoidData{rpt}; %unpack one more time to get each repeat cell - now this cell for each trial mean
                thisRptReactData = thisMouseReactData{rpt};
                thisRptOptoData = thisMouseOptoData{rpt};
                thisRptShuf1Data = thisMouseShuf1Data{rpt};

                % for each sampling repeat, get PCs and distances for each mouse
                % always use data without null removal for PC calculation, then only apply null data when projecting, if necessary
                data1 = allData.(probeName).(anatName).(pcaCondAvoid){i}{rpt};
                data2 = allData.(probeName).(anatName).(pcaCondReact){i}{rpt};
                if usePcContrast
                    % re-center data; already z-scored but need to re-center after making trial segment means
                    data1All = data1(:, pcaSamps(1):pcaSamps(2));
                    data1Cent = data1All - mean(data1All,2);
                    data2All = data2(:, pcaSamps(1):pcaSamps(2));
                    data2Cent = data2All - mean(data2All,2);
                    % contrastive PCA: subtract covariances of above matrices, then find PCs of that matrix to get linear combinations that maximize the differences in mean activity patterns between the cases
                    [PC,~,VarExplained] = pcaContrast(data1Cent, data2Cent);  %PCs in columns in decreasing VarExplained
                else
                    dataCombined = [data1 data2];
                    dataCent = dataCombined - mean(dataCombined,2);
%                     [PC,~,VarExplained] = pcaCov(dataCent);
                    % option to calculate on data instead of covariance:
                    [U,S,V] = svd(dataCent');
                    varExplained = 100*(diag(S).^2)./sum(diag(S).^2);
%                     if i==1
%                         fv(p) = figure; plot(cumsum(varExplained)); hold on; title(params(p).labelName);
%                     else
%                         figure(fv(p)); plot(cumsum(varExplained));
%                     end
                    PC = V;
                end

                avoidMeanProj = PC'*thisRptAvoidData(:,plotSamps(1):plotSamps(2)); %apply PCs to the mean responses across trials
                reactMeanProj = PC*thisRptReactData(:,plotSamps(1):plotSamps(2));
                optoMeanProj = PC'*thisRptOptoData(:,plotSamps(1):plotSamps(2));
                shuf1MeanProj = PC*thisRptShuf1Data(:,plotSamps(1):plotSamps(2));
                avoidMeanProj = smoothdata(avoidMeanProj, 2, 'gaussian', smoothingSamples);
                reactMeanProj = smoothdata(reactMeanProj, 2, 'gaussian', smoothingSamples);
                optoMeanProj = smoothdata(optoMeanProj, 2, 'gaussian', smoothingSamples); 
                shuf1MeanProj = smoothdata(shuf1MeanProj, 2, 'gaussian', smoothingSamples);

                euclDistanceOptoAvoid = zeros(size(avoidMeanProj,2)-1,1);
                euclDistanceOptoReact = zeros(size(avoidMeanProj,2)-1,1);
                euclDistanceOptoShuffle = zeros(size(avoidMeanProj,2)-1,1);

                for dataPoint = 1:size(avoidMeanProj,2)-1
                    euclDistanceOptoAvoid(dataPoint) = norm(optoMeanProj(1:numPCsToUse,dataPoint) - avoidMeanProj(1:numPCsToUse,dataPoint));
                    euclDistanceOptoReact(dataPoint) = norm(optoMeanProj(1:numPCsToUse,dataPoint) - reactMeanProj(1:numPCsToUse,dataPoint));
                    euclDistanceOptoShuffle(dataPoint) = norm(optoMeanProj(1:numPCsToUse,dataPoint) - shuf1MeanProj(1:numPCsToUse,dataPoint));
                end

                % calculate Euclidean distance between trajectories over time (normalized by sqrt(# subspace dimensions) for equal comparison with or without null removed)
                savedOptoAvoidEuclDistanceRpts(rpt, i, :) = euclDistanceOptoAvoid./sqrt(numPCsToUse); %save for bounded line below
                savedOptoReactEuclDistanceRpts(rpt, i, :) = euclDistanceOptoReact./sqrt(numPCsToUse);
                savedOptoShuffleEuclDistanceRpts(rpt, i, :) = euclDistanceOptoShuffle./sqrt(numPCsToUse);
            end %end sampling repeats
        end %end numMice

        % now take mean across repeats for each mouse:
        savedOptoAvoidEuclDistance = squeeze(mean(savedOptoAvoidEuclDistanceRpts,1));
        savedOptoReactEuclDistance = squeeze(mean(savedOptoReactEuclDistanceRpts,1));
        savedOptoShuffleEuclDistance = squeeze(mean(savedOptoShuffleEuclDistanceRpts,1));
    
       nexttile; hold on;
    %     % optionally plot individual mice to see variability:
    %     for i = 1:numMiceProbes
    %         plot(timescalePlot(1:end-1), savedOptoAvoidEuclDistance(i,:), 'Color', orangeCmap, 'LineWidth', 0.1)
    %         plot(timescalePlot(1:end-1), savedOptoReactEuclDistance(i,:), 'Color', blueCmap, 'LineWidth', 0.1)
    % %         plot(timescalePlot(1:end-1), savedOptoShuffleEuclDistance(i,:), 'Color', [0 0 0], 'LineWidth', 0.1)
    %     end
        [meanOptoAvoidEuclDistance, semOptoAvoidEuclDistance] = grpstats(savedOptoAvoidEuclDistance,[],{'mean' 'sem'}); %mean across columns/timepoints
        [h1(p),~] = boundedline(timescalePlot(1:end-1), meanOptoAvoidEuclDistance, semOptoAvoidEuclDistance, 'cmap', orangeCmap, 'alpha','y', 'transparency', 0.5, 'LineWidth', 3);
        [meanOptoReactEuclDistance, semOptoReactEuclDistance] = grpstats(savedOptoReactEuclDistance,[],{'mean' 'sem'}); %mean across columns/timepoints
        [h2(p),~] = boundedline(timescalePlot(1:end-1), meanOptoReactEuclDistance, semOptoReactEuclDistance, 'cmap', blueCmap, 'alpha','y', 'transparency', 0.5, 'LineWidth', 3);
    %     [meanOptoShuffleEuclDistance, semOptoShuffleEuclDistance] = grpstats(savedOptoShuffleEuclDistance,[],{'mean' 'sem'}); %mean across columns/timepoints
    %     [h3(p),~] = boundedline(timescalePlot(1:end-1), meanOptoShuffleEuclDistance, semOptoShuffleEuclDistance, 'cmap', [0 0 0], 'alpha','y', 'transparency', 0.5);

        ylim([0 0.7])
        xlim([timescalePlot(1)+0.5 timescalePlot(end)]) %adjust for smoothing edge effects
        zeroXDottedLine;
        title(titleStr);
        makepretty;
        hold off;
        if p==5
            ylabel('distance between PC trajectories', 'FontSize', fontSz)
        elseif p==8
            xlabel('time (sec)', 'FontSize', fontSz)
        end
    end %end params loop

    f.Position(3:4) = [1200 900]; %set figure size appropriate for plots
%     legend([h1(1) h2(1)], {'opto - avoid', 'opto - react'}, 'Box', 'off');
end

function plotShufMeanPcDistances(allShufPcDist, plotWindow, timescaleSeg)
% take means across all mice & regions & plot
    plotSamps = [find(plotWindow(1)>=timescaleSeg,1,'last') find(plotWindow(2)>=timescaleSeg,1,'last')];
    timescalePlot = timescaleSeg(plotSamps(1):plotSamps(2));
    fontSz = 16;

    figure; hold on;
    [meanEuclDistance, semEuclDistance] = grpstats(allShufPcDist,[],{'mean' 'sem'}); %mean across columns/timepoints
    [hl1,~] = boundedline(timescalePlot(1:end-1), meanEuclDistance, semEuclDistance, 'cmap', [0 0 0], 'alpha','y', 'transparency', 0.5, 'LineWidth', 3);
    ylim([0 0.7])
    xlim([timescalePlot(1)+0.5 timescalePlot(end)]) %adjust for smoothing edge effects
    zeroXDottedLine;
    hold off;
    legend([hl1], {'shuffle mean'}, 'Box', 'off');
    xlabel('time (sec)', 'FontSize', fontSz)
    ylabel('distance between PC trajectories (z-Hz/neuron)', 'FontSize', fontSz)
    makepretty;

end

function drawAnatColormapLegend(anatCmaps, anatNames)
% add anatomical colormap legend to current plot, regardless of what is plotted, using hidden scatterplots
    h=gobjects(length(anatNames),1);
    for i = 1:length(anatNames)
        h(i)=plot(NaN,NaN, 'Color', anatCmaps(i,:), 'MarkerFaceColor', anatCmaps(i,:), 'Marker', 'o', 'MarkerSize', 8, 'MarkerEdgeColor', 'k', 'LineStyle', 'none'); % control the marker type in the legend here 
    end
    legend(h, anatNames)
end

function [subsetIdx, notSubsetIdx] = getAnatSubset(anatHiers, anatName)
% for unit data with corresponding anatHiers, go through each unit to find indeces corresponding to anatomical specifiers
    subsetIdx = [];
    notSubsetIdx = [];
    for unitIdx = 1:length(anatHiers) % loop over all units to find matches to anatName
        currentAnatHiers = cell2mat(anatHiers{unitIdx});
        if contains(currentAnatHiers,anatName) 
            subsetIdx(end+1) = unitIdx;
        else
            notSubsetIdx(end+1) = unitIdx;
        end
    end
end

function plotUnitAnatomyPerProbe(ProbeX, titleStr)
% plot # units per anatomical area in list below on each probe
    numUnitsThisProbe = size(ProbeX.allUnitsMeanZFR,1);

    anatNames = {'HPF','TH', 'HY','CTX', 'HB', 'MB',  'CNU', 'BS'}; %based on Allen Reference taxonomy (see 'gui_data.st' in atlas script, https://atlas.brain-map.org/)
    anatCmaps = linspecer(length(anatNames)); %corresponding colormaps to index for plotting

    anatCounts = zeros(length(anatNames),1); 
    
    for unitIdx = 1:numUnitsThisProbe % loop over all units to find matches to anatName
        currentAnatHiers = cell2mat(ProbeX.anatHiers{unitIdx});
        
        currentAnatIdx = 1;
        for anatIdx = anatNames
            if contains(currentAnatHiers,anatIdx{1}) 
                anatCounts(currentAnatIdx) = anatCounts(currentAnatIdx) + 1;
                break; % go to next unit
            else
                currentAnatIdx = currentAnatIdx + 1;
                continue; % keep going down anatomical list, if nothing found current cmap stays default
            end
        end
    end

    figure; hold on;
    for barIdx = 1:length(anatNames)
        bar(barIdx, anatCounts(barIdx), 'FaceColor', anatCmaps(barIdx,:))
    end
    axis tight;
    xticklabels(anatNames)
    title(titleStr)
    drawAnatColormapLegend(anatCmaps, anatNames) %test if this matched labels
    hold off;
    set(gcf,'color','w'); %set figure background white
    makepretty
end

function [unitCmaps, unitCmapsIdx, anatCmaps, anatNames] = anatHiersToCmaps(allUnitsAnatHiers)
% take array of cells containing anatomical hierarchies for each unit and return an array of colormaps
    % anatNames here is ordered by hierarchical depth, deepest entries first, as inner loop will assign those colors first; others at the same hierarchical level are mutually exclusive
    % any units not in the list are assigned default map
%     anatNames = {'HPF','TH', 'HY','CTX', 'HB', 'MB',  'CNU', 'BS'}; %based on Allen Reference taxonomy (see 'gui_data.st' in atlas script, https://atlas.brain-map.org/)
    anatNames = {'ORB', 'DP', 'TT', 'MOp','CTX', 'BS'}; %alternative map to highlight cortex subregions
    anatCmaps = linspecer(length(anatNames)); %corresponding colormaps to index for plotting
    % abbr.     name                    atlas taxonomy hierarchy depth
    % -----------------------------------------------------------------
    % BS        brainstem               2
    % CTX       cerebral cortex         3   % also make "nonCTX" structure with compliment, as a special case
    % CNU       cerebral nuclei         3
    % HB        hindbrain               3
    % MB        midbrain                3
    % TH        thalamus                4
    % HY        hypothalamus            4
    % HPF       hippocampal formation   5
    
    defaultCmap = [0 0 0]; %default black map for neurons not matching any in anatomical list above
    unitCmaps = repmat(defaultCmap, length(allUnitsAnatHiers),1);
    unitCmapsIdx = zeros(length(allUnitsAnatHiers),1); % also return array of indeces into color list to change position, etc. of color easily

    for unitIdx = 1:length(allUnitsAnatHiers) % loop over all units to find matches to anatName
        currentAnatHiers = cell2mat(allUnitsAnatHiers{unitIdx});
        
        currentAnatIdx = 1;
        for anatIdx = anatNames
            if contains(currentAnatHiers,anatIdx{1}) 
                unitCmaps(unitIdx,:) = anatCmaps(currentAnatIdx,:);
                unitCmapsIdx(unitIdx) = currentAnatIdx;
                break; % go to next unit
            else
                currentAnatIdx = currentAnatIdx + 1;
                continue; % keep going down anatomical list, if nothing found current cmap stays default
            end
        end
    end

end

%% misc. PCA functions
function [PC, varExplained] = getVanillaPCsEphys(allUnitsMeanZFR, timescaleSeg, pcaWindow) 
%get vanilla PCs across all VHEx or cues
    pcaSamps = [find(pcaWindow(1)>=timescaleSeg,1,'last') find(pcaWindow(2)>=timescaleSeg,1,'last')];


    % first trim data to window & mean subtract, even though mean z-score data should be close to centered already
    allData = allUnitsMeanZFR(:, pcaSamps(1):pcaSamps(2));
    allDataCent = allData - mean(allData,2);

    % contrastive PCA: subtract covariances of above matrices, then find PCs of that matrix to get linear combinations that maximize the differences in mean activity patterns between the cases
    [PC,V,varExplained] = pcaCov(allDataCent); % option to just use standard PCs calcluated across all VHEx

end

function [pcSeg, varExplained] = getEphysPCs(dataCombined, data1, data2, numGoodUnits, numMice, numPcsToUse, pcaSamps, useContrast, alignPCs)
% function to get PCs organized by mouse session (i.e. each has a PC segment)
% 'contrast' input controls whether we perform contrastive PCA on data1 vs data2, or vanilla PCA on dataCombined 
% 'align' input controls whether PCs are aligned in common subspace or computed separately for each mouse/session
% numPcsToUse is empirically determined number of PCs to use with shared subspace (using varExplained below), which should not be more than than the smallest # neurons in a session/region
    
    if alignPCs 
        % re-center data; already z-scored but need to re-center after making trial segment means
        data1All = data1(:, pcaSamps(1):pcaSamps(2));
        data1Cent = data1All - mean(data1All,2);
        data2All = data2(:, pcaSamps(1):pcaSamps(2));
        data2Cent = data2All - mean(data2All,2);
        dataCombinedAll = dataCombined(:, pcaSamps(1):pcaSamps(2));
        allDataCent = dataCombinedAll - mean(dataCombinedAll,2);
        % contrastive PCA: subtract covariances of above matrices, then find PCs of that matrix to get linear combinations that maximize the differences in mean activity patterns between the cases
        if useContrast
            [PC,V,varExplained] = pcaContrast(data1Cent, data2Cent);  %PCs in columns in decreasing VarExplained
        else
            [PC,V,varExplained] = pcaCov(allDataCent);  %optionally try normal PCA 
        end
        PCtoUse = PC(:,1:numPcsToUse); %use only top PCs 
        
        %now re-orthogonalize for each mouse/session:
        currentEndUnit = 0;
        for i = 1:numMice
            currentStartUnit = currentEndUnit + 1;
            currentEndUnit = currentEndUnit + numGoodUnits(i);
            pcSegOrig{i} = PCtoUse(currentStartUnit:currentEndUnit,:); %take the subset corresponding to this mouse/session
            [Q,R] = mgs(pcSegOrig{i}); % modified Gram-Schmidt orthogonalization option; prevents rows and columns in Q and R from fliping their signs
            % [Q] = pcSegOrig{i}; %option for no re-orthogonalization
            pcSeg{i} = Q(:,1:numPcsToUse); % first numPcsToUse columns of Q are orthonormal basis of pcSeg{i}; can also use 'econ' option with qr() above and just compute first numPcsToUse
    %         diffQ = pcSegOrig{i} - pcSeg{i}; %optionally look at diff between original and re-orthogonalized
    %         figure; tiledlayout(1,3); 
    %         nexttile; imagesc(pcSegOrig{i}); nexttile; imagesc(pcSeg{i}); nexttile; imagesc(diffQ); %optionally view difference between original PC matric and orthogonalized version
        end
    else
        currentEndUnit = 0;
        for i = 1:numMice
            currentStartUnit = currentEndUnit + 1;
            currentEndUnit = currentEndUnit + numGoodUnits(i);

            % re-center data; already z-scored but need to re-center after making trial segment means
            data1All = data1(currentStartUnit:currentEndUnit, pcaSamps(1):pcaSamps(2));
            data1Cent = data1All - mean(data1All,2);
            data2All = data2(currentStartUnit:currentEndUnit, pcaSamps(1):pcaSamps(2));
            data2Cent = data2All - mean(data2All,2);
            dataCombinedAll = dataCombined(currentStartUnit:currentEndUnit, pcaSamps(1):pcaSamps(2));
            allDataCent = dataCombinedAll - mean(dataCombinedAll,2);
            % contrastive PCA: subtract covariances of above matrices, then find PCs of that matrix to get linear combinations that maximize the differences in mean activity patterns between the cases
            if useContrast
                [pcSeg{i},V{i},varExplained{i}] = pcaContrast(data1Cent, data2Cent);  %PCs in columns in decreasing VarExplained
            else
                [pcSeg{i},V{i},varExplained{i}] = pcaCov(allDataCent);  %optionally try normal PCA
            end
        end
    end

end

function [pcSeg, varExplained] = getEphysPCsFromData(dataCombined, numGoodUnits, numMice, numPcsToUse, alignPCs)
% function to get PCs organized by mouse/session (i.e. each has a PC segment)
% numPcsToUse is empirically determined number of PCs to use with shared subspace (using varExplained below), which should not be more than than the smallest # neurons in a session/region
    
% this version calculates PCs directly from trial means, rather than covariance

    if alignPCs
        [U,S,V] = svd(dataCombined'); %use svd directly to prevent sign flip in pca()...for comparison: [coeff,score,latent,tsquared,explained,mu] = pca(allDataCent');
%         covDataCombined = cov(dataCombined'); % optionally see effect if calculating from covariance instead
%         [U,S,V] = svd(covDataCombined);
        varExplained = (diag(S).^2)./sum(diag(S).^2);
        PCtoUse = V(:,1:numPcsToUse); %use only top PCs 
        
        %now re-orthogonalize for each mouse/session:
        currentEndUnit = 0;
        for i = 1:numMice
            currentStartUnit = currentEndUnit + 1;
            currentEndUnit = currentEndUnit + numGoodUnits(i);
            pcSegOrig{i} = PCtoUse(currentStartUnit:currentEndUnit,:); %take the subset corresponding to this mouse/session
            [Q,R] = mgs(pcSegOrig{i}); % modified Gram-Schmidt orthogonalization option; prevents rows and columns in Q and R from fliping their signs
            pcSeg{i} = Q(:,1:numPcsToUse); % first numPcsToUse columns of Q are orthonormal basis of pcSeg{i}
        end
    else
        currentEndUnit = 0;
        for i = 1:numMice
            currentStartUnit = currentEndUnit + 1;
            currentEndUnit = currentEndUnit + numGoodUnits(i);

            % re-center data; already z-scored but need to re-center after making trial segment means
            dataCombinedThisMouse = dataCombined(currentStartUnit:currentEndUnit, :);
            [U,S,V] = svd(dataCombinedThisMouse');
            pcSeg{i} = V;
            varExplained{i} = (diag(S).^2)./sum(diag(S).^2);
        end
    end
end

function [pcSeg, varExplained] = getEphysPCsWithNullspace(dataCombined, data1, data2, dataForNull, numGoodUnits, numMice, numPcsToUse, pcaSamps, nullSamps, useContrast)
% function to get PCs organized by mouse session (i.e. each has a PC segment)
% assume alignment into common subspace, but useContrast decides whether to do cPCA on data1-data2 or vPCA on dataCombined 

    % re-center data; already z-scored but need to re-center after making trial segment means
    dataCombinedAll = dataCombined(:, pcaSamps(1):pcaSamps(2));
    allDataCent = dataCombinedAll - mean(dataCombinedAll,2);
    covForAllData = cov(allDataCent');
    covForAllDataNullRemovedSpd = nearestSPD(covForAllData);

    data1All = data1(:, pcaSamps(1):pcaSamps(2));
    data1Cent = data1All - mean(data1All,2);
    covForData1 = cov(data1Cent');
    data2All = data2(:, pcaSamps(1):pcaSamps(2));
    data2Cent = data2All - mean(data2All,2);
    covForData2 = cov(data2Cent');
    covSubtracted = covForData1 - covForData2;
    covSubtractedSpd = nearestSPD(covSubtracted);

    dataNullAll = dataForNull(:, nullSamps(1):nullSamps(2));
    nullDataCent = dataNullAll - mean(dataNullAll,2);
    covForNull = cov(nullDataCent'); % use common subspace for both nullspace and PC calculation
    [PCnull,Vnull,varExplainedNull] = pcacov(covForNull);
%     covForNull = covSubtractedSpd; % sanity check option to use same cov as PC distance data, which should get rid of most PC distance
    nullForData = null(covForNull); %,0.0001); %use reasonable tolerance of 0.0001 (otherwise can end up with empty null matrix for cov matrix that is barely full rank); this will give a nullspace [# neurons x # null dimensions], where # null dimensions is less than the number of cells

    % ex.
%     [U, S, V] = svd(covForNull);
%     tol = max(size(covForNull)) * eps(norm(S, 'fro')); %0.0001; %
%     r = sum(diag(S) > tol); %rank
% %     nullForData = V(:, r+1:end);
%     nullForData = V(:, r+1:r+1+numPcsToUse);
% %     nullForData = V(:, end-100:end);

    if useContrast
        % can also project trial-averaged data stself into nullspace here, but working is shared subspace is simpler
        % calculate shared PC subspace based on original data, and then project that subspace onto the nullspace below
        [PC,V,varExplained] = pcacov(covSubtractedSpd);        
    else
        [PC,V,varExplained] = pcacov(covForAllDataNullRemovedSpd);
    end
%     PCNullRemoved = nullProj*PC; %project PCs onto both rows & columns of nullspace
%     PCtoUse = PCNullRemoved(:,1:numPcsToUse);
    PCtoUse = PC(:,1:numPcsToUse);

    %now re-orthogonalize for each mouse/session:
    currentEndUnit = 0;
    nullDimToUse = numPcsToUse; % # null components s/b <= numPcsToUse, so that we're not artificially re-injecting dimensionality (i.e. covariance patterns) that will counteract the intended effect of the null projection; using less makes the null projection stronger
    for i = 1:numMice
        currentStartUnit = currentEndUnit + 1;
        currentEndUnit = currentEndUnit + numGoodUnits(i);
        pcSegOrig{i} = PCtoUse(currentStartUnit:currentEndUnit,:); %take the subset corresponding to this mouse/session
        nullSeg = nullForData(currentStartUnit:currentEndUnit,1:nullDimToUse); %take matching nullspace subset; note that including more than numPcsToUse dimensions here will effectively create a projection matrix that expands the PC matrix into 

        [Qn,~] = mgs(nullSeg); %also re-orthogonalize null
        nullProj = Qn * Qn'; %create null prjection matrix

        [Q,R] = mgs(pcSegOrig{i}); % modified Gram-Schmidt orthogonalization option; prevents rows and columns in Q and R from fliping their signs
%         pcSeg{i} = pcSegOrig{i}; %option for no re-othogonalization
        pcSeg{i} = nullProj*Q; %project PCs onto nullspace AFTER re-othogonalization, so that QR decomposition doesn't interfere with null projection

%         diffQ = pcSegOrig{i} - pcSeg{i}; %optionally look at diff between original and re-orthogonalized
%         figure; tiledlayout(1,3); 
%         nexttile; imagesc(pcSegOrig{i}); nexttile; imagesc(pcSeg{i}); nexttile; imagesc(diffQ); %optionally view difference between original PC matric and orthogonalized version
    end

end

function [pcSeg, varExplained] = getKinPCs(dataCombined, data1, data2, numKpts, numMice, pcaSamps, useContrast, alignPCs) 
% function to get PCs organized by mouse session (i.e. each has a PC segment)
% 'contrast' input controls whether we perform contrastive PCA on data1 vs data2, or vanilla PCA on dataCombined 
% 'align' input controls whether PCs are aligned in common subspace or computed separately for each mouse/session

    if alignPCs
        numPcsToUse = 10; % empirically determined number of PCs to use with shared subspace (using varExplained below); since there are always same # keypoints just use that many to match the independent PC case
        % kin data is centered already
        data1Cent = data1(:, pcaSamps(1):pcaSamps(2));
        data2Cent = data2(:, pcaSamps(1):pcaSamps(2));
        allDataCent = dataCombined(:, pcaSamps(1):pcaSamps(2));
        % contrastive PCA: subtract covariances of above matrices, then find PCs of that matrix to get linear combinations that maximize the differences in mean activity patterns between the cases
        if useContrast
            [PC,V,varExplained] = pcaContrast(data1Cent, data2Cent);  %PCs in columns in decreasing VarExplained
        else
            [PC,V,varExplained] = pcaCov(allDataCent);  %optionally try normal PCA 
        end
        PCtoUse = PC(:,1:numPcsToUse); %use only top PCs 
        
        %now re-orthogonalize for each mouse/session:
        currentEndUnit = 0;
        for i = 1:numMice
            currentStartUnit = currentEndUnit + 1;
            currentEndUnit = currentEndUnit + numKpts(i);
            pcSegOrig{i} = PCtoUse(currentStartUnit:currentEndUnit,:); %take the subset corresponding to this mouse/session
            [Q,R] = mgs(pcSegOrig{i}); % modified Gram-Schmidt orthogonalization option; prevents rows and columns in Q and R from fliping their signs
            pcSeg{i} = Q(:,1:numPcsToUse); % first numPcsToUse columns of Q are orthonormal basis of pcSeg{i}; can also use 'econ' option with qr() above and just compute first numPcsToUse
    %         diffQ = pcSegOrig{i} - pcSeg{i}; %optionally look at diff between original and re-orthogonalized
    %         figure; tiledlayout(1,3); 
    %         nexttile; imagesc(pcSegOrig{i}); nexttile; imagesc(pcSeg{i}); nexttile; imagesc(diffQ); %optionally view difference between original PC matric and orthogonalized version
        end
    else
        currentEndUnit = 0;
        for i = 1:numMice
            currentStartUnit = currentEndUnit + 1;
            currentEndUnit = currentEndUnit + numKpts(i);

            % kin data is centered already
            data1Cent = data1(currentStartUnit:currentEndUnit, pcaSamps(1):pcaSamps(2));
            data2Cent = data2(currentStartUnit:currentEndUnit, pcaSamps(1):pcaSamps(2));
            allDataCent = dataCombined(currentStartUnit:currentEndUnit, pcaSamps(1):pcaSamps(2));
            % contrastive PCA: subtract covariances of above matrices, then find PCs of that matrix to get linear combinations that maximize the differences in mean activity patterns between the cases
            if useContrast
                [pcSeg{i},V{i},varExplained{i}] = pcaContrast(data1Cent, data2Cent);  %PCs in columns in decreasing VarExplained
            else
                [pcSeg{i},V{i},varExplained{i}] = pcaCov(allDataCent);  %optionally try normal PCA
            end
        end
    end

end

function [PC, V, varExplained] = pcaCov(data)
    % pcaCov: Perform PCA using SVD
    %
    %  data         - NxT matrix of input data (already mean subtracted & normalized/z-scored)
    %                 (N neurons/dimensions, T trials/timepoints)
    %  PC           - each column is a principal component / singular vector
    %  V            - Nx1 matrix of singular values in decreasing order
    %  varExplained - V ./ sum(V) *100;
    
    %[N, T] = size(data);

    % calculate the NxN covariance matrix
    % since the covariance matrix is symmetric and thus transforms a
    % vector by stretching or shrinking it along its eigenvectors, proportional to the eigenvalues;
    % so, computing the eigen/singular vectors are the linear combinations
    % of neurons that contribute to the directions of greatest variance in
    % the firing rate covariance pattern 
    covariance = cov(data');
%     figure; imagesc(covariance);
%     title('standard PCA cov matrix')
%     colorbar; %clim([0 0.05]);
    % find the N singular vectors and values; use pcacov() (uses SVD) instead of eig() to take care of possibility of small complex eigenvalues
    [PC,V,varExplained] = pcacov(covariance);
end

function [PC, V, varExplained] = pcaContrast(data1, data2)
    % pcaContrast: Perform contrastive PCA using SVD of subtracted covariance matrices (see Abid et al. 2018)
    % NOTE: this will create a DC offset in the PC loadings because the
    % subtracted covariance matrix will no longer be centered at ~0
    %
    %  data1,2      - NxT matrices of input data (already mean subtracted & normalized/z-scored)
    %                 (N neurons/dimensions, T trials/timepoints)
    %  PC           - each column is a principal component / singular vector
    %  V            - Nx1 matrix of singular values in decreasing order
    %  varExplained - % of variance explained by each vector/value
    
    %[N1, T1] = size(data1);
    %[N2, T2] = size(data2);

    %"The contrast parameter α represents the trade-off between having the 
    % high target variance and the low background variance. When α = 0, cPCA 
    % selects the directions that only maximize the target variance, and hence 
    % reduces to PCA applied on the target data {x i }. As α increases, 
    % directions with smaller background variance become more important and the 
    % cPCs are driven towards the null space of the background data {y i }. In 
    % the limiting case α = ∞, any direction not in the null space of {y i } 
    % receives an infinite penalty. In this case, cPCA corresponds to first 
    % projecting the target data onto the null space of the background data, 
    % and then performing PCA on the projected data"
    alpha = 1; % keep at 1 for simplicity, interperability, and symmetry between data1 & data2

    % calculate the NxN covariance matrix
    covariance1 = cov(data1'); 
    covariance2 = cov(data2'); 
    covSubtracted = covariance1 - covariance2.*alpha;
    covSubtractedSpd = nearestSPD(covSubtracted); %get nearest positive semi-definite version for pcacov
%     figure; imagesc(covSubtractedSpd);
%     colorbar; clim([0 0.05]);
%     title('contrastive PCA cov matrix')
    % find the N singular vectors and values; use pcacov() (uses SVD) instead of eig() to take care of possibility of small complex eigenvalues
    [PC,V,varExplained] = pcacov(covSubtractedSpd);

%     [V,varExplainedStruct,PC] = gcPCA(data1', data2', 4.1); %option to use generalized contrastive PCA algorithm to normalize by added cov matrices (see Sjulson et al. 2024)
%     varExplained = varExplainedStruct.a;
end

function [Q,R] = mgs(X)
% see https://blogs.mathworks.com/cleve/2016/07/25/compare-gram-schmidt-and-householder-orthogonalization-algorithms/?s_tid=answers_rc2-1_p4_BOTH#56bc9f27-5693-4279-a435-f45848422d17
    % Modified Gram-Schmidt.  [Q,R] = mgs(X); more numerically-stable (better orthogonality) than classic GS
    % G. W. Stewart, "Matrix Algorithms, Volume 1", SIAM, 1998.
    [n,p] = size(X);
    Q = zeros(n,p);
    R = zeros(p,p);
    for k = 1:p
        Q(:,k) = X(:,k);
        for i = 1:k-1
            R(i,k) = Q(:,i)'*Q(:,k);
            Q(:,k) = Q(:,k) - R(i,k)*Q(:,i);
        end
        R(k,k) = norm(Q(:,k))';
        Q(:,k) = Q(:,k)/R(k,k);
    end
end







