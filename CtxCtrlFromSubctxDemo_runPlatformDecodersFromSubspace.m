% Jason Keller
% 2024
% load saved decoder data in structure and run decoders, save output in MAT file

rootFolder = 'C:\Users\labadmin\Desktop\';
decodeAndSave = 1;
kinematicsData = 1;
ephysData = 0;
useNetOrLSTM = 0; % (1) neuralNet/LSTM or (0) SVRegression; neuralNet/LSTM not currently used
numDecoderSampRepeats = 10;

medGreyCmap = [0.4 0.4 0.4];
lightGreyCmap = [0.7 0.7 0.7];
blueCmap = [0 0.4470 0.7410];
orangeCmap = [0.8500 0.3250 0.0980];
fontSz = 16;

if decodeAndSave
    decoderMatNameEphys = 'dataToDecodePlatform.mat';
    decoderMatNameKinematics = 'kinDataToDecodePlatform.mat';
    
    if ephysData
        load([rootFolder decoderMatNameEphys]);
    end
    if kinematicsData
        load([rootFolder decoderMatNameKinematics]); %#ok<*UNRCH> 
    end
    
    if ephysData
        probeToPlot = {'Probe0','Probe1','Probe0','Probe2'}; %matched pairs with below
        anatToPlot =  {'CTX',   'CTX',   'BS',    'HB'};
        labels =      {'M1',    'PFC',   'TH/HY', 'HB'}; 

        for rpt = 1:numDecoderSampRepeats
            disp(['Processing  ephys sampling repeat ' num2str(rpt) ' of ' num2str(numDecoderSampRepeats)])
            % first run ephys data:
            for probeAnatIdx = 1:length(probeToPlot)
                disp(['Processing ' probeToPlot{probeAnatIdx} ' ' anatToPlot{probeAnatIdx}])
                if useNetOrLSTM
                    %not used                              
                else
                    d = dataToDecodePlatform.(probeToPlot{probeAnatIdx}).(anatToPlot{probeAnatIdx});
%                     mixedPredictors = cat(2, d.allUnitsTrialsReact, d.allUnitsTrialsAvoid, d.allUnitsTrialsOpto);
%                     mixedResponses = cat(2, d.platDistZReact, d.platDistZAvoid, d.platDistZOpto);
%                     concatReactP = cat(2, d.allUnitsTrialsReact, d.allUnitsTrialsReact, d.allUnitsTrialsReact);%to match partitioning indeces with above
%                     concatReactR = cat(2, d.platDistZReact, d.platDistZReact, d.platDistZReact);
%                     concatAvoidP = cat(2, d.allUnitsTrialsAvoid, d.allUnitsTrialsAvoid, d.allUnitsTrialsAvoid);%to match partitioning indeces with above
%                     concatAvoidR = cat(2, d.platDistZAvoid, d.platDistZAvoid, d.platDistZAvoid);
%                     concatOptoP = cat(2, d.allUnitsTrialsOpto, d.allUnitsTrialsOpto, d.allUnitsTrialsOpto);%to match partitioning indeces with above
%                     concatOptoR = cat(2, d.platDistZOpto, d.platDistZOpto, d.platDistZOpto);
%                     disp('train on all')
%                     [dd.rSquaredMixedMixed, dd.rSquaredMixedAvoid, ~, ~, ~, ~] = getPlatformDecoderSVR(mixedResponses(rpt,:)', squeeze(mixedPredictors(rpt,:,:)), concatAvoidR(rpt,:)', squeeze(concatAvoidP(rpt,:,:)), concatOptoR(rpt,:)', squeeze(concatOptoP(rpt,:,:)));
%                     [~, dd.rSquaredMixedReact, ~, ~, ~, ~] = getPlatformDecoderSVR(mixedResponses(rpt,:)', squeeze(mixedPredictors(rpt,:,:)), concatReactR(rpt,:)', squeeze(concatReactP(rpt,:,:)), concatOptoR(rpt,:)', squeeze(concatOptoP(rpt,:,:)));
%                     [~, dd.rSquaredMixedOpto, ~, ~, ~, ~] = getPlatformDecoderSVR(mixedResponses(rpt,:)', squeeze(mixedPredictors(rpt,:,:)), concatOptoR(rpt,:)', squeeze(concatOptoP(rpt,:,:)), concatOptoR(rpt,:)', squeeze(concatOptoP(rpt,:,:)));
%                     decodedDataAll.(probeToPlot{probeAnatIdx}).(anatToPlot{probeAnatIdx})(rpt) = dd;
%                     clear dd;

                    disp('train on react')
                    [dd.rSquaredReact, dd.rSquaredReactAvoid, dd.rSquaredReactOpto, dd.lagReact, dd.lagReactAvoid, dd.lagReactOpto] = getPlatformDecoderSVR(d.platDistZReact(rpt,:)', squeeze(d.allUnitsTrialsReact(rpt,:,:)), d.platDistZAvoid(rpt,:)', squeeze(d.allUnitsTrialsAvoid(rpt,:,:)), d.platDistZOpto(rpt,:)', squeeze(d.allUnitsTrialsOpto(rpt,:,:)));
                    disp('train on avoid')
                    [dd.rSquaredAvoid, dd.rSquaredAvoidReact, dd.rSquaredAvoidOpto, dd.lagAvoid, dd.lagAvoidReact, dd.lagAvoidOpto] = getPlatformDecoderSVR(d.platDistZAvoid(rpt,:)', squeeze(d.allUnitsTrialsAvoid(rpt,:,:)), d.platDistZReact(rpt,:)', squeeze(d.allUnitsTrialsReact(rpt,:,:)), d.platDistZOpto(rpt,:)', squeeze(d.allUnitsTrialsOpto(rpt,:,:)));
                    disp('train on opto')
                    [dd.rSquaredOpto, dd.rSquaredOptoReact, dd.rSquaredOptoAvoid, dd.lagOpto, dd.lagOptoReact, dd.lagOptoAvoid] = getPlatformDecoderSVR(d.platDistZOpto(rpt,:)', squeeze(d.allUnitsTrialsOpto(rpt,:,:)), d.platDistZReact(rpt,:)', squeeze(d.allUnitsTrialsReact(rpt,:,:)), d.platDistZAvoid(rpt,:)', squeeze(d.allUnitsTrialsAvoid(rpt,:,:)));
                    decodedData.(probeToPlot{probeAnatIdx}).(anatToPlot{probeAnatIdx})(rpt) = dd;
    
                end
            end %end main probe loop
        end %end repeats
        save('decodedDataPlatform.mat', 'decodedData', 'decodedDataAll', 'probeToPlot', 'anatToPlot', 'timescaleSegPlatDec');
    end
    
    if kinematicsData
    % then kinematics data:
        for rpt = 1:numDecoderSampRepeats
            disp(['Processing kinematics sampling repeat ' num2str(rpt) ' of ' num2str(numDecoderSampRepeats)])
            if useNetOrLSTM
                %not used
            else
                d = dataToDecodeKinPlatformData;
%                 mixedPredictors = cat(2, d.kinKptsTrialsReact, d.kinKptsTrialsAvoid, d.kinKptsTrialsOpto);
%                 mixedResponses = cat(2, d.platDistZReact, d.platDistZAvoid, d.platDistZOpto);
%                 concatReactP = cat(2, d.kinKptsTrialsReact, d.kinKptsTrialsReact, d.kinKptsTrialsReact);%to match partitioning indeces with above
%                 concatReactR = cat(2, d.platDistZReact, d.platDistZReact, d.platDistZReact);
%                 concatAvoidP = cat(2, d.kinKptsTrialsAvoid, d.kinKptsTrialsAvoid, d.kinKptsTrialsAvoid);%to match partitioning indeces with above
%                 concatAvoidR = cat(2, d.platDistZAvoid, d.platDistZAvoid, d.platDistZAvoid);
%                 concatOptoP = cat(2, d.kinKptsTrialsOpto, d.kinKptsTrialsOpto, d.kinKptsTrialsOpto);%to match partitioning indeces with above
%                 concatOptoR = cat(2, d.platDistZOpto, d.platDistZOpto, d.platDistZOpto);
%                 [dd.rSquaredMixedMixed, dd.rSquaredMixedAvoid, ~, ~, ~, ~] = getPlatformDecoderSVR(mixedResponses(rpt,:)', squeeze(mixedPredictors(rpt,:,:)), concatAvoidR(rpt,:)', squeeze(concatAvoidP(rpt,:,:)), concatOptoR(rpt,:)', squeeze(concatOptoP(rpt,:,:)));
%                 [~, dd.rSquaredMixedReact, ~, ~, ~, ~] = getPlatformDecoderSVR(mixedResponses(rpt,:)', squeeze(mixedPredictors(rpt,:,:)), concatReactR(rpt,:)', squeeze(concatReactP(rpt,:,:)), concatOptoR(rpt,:)', squeeze(concatOptoP(rpt,:,:)));
%                 [~, dd.rSquaredMixedOpto, ~, ~, ~, ~] = getPlatformDecoderSVR(mixedResponses(rpt,:)', squeeze(mixedPredictors(rpt,:,:)), concatOptoR(rpt,:)', squeeze(concatOptoP(rpt,:,:)), concatOptoR(rpt,:)', squeeze(concatOptoP(rpt,:,:)));
%                 decodedKinDataAll(rpt) = dd; %#ok<*SAGROW> 
%                 clear dd;

                disp('train on react')
                [dd.rSquaredReact, dd.rSquaredReactAvoid, dd.rSquaredReactOpto, dd.lagReact, dd.lagReactAvoid, dd.lagReactOpto] = getPlatformDecoderSVR(d.platDistZReact(rpt,:)', squeeze(d.kinKptsTrialsReact(rpt,:,:)), d.platDistZAvoid(rpt,:)',  squeeze(d.kinKptsTrialsAvoid(rpt,:,:)), d.platDistZOpto(rpt,:)',  squeeze(d.kinKptsTrialsOpto(rpt,:,:)), timescaleSegKinPlatDec);
                disp('train on avoid')
                [dd.rSquaredAvoid, dd.rSquaredAvoidReact, dd.rSquaredAvoidOpto, dd.lagAvoid, dd.lagAvoidReact, dd.lagAvoidOpto] = getPlatformDecoderSVR(d.platDistZAvoid(rpt,:)',  squeeze(d.kinKptsTrialsAvoid(rpt,:,:)), d.platDistZReact(rpt,:)',  squeeze(d.kinKptsTrialsReact(rpt,:,:)), d.platDistZOpto(rpt,:)',  squeeze(d.kinKptsTrialsOpto(rpt,:,:)), timescaleSegKinPlatDec);
                disp('train on opto')
                [dd.rSquaredOpto, dd.rSquaredOptoReact, dd.rSquaredOptoAvoid, dd.lagOpto, dd.lagOptoReact, dd.lagOptoAvoid] = getPlatformDecoderSVR(d.platDistZOpto(rpt,:)',  squeeze(d.kinKptsTrialsOpto(rpt,:,:)), d.platDistZReact(rpt,:)',  squeeze(d.kinKptsTrialsReact(rpt,:,:)), d.platDistZAvoid(rpt,:)',  squeeze(d.kinKptsTrialsAvoid(rpt,:,:)), timescaleSegKinPlatDec);
                decodedKinData(rpt) = dd;

            end
        end %repeats
        save('decodedDataKinPlatform.mat', 'decodedKinData', 'decodedKinDataAll', 'timescaleSegKinPlatDec');
    end

else
    %load saved results & plot
    decodedDataMat = 'decodedDataPlatform.mat';
    load([rootFolder decodedDataMat]);
    decodedDataKinMat = 'decodedDataKinPlatform.mat';
    load([rootFolder decodedDataKinMat]);
    
    % pass anatNames & plot all together:
    plotDecoderAccuracyTrainAll(decodedDataAll, decodedKinDataAll, decodedData, decodedKinData, numDecoderSampRepeats)

%     plotDecoderAccuracyCrossHeatmap(decodedData, decodedKinData, numDecoderSampRepeats)
%     plotDecoderAccuracySameHeatmap(decodedData, decodedKinData, numDecoderSampRepeats) % for use with mixedPredictors etc. above

end

% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function zeroXDottedLine
% for current axes, draw a dotted line at X=0 extending to yLimits
    x = [0 0]; %dotted line at event = time 0
    yLims = get(gca, 'YLim');
    y = [yLims(1) yLims(2)];
    plot(x, y, 'k--', 'LineWidth', 0.1);
end

function [rSquaredSameMean, rSquaredCrossMean1, rSquaredCrossMean2, lagSameMean, lagCrossMean1, lagCrossMean2] = getLinearPlatformDecoder(trainDataPlatform, trainDataPredictors, crossDataPlatform1, crossDataPredictors1, crossDataPlatform2, crossDataPredictors2, timescaleSeg)
% function to get linear decoder weights using pseudoinverse (w/ averaging to prevent overfitting), based on TONIC repository from Josh
%
% behaviorData input is z-scored behavior data (platform distance here), rows are windowed trials around events of interest (extension times)
% predictorData input is z-scored ephys/kinematics data at same windows (time dimensions need to match between behavior and neural/kinematic data)
%
% the basic premise is:
% --------------------
%       behavior = beta*predictor
%       beta = behavior / predictor = behavior * (1/predictor)
% --------------------
% the pseudoinverse of predictor is closest matrix to inverse guaranteed to exist, and corresponds to the least squares solution
% regular least squares can lead to overfitting, so take the average of several batches of trials
% then, multiply computed beta coefficients for each neuron by the spikes to get predicted behavior, and compare to observed behavior (quantify fit using R^2) 

    numTimepoints = length(trainDataPlatform);
    numFolds = 10;
    rSquaredSame = NaN(numFolds,1);
    rSquaredCross1 = NaN(numFolds,1);

    % Define cross-validation options
    cv = cvpartition(numTimepoints,'KFold', numFolds);     

    % loop over folds
    for j = 1:cv.NumTestSets
%         waitbar(j/numFolds, hWait); %update progress bar
        % Split data into training and testing sets
        trainIdx = training(cv,j);
        testIdx = test(cv,j);
        XTrain = trainDataPredictors(trainIdx,:);
        YTrain = trainDataPlatform(trainIdx);

        %train linear model using pseudoinverse linear regression
        betas(:,j) = pinv(XTrain)*YTrain; 

        % save test sets for corss validation in loop below, after taking meanBeta
        XTestSame(j,:,:) = trainDataPredictors(testIdx,:);
        YTestSame(j,:) = trainDataPlatform(testIdx);
        XTestCross1(j,:,:) = crossDataPredictors1(testIdx,:); %just use same random fold indeces
        YTestCross1(j,:) = crossDataPlatform1(testIdx);
        XTestCross2(j,:,:) = crossDataPredictors2(testIdx,:); %just use same random fold indeces
        YTestCross2(j,:) = crossDataPlatform2(testIdx);

    end
    % average betas to get consensus KN dimension
    meanBetas = mean(betas,2);

    % now loop over folds once again to test
    for j = 1:cv.NumTestSets

        % Predict on test sets
        YPredSame = squeeze(XTestSame(j,:,:))*meanBetas;
        YPredCross1 = squeeze(XTestCross1(j,:,:))*meanBetas;
        YPredCross2 = squeeze( XTestCross2(j,:,:))*meanBetas;

        % optionally use xcorr to find best lag before calculating accuracy:
        [rSquaredSame(j), xCorrSame(j), lagSame(j)] = rSquaredXcorr(YTestSame(j,:)',YPredSame);
        %rCorrCoeff = corrcoef(YTestSame,YPredSame)
        [rSquaredCross1(j), xCorrCross1(j), lagCross1(j)] = rSquaredXcorr(YTestCross1(j,:)',YPredCross1);
        [rSquaredCross2(j), xCorrCross2(j), lagCross2(j)] = rSquaredXcorr(YTestCross2(j,:)',YPredCross2);

        % optionally plot overlays %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if j==1
            endSamp = 500;
            timestep = timescaleSeg(2) - timescaleSeg(1);
            timescaleFull = 0:timestep:(timestep*length(YTestSame) - timestep);
            figure; hold on;
            plot(timescaleFull(1:endSamp), YTestSame(j,1:endSamp), 'b-', 'LineWidth', 2)
            plot(timescaleFull(1:endSamp), YPredSame(1:endSamp), 'k--', 'LineWidth', 2)
            title(['same R^2 = ' num2str(rSquaredSame(j)) ', xcorr = ' num2str(xCorrSame(j))])
            hold off;
            legend({'actual','predicted'})
    %         legend({'data same','prediction same', 'data cross','prediction cross'})
            makepretty;

            figure; hold on;
            plot(timescaleFull(1:endSamp), YTestCross1(j,1:endSamp), 'm-', 'LineWidth', 2)
            plot(timescaleFull(1:endSamp), YPredCross1(1:endSamp), 'k--', 'LineWidth', 2)
            title(['cross1 R^2 = ' num2str(rSquaredCross1(j)) ', xcorr = ' num2str(xCorrCross1(j))])
            hold off;
            legend({'actual','predicted'})
    %         legend({'data same','prediction same', 'data cross','prediction cross'})
            makepretty;

            figure; hold on;
            plot(timescaleFull(1:endSamp), YTestCross2(j,1:endSamp), 'g-', 'LineWidth', 2)
            plot(timescaleFull(1:endSamp), YPredCross2(1:endSamp), 'k--', 'LineWidth', 2)
            title(['cross2 R^2 = ' num2str(rSquaredCross2(j)) ', xcorr = ' num2str(xCorrCross2(j))])
            hold off;
            legend({'actual','predicted'})
    %         legend({'data same','prediction same', 'data cross','prediction cross'})
            makepretty;
            keyboard %option to stop and look
        end
    end
%     close(hWait);
        
    rSquaredSameMean = mean(rSquaredSame);
    rSquaredCrossMean1 = mean(rSquaredCross1);
    rSquaredCrossMean2 = mean(rSquaredCross2);
    lagSameMean = mean(lagSame);
    lagCrossMean1 = mean(lagCross1);
    lagCrossMean2 = mean(lagCross2);

end

function [rSquaredSameMean, rSquaredCrossMean1, rSquaredCrossMean2, lagSameMean, lagCrossMean1, lagCrossMean2] = getPlatformDecoderSVR(trainDataPlatform, trainDataPredictors, crossDataPlatform1, crossDataPredictors1, crossDataPlatform2, crossDataPredictors2, timescaleSeg)
% use support vector regression to train model to predict platfrom distance using kinematics or ephys data
% uses 5-fold cross validation to get mean rSquared: folding and training on same set is rSquaredSame output, whereas training on trainData and cross-testing model on testData is rSquaredCross

    numTimepoints = length(trainDataPlatform);
    numFolds = 5;
    rSquaredSame = NaN(numFolds,1);
    rSquaredCross1 = NaN(numFolds,1);

    % Define cross-validation options
    cv = cvpartition(numTimepoints,'KFold', numFolds);     

    % loop over folds
%     hWait = waitbar(1/numFolds,'Creating SVRegression for platform distance...'); %set up progress bar
    % Perform cross-validation
    for j = 1:cv.NumTestSets
%         waitbar(j/numFolds, hWait); %update progress bar
        % Split data into training and testing sets
        trainIdx = training(cv,j);
        testIdx = test(cv,j);
        XTrain = trainDataPredictors(trainIdx,:);
        YTrain = trainDataPlatform(trainIdx);
        XTestSame = trainDataPredictors(testIdx,:);
        YTestSame = trainDataPlatform(testIdx);
        XTestCross1 = crossDataPredictors1(testIdx,:); %just use same random fold indeces
        YTestCross1 = crossDataPlatform1(testIdx);
        XTestCross2 = crossDataPredictors2(testIdx,:); %just use same random fold indeces
        YTestCross2 = crossDataPlatform2(testIdx);

        %train SVR with nonlinear kernel
%         rng(1); %specify random number seed to reproducible results with KernelScale=auto below
%         tic
        SVMModel = fitrsvm(XTrain, YTrain,'KernelFunction','gaussian','Standardize', false, 'BoxConstraint',1,'KernelScale',1,'CacheSize','maximal','Verbose', 0,'NumPrint',1000); %set Verbose=1 to print convergence & debug
%         toc
        % BoxConstraint is regularization parameter C to prevent overfitting, but trades off with computation time

        % Predict on test set
        YPredSame = predict(SVMModel, XTestSame);
        YPredCross1 = predict(SVMModel, XTestCross1);
        YPredCross2 = predict(SVMModel, XTestCross2);

        % optionally check data:
%         figure; hold on;
%         plot(YPredSame, 'k--', 'LineWidth', 2)
%         plot(YTestSame, 'b--', 'LineWidth', 2)
%         plot(XTestSame(:,1))
%         plot(XTestSame(:,2))
%         plot(XTestSame(:,3))
%         plot(XTestSame(:,4))
%         plot(XTestSame(:,5))
        
        % Calculate accuracy
%         rSquaredSame(j) = rSquared(YTestSame,YPredSame);
%         %rCorrCoeff = corrcoef(YTestSame,YPredSame)
%         rSquaredCross1(j) = rSquared(YTestCross1,YPredCross1);
%         rSquaredCross2(j) = rSquared(YTestCross2,YPredCross2);

        % optionally use xcorr to find best lag before calculating accuracy:
        [rSquaredSame(j), xCorrSame(j), lagSame(j)] = rSquaredXcorr(YTestSame,YPredSame);
        %rCorrCoeff = corrcoef(YTestSame,YPredSame)
        [rSquaredCross1(j), xCorrCross1(j), lagCross1(j)] = rSquaredXcorr(YTestCross1,YPredCross1);
        [rSquaredCross2(j), xCorrCross2(j), lagCross2(j)] = rSquaredXcorr(YTestCross2,YPredCross2);

        % optionally plot overlays %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if j==1
            timestep = 0.02; %20ms bins
            timescaleFull = 0:timestep:(timestep*length(YTestSame) - timestep);
            figure; hold on;
            plot(timescaleFull(1:500), YTestSame(1:500), 'b-', 'LineWidth', 2)
            plot(timescaleFull(1:500), YPredSame(1:500), 'k--', 'LineWidth', 2)
            title(['same R^2 = ' num2str(rSquaredSame(j)) ', xcorr = ' num2str(xCorrSame(j))])
            hold off;
            legend({'actual','predicted'})
    %         legend({'data same','prediction same', 'data cross','prediction cross'})
            makepretty;

            figure; hold on;
            plot(timescaleFull(1:500), YTestCross1(1:500), 'm-', 'LineWidth', 2)
            plot(timescaleFull(1:500), YPredCross1(1:500), 'k--', 'LineWidth', 2)
            title(['cross1 R^2 = ' num2str(rSquaredCross1(j)) ', xcorr = ' num2str(xCorrCross1(j))])
            hold off;
            legend({'actual','predicted'})
    %         legend({'data same','prediction same', 'data cross','prediction cross'})
            makepretty;

            figure; hold on;
            plot(timescaleFull(1:500), YTestCross2(1:500), 'g-', 'LineWidth', 2)
            plot(timescaleFull(1:500), YPredCross2(1:500), 'k--', 'LineWidth', 2)
            title(['cross2 R^2 = ' num2str(rSquaredCross2(j)) ', xcorr = ' num2str(xCorrCross2(j))])
            hold off;
            legend({'actual','predicted'})
    %         legend({'data same','prediction same', 'data cross','prediction cross'})
            makepretty;
            keyboard %option to stop and look
        end
    end
%     close(hWait);
        
    rSquaredSameMean = mean(rSquaredSame);
    rSquaredCrossMean1 = mean(rSquaredCross1);
    rSquaredCrossMean2 = mean(rSquaredCross2);
    lagSameMean = mean(lagSame);
    lagCrossMean1 = mean(lagCross1);
    lagCrossMean2 = mean(lagCross2);
    
end

function plotDecoderAccuracyTrainAll(decodedDataTrainedAll, decodedKinDataTrainedAll, decodedData, decodedKinData, numDecoderSampRepeats)
 % plot accuracies where training data corresponds to all data mixed

    probeToPlot = {'Probe1','Probe0','Probe0','Probe2'}; %matched pairs with below
    anatToPlot =  {'CTX',   'CTX',   'BS',    'HB'};
    labels =      {'kinematics', 'PFC',    'M1',   'TH/HY', 'HB'}; 
    lineCmaps = linspecer(size(probeToPlot,2));
    anatCmaps =   [lineCmaps(2,:); lineCmaps(1,:); lineCmaps(3,:); lineCmaps(4,:)]; %reorder to put PFC first
    blueCmap = [0 0.4470 0.7410];
    purpleCmap = [0.4940 0.1840 0.5560]; %purple for kinematics
    orangeCmap = [0.8500 0.3250 0.0980];
    cmapKin = [1 0 1];

    f = figure;
    tiledlayout(1,4,"TileSpacing","compact")

    % first collect data over repeats
    for rpt = 1:numDecoderSampRepeats
        for i = 1:size(probeToPlot,2) %main loop over probe structures

            if i == 1 % kinematics once only
                    dataAKin(rpt) = decodedKinDataTrainedAll(rpt).rSquaredMixedAvoid;
                    dataRKin(rpt) = decodedKinDataTrainedAll(rpt).rSquaredMixedReact;
                    dataOKin(rpt) = decodedKinDataTrainedAll(rpt).rSquaredMixedOpto;
            end
            dataA(i,rpt) = decodedDataTrainedAll(rpt).(probeToPlot{1,i}).(anatToPlot{1,i}).rSquaredMixedAvoid;
            dataR(i,rpt) = decodedDataTrainedAll(rpt).(probeToPlot{1,i}).(anatToPlot{1,i}).rSquaredMixedReact;
            dataO(i,rpt) = decodedDataTrainedAll(rpt).(probeToPlot{1,i}).(anatToPlot{1,i}).rSquaredMixedOpto;

        end
    end

    columnData = [];
    % first plot trained on ALL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%:
    nexttile; hold on;
    for i = 1:size(probeToPlot,2)

        if i == 1 %plot kinematics first
            x1 = ones(size(dataAKin,1)) - 0.3;
            x2 = ones(size(dataAKin,1));
            x3 = ones(size(dataAKin,1)) + 0.3;
            h1k(i) = bar(mean(x1),mean(dataAKin),0.2,'FaceColor',purpleCmap,'FaceAlpha',0.2,'EdgeColor',purpleCmap);
            h2k(i) = bar(mean(x2),mean(dataRKin),0.2,'w','EdgeColor',purpleCmap); 
            h3k(i) = bar(mean(x3),mean(dataOKin),0.2,'w','EdgeColor',purpleCmap);
            hatchfill2(h3k(i),'single','HatchAngle',45,'HatchColor',purpleCmap);
            swarmchart(x1,dataAKin, 50, purpleCmap, 'filled','MarkerFaceAlpha',0.5, 'XJitter', 'rand', 'XJitterWidth', 1);
            swarmchart(x2,dataRKin, 50, purpleCmap, 'filled','MarkerFaceAlpha',0.5, 'XJitter', 'rand', 'XJitterWidth', 1);
            swarmchart(x3,dataOKin, 50, purpleCmap, 'filled','MarkerFaceAlpha',0.5, 'XJitter', 'rand', 'XJitterWidth', 1);
            columnData = [dataAKin' dataRKin' dataOKin'];
        end
        
        x1 = ones(size(dataA,2),1) + i*2 - 0.3;
        x2 = ones(size(dataA,2),1) + i*2;
        x3 = ones(size(dataA,2),1) + i*2 + 0.3;
        h1(i) = bar(mean(x1),mean(dataA(i,:)),0.2,'FaceColor',anatCmaps(i,:),'FaceAlpha',0.2,'EdgeColor',anatCmaps(i,:));
        h2(i) = bar(mean(x2),mean(dataR(i,:)),0.2,'w','EdgeColor',anatCmaps(i,:));
        h3(i) = bar(mean(x3),mean(dataO(i,:)),0.2,'w','EdgeColor',anatCmaps(i,:));
        hatchfill2(h3(i),'single','HatchAngle',45,'HatchColor',anatCmaps(i,:));
        swarmchart(x1,dataA(i,:)', 50, anatCmaps(i,:), 'filled','MarkerFaceAlpha',0.5, 'XJitter', 'rand', 'XJitterWidth', 0.2);
        swarmchart(x2,dataR(i,:)', 50, anatCmaps(i,:), 'filled','MarkerFaceAlpha',0.5, 'XJitter', 'rand', 'XJitterWidth', 0.2);
        swarmchart(x3,dataO(i,:)', 50, anatCmaps(i,:), 'filled','MarkerFaceAlpha',0.5, 'XJitter', 'rand', 'XJitterWidth', 0.2);
        columnData = [columnData dataA(i,:)' dataR(i,:)' dataO(i,:)'];
    end 
    ylabel('extension distance decoder accuracy');
    ylim([0 1])
    xticks([1:2:10])
    xticklabels(labels)
    hold off;
    makepretty;

    % now plot trained on avoid /react / opto the same way %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%:
    for avoidReactOpto = 1:3
        nexttile; hold on;
        % first collect data over repeats
        for rpt = 1:numDecoderSampRepeats
            for i = 1:size(probeToPlot,2) %main loop over probe structures
    
                if i == 1 % kinematics once only
                    if avoidReactOpto==1
                        dataAKin(rpt) = decodedKinData(rpt).rSquaredAvoid;
                        dataRKin(rpt) = decodedKinData(rpt).rSquaredAvoidReact;
                        dataOKin(rpt) = decodedKinData(rpt).rSquaredAvoidOpto;
                    elseif avoidReactOpto==2
                        dataAKin(rpt) = decodedKinData(rpt).rSquaredReactAvoid;
                        dataRKin(rpt) = decodedKinData(rpt).rSquaredReact;
                        dataOKin(rpt) = decodedKinData(rpt).rSquaredReactOpto;
                    elseif avoidReactOpto==3
                        dataAKin(rpt) = decodedKinData(rpt).rSquaredOptoAvoid;
                        dataRKin(rpt) = decodedKinData(rpt).rSquaredOptoReact;
                        dataOKin(rpt) = decodedKinData(rpt).rSquaredOpto;
                    end
                end
    
                if avoidReactOpto==1
                    dataA(i,rpt) = decodedData(rpt).(probeToPlot{1,i}).(anatToPlot{1,i}).rSquaredAvoid;
                    dataR(i,rpt) = decodedData(rpt).(probeToPlot{1,i}).(anatToPlot{1,i}).rSquaredAvoidReact;
                    dataO(i,rpt) = decodedData(rpt).(probeToPlot{1,i}).(anatToPlot{1,i}).rSquaredAvoidOpto;
                    titleStr = 'trained on avoid trials';
                elseif avoidReactOpto==2
                    dataA(i,rpt) = decodedData(rpt).(probeToPlot{1,i}).(anatToPlot{1,i}).rSquaredReactAvoid;
                    dataR(i,rpt) = decodedData(rpt).(probeToPlot{1,i}).(anatToPlot{1,i}).rSquaredReact;
                    dataO(i,rpt) = decodedData(rpt).(probeToPlot{1,i}).(anatToPlot{1,i}).rSquaredReactOpto;
                    titleStr = 'trained on react trials';
                elseif avoidReactOpto==3
                    dataA(i,rpt) = decodedData(rpt).(probeToPlot{1,i}).(anatToPlot{1,i}).rSquaredOptoAvoid;
                    dataR(i,rpt) = decodedData(rpt).(probeToPlot{1,i}).(anatToPlot{1,i}).rSquaredOptoReact;
                    dataO(i,rpt) = decodedData(rpt).(probeToPlot{1,i}).(anatToPlot{1,i}).rSquaredOpto;
                    titleStr = 'trained on opto trials';
                else
                    disp('bad avoidReactOpto input')
                    return
                end
    
            end
        end
    
        columnData = [];
        for i = 1:size(probeToPlot,2)
            if i == 1 %plot kinematics first
                x1 = ones(size(dataAKin,1)) - 0.3;
                x2 = ones(size(dataAKin,1));
                x3 = ones(size(dataAKin,1)) + 0.3;
                h1k(i) = bar(mean(x1),mean(dataAKin),0.2,'FaceColor',purpleCmap,'FaceAlpha',0.2,'EdgeColor',purpleCmap);
                h2k(i) = bar(mean(x2),mean(dataRKin),0.2,'w','EdgeColor',purpleCmap); 
                h3k(i) = bar(mean(x3),mean(dataOKin),0.2,'w','EdgeColor',purpleCmap);
                hatchfill2(h3k(i),'single','HatchAngle',45,'HatchColor',purpleCmap);
                swarmchart(x1,dataAKin, 50, purpleCmap, 'filled','MarkerFaceAlpha',0.5, 'XJitter', 'rand', 'XJitterWidth', 1);
                swarmchart(x2,dataRKin, 50, purpleCmap, 'filled','MarkerFaceAlpha',0.5, 'XJitter', 'rand', 'XJitterWidth', 1);
                swarmchart(x3,dataOKin, 50, purpleCmap, 'filled','MarkerFaceAlpha',0.5, 'XJitter', 'rand', 'XJitterWidth', 1);
                columnData = [dataAKin' dataRKin' dataOKin'];
            end
            
            x1 = ones(size(dataA,2),1) + i*2 - 0.3;
            x2 = ones(size(dataA,2),1) + i*2;
            x3 = ones(size(dataA,2),1) + i*2 + 0.3;
            h1(i) = bar(mean(x1),mean(dataA(i,:)),0.2,'FaceColor',anatCmaps(i,:),'FaceAlpha',0.2,'EdgeColor',anatCmaps(i,:));
            h2(i) = bar(mean(x2),mean(dataR(i,:)),0.2,'w','EdgeColor',anatCmaps(i,:));
            h3(i) = bar(mean(x3),mean(dataO(i,:)),0.2,'w','EdgeColor',anatCmaps(i,:));
            hatchfill2(h3(i),'single','HatchAngle',45,'HatchColor',anatCmaps(i,:));
            swarmchart(x1,dataA(i,:)', 50, anatCmaps(i,:), 'filled','MarkerFaceAlpha',0.5, 'XJitter', 'rand', 'XJitterWidth', 0.2);
            swarmchart(x2,dataR(i,:)', 50, anatCmaps(i,:), 'filled','MarkerFaceAlpha',0.5, 'XJitter', 'rand', 'XJitterWidth', 0.2);
            swarmchart(x3,dataO(i,:)', 50, anatCmaps(i,:), 'filled','MarkerFaceAlpha',0.5, 'XJitter', 'rand', 'XJitterWidth', 0.2);
            columnData = [columnData dataA(i,:)' dataR(i,:)' dataO(i,:)'];
        end 
        ylim([0 1])
        set(gca,'ytick',[]);
        xticks([1:2:10])
        xticklabels(labels)
        hold off;
        makepretty;

    end

end

function plotDecoderAccuracyCrossHeatmap(decodedData, decodedKinData, numDecoderSampRepeats)
 % plot cross-decoding accuracies heatmaps

    probeToPlot = {'Probe1','Probe0','Probe0','Probe2'}; %matched pairs with below
    anatToPlot =  {'CTX',   'CTX',   'BS',    'HB'};
%     labels =      {'kinematics', 'PFC',    'M1',   'TH/HY', 'HB'}; 
    labels =      {'PFC',    'M1',   'TH/HY', 'HB'}; 


    % first collect mean data over repeats
    for i = 1:size(probeToPlot,2) %main loop over probe structures
        for rpt = 1:numDecoderSampRepeats
        
%             if i == 1 %plot kinematics first
%                 dataAKin(rpt) = decodedKinData(rpt).rSquaredReactAvoid;
%                 dataRKin(rpt) = decodedKinData(rpt).rSquaredReact;
%                 dataOKin(rpt) = decodedKinData(rpt).rSquaredReactOpto;
%             end
            
            %trained on avoid
            dataAA(i,rpt) = decodedData.(probeToPlot{1,i}).(anatToPlot{1,i})(rpt).rSquaredAvoid;
            dataAR(i,rpt) = decodedData.(probeToPlot{1,i}).(anatToPlot{1,i})(rpt).rSquaredAvoidReact;
            dataAO(i,rpt) = decodedData.(probeToPlot{1,i}).(anatToPlot{1,i})(rpt).rSquaredAvoidOpto;
            %trained on react
            dataRA(i,rpt) = decodedData.(probeToPlot{1,i}).(anatToPlot{1,i})(rpt).rSquaredReactAvoid;
            dataRR(i,rpt) = decodedData.(probeToPlot{1,i}).(anatToPlot{1,i})(rpt).rSquaredReact;
            dataRO(i,rpt) = decodedData.(probeToPlot{1,i}).(anatToPlot{1,i})(rpt).rSquaredReactOpto;
            %trained on opto
            dataOA(i,rpt) = decodedData.(probeToPlot{1,i}).(anatToPlot{1,i})(rpt).rSquaredOptoAvoid;
            dataOR(i,rpt) = decodedData.(probeToPlot{1,i}).(anatToPlot{1,i})(rpt).rSquaredOptoReact;
            dataOO(i,rpt) = decodedData.(probeToPlot{1,i}).(anatToPlot{1,i})(rpt).rSquaredOpto;

        end
        %take mean over repeats:
        if numDecoderSampRepeats==1
            dataAAmean(i) = dataAA(i);
            dataARmean(i) = dataAR(i);
            dataAOmean(i) = dataAO(i);
            dataRAmean(i) = dataRA(i);
            dataRRmean(i) = dataRR(i);
            dataROmean(i) = dataRO(i);
            dataOAmean(i) = dataOA(i);
            dataORmean(i) = dataOR(i);
            dataOOmean(i) = dataOO(i);
        else
            dataAAmean(i) = mean(dataAA(i,:),2);
            dataARmean(i) = mean(dataAR(i,:),2);
            dataAOmean(i) = mean(dataAO(i,:),2);
            dataRAmean(i) = mean(dataRA(i,:),2);
            dataRRmean(i) = mean(dataRR(i,:),2);
            dataROmean(i) = mean(dataRO(i,:),2);
            dataOAmean(i) = mean(dataOA(i,:),2);
            dataORmean(i) = mean(dataOR(i,:),2);
            dataOOmean(i) = mean(dataOO(i,:),2);
        end
    end



    f = figure; 


    heatData = [];
    % then plot
    for i = 1:size(probeToPlot,2)
%         if i == 1 %plot kinematics first
%             heatData = [dataAKin' dataRKin' dataOKin'];
%         end
        
        nexttile; hold on;
        heatData = [dataAAmean(i) dataARmean(i) dataAOmean(i);...
                dataRAmean(i) dataRRmean(i) dataROmean(i);...
                dataOAmean(i) dataORmean(i) dataOOmean(i)];
        imagesc(heatData)
        set(gca,'YDir','reverse')
        colormap("parula")
        xticks([1 2 3])
        yticks([1 2 3])
        xticklabels({'avoid','react','opto'});
        yticklabels({'opto','react','avoid'});
        ylabel('trained on');
        xlabel('tested on');
        clim([0 1])
        colorbar
        title(labels{1,i})
        hold off;
        makepretty;

    end 


    %check statistics:
%     [pFriedman,tblFriedman,statsFriedman] = friedman(columnData,1); %overall p-value against null that all columns same
%     figure
%     cFriedman = multcompare(statsFriedman, 'estimate', 'friedman', 'ctype', 'dunn-sidak'); %table of p-values across pairs

end

function plotDecoderAccuracySameHeatmap(decodedData, decodedKinData, numDecoderSampRepeats)
 % plot decoding accuracy heatmaps when trained on all types of trials

    probeToPlot = {'Probe1','Probe0','Probe0','Probe2'}; %matched pairs with below
    anatToPlot =  {'CTX',   'CTX',   'BS',    'HB'};
%     labels =      {'kinematics', 'PFC',    'M1',   'TH/HY', 'HB'}; 
    labels =      {'PFC',    'M1',   'TH/HY', 'HB'}; 


    % first collect mean data over repeats
    for i = 1:size(probeToPlot,2) %main loop over probe structures
        for rpt = 1:numDecoderSampRepeats
        
%             if i == 1 %plot kinematics first
%                 dataAKin(rpt) = decodedKinData(rpt).rSquaredReactAvoid;
%                 dataRKin(rpt) = decodedKinData(rpt).rSquaredReact;
%                 dataOKin(rpt) = decodedKinData(rpt).rSquaredReactOpto;
%             end
            
            %trained on all
%             dataA(i,rpt) = decodedData.(probeToPlot{1,i}).(anatToPlot{1,i})(rpt).rSquaredMixedAvoid;
%             dataR(i,rpt) = decodedData.(probeToPlot{1,i}).(anatToPlot{1,i})(rpt).rSquaredMixedReact;
%             dataO(i,rpt) = decodedData.(probeToPlot{1,i}).(anatToPlot{1,i})(rpt).rSquaredMixedOpto;

            dataA(i,rpt) = decodedData.(probeToPlot{1,i}).(anatToPlot{1,i})(rpt).rSquaredReactAvoid;
            dataR(i,rpt) = decodedData.(probeToPlot{1,i}).(anatToPlot{1,i})(rpt).rSquaredAvoidReact;
            dataO(i,rpt) = decodedData.(probeToPlot{1,i}).(anatToPlot{1,i})(rpt).rSquaredOptoReact;


        end
        %take mean over repeats:
        if numDecoderSampRepeats==1
            dataAmean(i) = dataA(i);
            dataRmean(i) = dataR(i);
            dataOmean(i) = dataO(i);
        else
            dataAmean(i) = mean(dataA(i,:),2);
            dataRmean(i) = mean(dataR(i,:),2);
            dataOmean(i) = mean(dataO(i,:),2);
        end
    end

    f = figure; 
    heatData = [];
    % then plot
    for i = 1:size(probeToPlot,2)
%         if i == 1 %plot kinematics first
%             heatData = [dataAKin' dataRKin' dataOKin'];
%         end
        
        nexttile; hold on;
        heatData = [dataAmean(i) dataRmean(i) dataOmean(i)];
        imagesc(heatData)
        set(gca,'YDir','reverse')
        colormap("parula")
        xticks([1 2 3])
        xticklabels({'avoid','react','opto'});
        ylabel('trained on all');
        xlabel('tested on');
        clim([0 1])
        colorbar
        title(labels{1,i})
        hold off;
        makepretty;

    end 


    %check statistics:
%     [pFriedman,tblFriedman,statsFriedman] = friedman(columnData,1); %overall p-value against null that all columns same
%     figure
%     cFriedman = multcompare(statsFriedman, 'estimate', 'friedman', 'ctype', 'dunn-sidak'); %table of p-values across pairs

end

function plotDecoderLag(decodedData, decodedKinData, numDecoderSampRepeats,avoidReactOpto)
 % plot accuracies where training data corresponds to
 % combined ctx or non-ctx units avoid vs react mean responses
  % 'avoidReactOpto' input is whether to plot data when trained on avoid=1, trained on react=2, trained on opto=3


    probeToPlot = {'Probe0','Probe1','Probe0','Probe2'}; %matched pairs with below
    anatToPlot =  {'CTX',   'CTX',   'BS',    'HB'};
    labels =      {'M1',    'PFC',   'TH/HY', 'HB'}; 
    blueCmap = [0 0.4470 0.7410];
    orangeCmap = [0.8500 0.3250 0.0980];
    lineCmaps = linspecer(size(probeToPlot,2));
    cmapKin = [1 0 1];

    f = figure; hold on;

    % first collect data over repeats
    for rpt = 1:numDecoderSampRepeats
        for i = 1:size(probeToPlot,2) %main loop over probe structures

%             if i == 1 %plot kinematics first
%                 dataAKin(rpt) = decodedKinData(rpt).rSquaredReactAvoid;
%                 dataRKin(rpt) = decodedKinData(rpt).rSquaredReact;
%                 dataOKin(rpt) = decodedKinData(rpt).rSquaredReactOpto;
%             end

            if avoidReactOpto==1
                dataA(i,rpt) = decodedData.(probeToPlot{1,i}).(anatToPlot{1,i})(rpt).lagAvoid*1000; %convert to milleseconds;
                dataR(i,rpt) = decodedData.(probeToPlot{1,i}).(anatToPlot{1,i})(rpt).lagAvoidReact*1000;
                dataO(i,rpt) = decodedData.(probeToPlot{1,i}).(anatToPlot{1,i})(rpt).lagAvoidOpto*1000;
                titleStr = 'lag trained on avoid trials';
            elseif avoidReactOpto==2
                dataA(i,rpt) = decodedData.(probeToPlot{1,i}).(anatToPlot{1,i})(rpt).lagReactAvoid*1000;
                dataR(i,rpt) = decodedData.(probeToPlot{1,i}).(anatToPlot{1,i})(rpt).lagReact*1000;
                dataO(i,rpt) = decodedData.(probeToPlot{1,i}).(anatToPlot{1,i})(rpt).lagReactOpto*1000;
                titleStr = 'lag trained on react trials';
            elseif avoidReactOpto==3
                dataA(i,rpt) = decodedData.(probeToPlot{1,i}).(anatToPlot{1,i})(rpt).lagOptoAvoid*1000;
                dataR(i,rpt) = decodedData.(probeToPlot{1,i}).(anatToPlot{1,i})(rpt).lagOptoReact*1000;
                dataO(i,rpt) = decodedData.(probeToPlot{1,i}).(anatToPlot{1,i})(rpt).lagOpto*1000;
                titleStr = 'lag trained on opto trials';
            else
                disp('bad avoidReactOpto input')
                return
            end

        end
    end

    columnData = [];
    % then plot
    for i = 1:size(probeToPlot,2)
%         if i == 1 %plot kinematics first
%             x1 = ones(size(dataAKin,1)) - 0.3;
%             x2 = ones(size(dataAKin,1));
%             x3 = ones(size(dataAKin,1)) + 0.3;
%             bar(mean(x1),mean(dataAKin),0.2,'w','EdgeColor',orangeCmap);
%             bar(mean(x2),mean(dataRKin),0.2,'w','EdgeColor',blueCmap);
%             bar(mean(x3),mean(dataOKin),0.2,'w','EdgeColor',[0 0 0]);
%             swarmchart(x1,dataAKin, 100, orangeCmap, 'filled','MarkerFaceAlpha',0.4, 'XJitter', 'rand', 'XJitterWidth', 1);
%             swarmchart(x2,dataRKin, 100, blueCmap, 'filled','MarkerFaceAlpha',0.4, 'XJitter', 'rand', 'XJitterWidth', 1);
%             swarmchart(x3,dataOKin, 100, [0 0 0], 'filled','MarkerFaceAlpha',0.4, 'XJitter', 'rand', 'XJitterWidth', 1);
%             columnData = [dataAKin' dataRKin' dataOKin'];
% 
%         end
        
        x1 = ones(size(dataA,2),1) + i*2 - 0.3;
        x2 = ones(size(dataA,2),1) + i*2;
        x3 = ones(size(dataA,2),1) + i*2 + 0.3;
        bar(mean(x1),mean(dataA(i,:)),0.2,'w','EdgeColor',orangeCmap);
        bar(mean(x2),mean(dataR(i,:)),0.2,'w','EdgeColor',blueCmap);
        bar(mean(x3),mean(dataO(i,:)),0.2,'w','EdgeColor',[0 0 0]);
        swarmchart(x1,dataA(i,:)', 100, orangeCmap, 'filled','MarkerFaceAlpha',0.4, 'XJitter', 'rand', 'XJitterWidth', 0.2);
        swarmchart(x2,dataR(i,:)', 100, blueCmap, 'filled','MarkerFaceAlpha',0.4, 'XJitter', 'rand', 'XJitterWidth', 0.2);
        swarmchart(x3,dataO(i,:)', 100, [0 0 0], 'filled','MarkerFaceAlpha',0.4, 'XJitter', 'rand', 'XJitterWidth', 0.2);
        columnData = [columnData dataA(i,:)' dataR(i,:)' dataO(i,:)'];

    end 
    ylabel('extension distance decoder lag (msec)');
%     ylim([0 1])
    xticks([3 5 7 9])
    xticklabels(labels)
    title(titleStr)
    hold off;
    makepretty;

    %check statistics:
%     [pFriedman,tblFriedman,statsFriedman] = friedman(columnData,1); %overall p-value against null that all columns same
%     figure
%     cFriedman = multcompare(statsFriedman, 'estimate', 'friedman', 'ctype', 'dunn-sidak'); %table of p-values across pairs

end

function [coeffDet] = rSquared(data, prediction)
% Coefficient of determination (R-squared) indicates the proportionate amount of variation in the data explained by the model prediction
% The larger the R-squared is, the more variability is explained by the linear regression model.

    % Sum of squared residuals
    SSR = sum((prediction - data).^2);
    % Total sum of squares
    SST = sum(((data - mean(data)).^2));
    % R squared
    coeffDet = 1 - SSR/SST;

end

function [maxRsquared, maxXcorr, lagAtMax] = rSquaredXcorr(data, prediction)
% instead of zero-lag R-squared, calculate best xcorr lag first and then do coeff of determination

    numLags = 20; %make reasonable range here; with 20ms bins, 10 lags = 200ms
    timeBin = 0.02;
%     [r,~] = xcorr(data, prediction, numLags,'biased'); %biased normalization just divides by # elements
    r =  myXcorr(data, prediction, numLags);
    bestCorrLagIdx = find(r == max(r), 1, 'first') - numLags - 1; %
    lagAtMax = bestCorrLagIdx*timeBin; %negative lag means prediction leads, positive lag means data leads

    % option to return simple Pearson correlation coefficient:
%     [R,P] = corrcoef(data,prediction);
%     maxRsquared = R(2);

    % option to return simple normalized cross correlation value:
    maxXcorr = max(r);

    % also get coefficient of determination:
    maxRsquared = rSquared(data, prediction);

    % first check if correlation is even positive; else shifting won't really help
%     origRsquared = rSquared(data, prediction);
%     if origRsquared < 0
%         maxRsquared = origRsquared;
%         lagAtMax = NaN;
%     else
%         % if correlation is reasonable, then try xcorr to find the best shift
%         if bestCorrLagIdx == 0
%             maxRsquared = rSquared(data, prediction);
%         elseif bestCorrLagIdx > 0
%             dataShifted = circshift(data,bestCorrLagIdx);  %positive shift to right, forward in time
%             dataShiftedTruncated = dataShifted(1+bestCorrLagIdx:end);
%             predictionTruncated = prediction(1+bestCorrLagIdx:end);
%             maxRsquared = rSquared(dataShiftedTruncated, predictionTruncated);
%         elseif bestCorrLagIdx < 0
%             dataShifted = circshift(data,bestCorrLagIdx); %negative shift to left, backwards in time
%             dataShiftedTruncated = dataShifted(1:end+bestCorrLagIdx);
%             predictionTruncated = prediction(1:end+bestCorrLagIdx);
%             maxRsquared = rSquared(dataShiftedTruncated, predictionTruncated);
%         end
%     end

end

function r = myXcorr(X,Y,maxLag)
% normalized by vector norms

    % Preallocate array to store normalized cross-correlation values
    r = zeros(1, 2 * maxLag + 1);
    
    % Loop over each lag from -maxLag to +maxLag
    for lag = -maxLag:maxLag
        if lag < 0
            % Lag Y negatively (shift Y to the left)
            X_shifted = X(1:end+lag);
            Y_shifted = Y(1-lag:end);
        elseif lag > 0
            % Lag Y positively (shift Y to the right)
            X_shifted = X(1+lag:end);
            Y_shifted = Y(1:end-lag);
        else
            % No lag
            X_shifted = X;
            Y_shifted = Y;
        end
        
        % Compute the dot product of shifted X and Y
        dotProduct = dot(X_shifted, Y_shifted);
        
        % Compute the Euclidean norms of the shifted vectors
        normX = norm(X_shifted);
        normY = norm(Y_shifted);
        
        % Compute the normalized cross-correlation for the current lag
        r(lag + maxLag + 1) = dotProduct / (normX * normY);
    end
    
end
