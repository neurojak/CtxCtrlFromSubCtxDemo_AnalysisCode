% Jason Keller
% 2024
% load saved decoder data in structure and run decoders, save output in MAT file
% note that classifiers are always run from shared subspace PCs, so that we can combine trials across mice

rootFolder = 'C:\';
decodeAndSave = 0; % 1=train, 0=plot results
kinematicsData = 1;
ephysData = 1;
useNetOrLSTM = 0; % (1) neural net / LSTM or (0) SVM; neural net option not currently used
numSampRepeats = 10;
fontSz = 16;

if decodeAndSave
    decoderMatNameEphys = 'dataToClassify.mat';
    decoderMatNameKinematics = 'kinDataToClassify.mat';
    
    if ephysData
        load([rootFolder decoderMatNameEphys]);
        numDecTimepoints = length(timescaleSegClassify); % trialData is one longer
    end
    if kinematicsData
        load([rootFolder decoderMatNameKinematics]);
        numDecKinTimepoints = length(timescaleSegKinClassify);
    end

    if ephysData
        % first run ephys data:
        for rpt = 1:numSampRepeats
            disp(['Processing  ephys sampling repeat ' num2str(rpt) ' of ' num2str(numSampRepeats)])
            for probeAnatIdx = 1:length(probesSub) %main loop over probe structures; ex. {'Probe1'} 
                disp(['Processing ' labelsSub{probeAnatIdx}])

                trialLabels = dataToClassify.(probesSub{probeAnatIdx}).(anatsSub{probeAnatIdx}).classifierTrialCategories{rpt};
                trialLabelsShuffled = dataToClassify.(probesSub{probeAnatIdx}).(anatsSub{probeAnatIdx}).classifierTrialCategoriesShuf{rpt};
                trialData = dataToClassify.(probesSub{probeAnatIdx}).(anatsSub{probeAnatIdx}).classifierTrialData{rpt};

                if useNetOrLSTM
                    classifiedData(rpt).(probesSub{probeAnatIdx}).(anatsSub{probeAnatIdx}).classifierAccuracyOverTime = getNonlinearTrialTimeClassifierLSTM(trialLabels, trialData);
                    classifiedData(rpt).(probesSub{probeAnatIdx}).(anatsSub{probeAnatIdx}).classifierAccuracyOverTimeShuffle = getNonlinearTrialTimeClassifierLSTM(trialLabelsShuffled, trialData);
                    %classifiedData(rpt).(probesSub{probeAnatIdx}).(anatsSub{probeAnatIdx}).classifierAccuracyOverTime = getNonlinearTrialTimeClassifierNet(trialLabels, trialData); %#ok<*UNRCH> 
                    %classifiedData(rpt).(probesSub{probeAnatIdx}).(anatsSub{probeAnatIdx}).classifierAccuracyOverTimeShuffle = getNonlinearTrialTimeClassifierNet(trialLabelsShuffled, trialData);                               
                else
                    classifiedData(rpt).(probesSub{probeAnatIdx}).(anatsSub{probeAnatIdx}).classifierAccuracyOverTime = getNonlinearTrialTimeClassifierSVM(trialLabels, trialData);
                    classifiedData(rpt).(probesSub{probeAnatIdx}).(anatsSub{probeAnatIdx}).classifierAccuracyOverTimeShuffle = getNonlinearTrialTimeClassifierSVM(trialLabelsShuffled, trialData);
                end


            end %end main current probe loop
        end %end rpt
        save('classifiedData.mat', 'classifiedData', 'timescaleSegClassify');
    end
    
    if kinematicsData
        % then kinematics data:
        for rpt = 1:numSampRepeats
            disp(['Processing  kinematics sampling repeat ' num2str(rpt) ' of ' num2str(numSampRepeats)])
            numMice = size(dataToClassifyKinData.trialLabels,2);

                trialLabels = dataToClassifyKinData.trialLabels{rpt};
                trialLabelsShuffled = dataToClassifyKinData.trialLabelsShuffled{rpt};
                trialData = dataToClassifyKinData.trialData{rpt};

                if useNetOrLSTM
                    classifiedKinData(rpt).classifierAccuracyOverTime = getNonlinearTrialTimeClassifierLSTM(trialLabels, trialData);
                    classifiedKinData(rpt).classifierAccuracyOverTimeShuffle = getNonlinearTrialTimeClassifierLSTM(trialLabelsShuffled, trialData);
                    %classifiedKinData(rpt).classifierAccuracyOverTime(k,:) = getNonlinearTrialTimeClassifierNet(trialLabels, trialData);
                    %classifiedKinData(rpt).classifierAccuracyOverTimeShuffle(k,:) = getNonlinearTrialTimeClassifierNet(trialLabelsShuffled, trialData);
                else
                    classifiedKinData(rpt).classifierAccuracyOverTime = getNonlinearTrialTimeClassifierSVM(trialLabels, trialData);
                    classifiedKinData(rpt).classifierAccuracyOverTimeShuffle = getNonlinearTrialTimeClassifierSVM(trialLabelsShuffled, trialData);
                end
        end
        save('classifiedDataKin.mat', 'classifiedKinData', 'timescaleSegKinClassify');
    end

else
    %load saved results & plot
    classifiedDataMat = 'classifiedData.mat';
    load([rootFolder classifiedDataMat]);
    classifiedDataKinMat = 'classifiedDataKin.mat';
    load([rootFolder classifiedDataKinMat]);

    % pass anatNames & plot all together:
    plotClassifierAccuracyOverTime(classifiedData, classifiedKinData, timescaleSegClassify, timescaleSegKinClassify, numSampRepeats)
end

% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function zeroXDottedLine
% for current axes, draw a dotted line at X=0 extending to yLimits
    x = [0 0]; %dotted line at event = time 0
    yLims = get(gca, 'YLim');
    y = [yLims(1) yLims(2)];
    plot(x, y, 'k--', 'LineWidth', 0.1);
end

function [accuracyOverTime] = getNonlinearTrialTimeClassifierLSTM(trialLabels, trialData)
% function to take trial labels (whether avoid or react) and train nonlinear classifier using neural or kinematics data
% train decoder at each time step, and return average cross-validated accuracy using 5 folds
% uses LSTM deep learning architecture w/ backprop gradient descent

    %input is trialLabels (binary, 1=avoid), and neuralData w/ associated timescaleSeg is time series leading up to extension
    numTimepoints = size(trialData{1},2);
    numTrials = size(trialData,2);
    numUnits = size(trialData{1},1);

    accuracyOverTime = zeros(numTimepoints,1);
    
    % loop over timepoints
    hWait = waitbar(1/numTimepoints,'Creating LSTM model for each timepoint...'); %set up progress bar

    for i = 1:numTimepoints
        waitbar(i/numTimepoints, hWait); %update progress bar
        for t = 1:numTrials
            X{t} = trialData{t}(:,i); %get single timepoint
            Y(t) = trialLabels{t};
        end
        
        % Define LSTM architecture
        numHiddenUnits = 100;
        layers = [
            sequenceInputLayer(numUnits)
            lstmLayer(numHiddenUnits,OutputMode="last")
            fullyConnectedLayer(2) %this is the number of classes to predict
            softmaxLayer % for binary fullyConnectedLayer this is logistic function; normalizes output so it can be interpreted as probability
            classificationLayer];
        
        % Define cross-validation options
        cv = cvpartition(numTrials,'KFold',5);
        
        % Initialize variables to store cross-validation results
        meanAccuracy = zeros(cv.NumTestSets,1);

        % Perform cross-validation
        for j = 1:cv.NumTestSets
            % Split data into training and testing sets
            trainIdx = training(cv,j);
            testIdx = test(cv,j);
            XTrain = X(:,trainIdx);
            YTrain = Y(trainIdx);
            XTest = X(:,testIdx);
            YTest = Y(testIdx);

            options = trainingOptions("adam", ...
                ExecutionEnvironment="gpu", ...
                GradientThreshold=1, ...
                MaxEpochs=20, ...
                SequenceLength="longest", ...
                Shuffle="never", ...
                Verbose=0, ...
                Plots="none"); %"training-progress"

            %train LSTM model
            net = trainNetwork(XTrain,YTrain,layers,options);
            
            % Predict on test set
            YPred = classify(net,XTest');
            
            % Calculate accuracy
            meanAccuracy(j) = sum(YPred' == YTest) / numel(YTest);
        end
        
        % Display cross-validation results
        mean_accuracy = mean(meanAccuracy);
%         fprintf('Mean Accuracy: %.2f%%\n', mean_accuracy * 100);   
        accuracyOverTime(i) = mean_accuracy;
    end %end loop over timepoints
    accuracyOverTime = accuracyOverTime'; %make column vector for consistent concatenation
    close(hWait);
end

function [accuracyOverTime] = getNonlinearTrialTimeClassifierNet(trialLabels, trialData)
% function to take trial labels (whether avoid or react) and train nonlinear classifier using neural or kinematics data
% train decoder at each time step, and return average cross-validated accuracy using 5 folds
% uses non-sequence/LSTM deep learning architecture w/ backprop gradient descent, with just a fully-connected network instead of LSTM

    %input is trialLabels (binary, 1=avoid), and neuralData w/ associated timescaleSeg is time series leading up to extension
    numTimepoints = size(trialData{1},2);
    numTrials = size(trialData,2);
    numUnits = size(trialData{1},1);

    accuracyOverTime = zeros(numTimepoints,1);
    
    % loop over timepoints
    hWait = waitbar(1/numTimepoints,'Creating neural network model for each timepoint...'); %set up progress bar

    for i = 1:numTimepoints
        waitbar(i/numTimepoints, hWait); %update progress bar
        for t = 1:numTrials
            X{t} = trialData{t}(:,i); %get single timepoint
        end
        Y = trialLabels;
        
        % Define network architecture
        numHiddenUnits = 100;
        layers = [
            featureInputLayer(numUnits) %use featureInput layer if no time or spatial dimensions (which is true for time-dependent classifier)
            fullyConnectedLayer(numHiddenUnits) %instead of sequence layer for time-dependent classifier, just use a fully connected layer
            fullyConnectedLayer(2) %this is the number of classes to predict
            softmaxLayer % for binary fullyConnectedLayer this is logistic function; normalizes output so it can be interpreted as probability
            classificationLayer];

        %NOTE: tried batchNormalizationLayer & reluLayer between fullConnectedLayers, but no appreciable speed up;  % "To speed up training of the neural network and reduce the sensitivity to network initialization, use batch normalization layers and nonlinearities, such as ReLU layers"
            
        % Define cross-validation options
        cv = cvpartition(numTrials,'KFold',5);
        
        % Initialize variables to store cross-validation results
        meanAccuracy = zeros(cv.NumTestSets,1);

        % Perform cross-validation
        for j = 1:cv.NumTestSets
            % Split data into training and testing sets
            trainIdx = training(cv,j);
            testIdx = test(cv,j);
            XTrain = X(:,trainIdx);
            YTrain = Y(trainIdx);
            XTest = X(:,testIdx);
            YTest = Y(testIdx);

            options = trainingOptions("adam", ...
                ExecutionEnvironment="gpu", ...
                GradientThreshold=1, ... %regularize gradient/weights by default L2-norm method (i.e L2norm = GradientThreshold)
                MaxEpochs=20, ... % max epochs used for training; an epoch is the full pass of the training algorithm over the entire training set.
                Shuffle="never", ...
                Verbose=0, ...
                Plots="none"); %"training-progress"

            %train LSTM model
            XTrain = cell2mat(XTrain)'; % non-sequence network requires matrix input instead of cell, with # rows = length(Y)
            net = trainNetwork(XTrain,YTrain,layers,options);
            
            % Predict on test set
            XTest = cell2mat(XTest)';
            YPred = classify(net,XTest);
            
            % Calculate accuracy
            meanAccuracy(j) = sum(YPred' == YTest) / numel(YTest);
        end
        
        % Display cross-validation results
        mean_accuracy = mean(meanAccuracy);
%         fprintf('Mean Accuracy: %.2f%%\n', mean_accuracy * 100);   
        accuracyOverTime(i) = mean_accuracy;
    end %end loop over timepoints
    accuracyOverTime = accuracyOverTime'; %make column vector for consistent concatenation
    close(hWait);
end

function [accuracyOverTime] = getNonlinearTrialTimeClassifierSVM(trialLabels, trialData)
% function to take trial labels (whether avoid or react) and train nonlinear classifier using neural or kinematics data
% train decoder at each time step, and return average cross-validated accuracy using 5 folds
% uses SVM

    %input is trialLabels (binary, 1=avoid), and neuralData w/ associated timescaleSeg is time series leading up to extension
    numTimepoints = size(trialData{1},2);
    numTrials = size(trialData,2);
    numUnits = size(trialData{1},1);

    accuracyOverTime = zeros(numTimepoints,1);
    
    % loop over timepoints
%     hWait = waitbar(1/numTimepoints,'Creating SVM model for each timepoint...'); %set up progress bar

    for i = 1:numTimepoints
%         waitbar(i/numTimepoints, hWait); %update progress bar
        for t = 1:numTrials
            X{t} = trialData{t}(:,i); %get single timepoint with 10 PCs from shared subspace
            Y(t) = trialLabels{t};
        end
        
        % Define cross-validation options
        cv = cvpartition(numTrials,'KFold',5);
        
        % Initialize variables to store cross-validation results
        meanAccuracy = zeros(cv.NumTestSets,1);

        % Perform cross-validation
        for j = 1:cv.NumTestSets
            % Split data into training and testing sets
            trainIdx = training(cv,j);
            testIdx = test(cv,j);
            XTrain = X(:,trainIdx);
            YTrain = Y(trainIdx);
            XTest = X(:,testIdx);
            YTest = Y(testIdx);

            %train SVM with nonlinear kernel
            XTrain = cell2mat(XTrain)'; % SVM requires matrix input instead of cell, with # rows = length(Y)
            %rng(1); %specify random number seed to reproducible results with KernelScale=auto below
            SVMModel = fitcsvm(XTrain, YTrain,'KernelFunction','gaussian','BoxConstraint',1,'KernelScale',1,'CacheSize','maximal','Verbose', 0,'NumPrint',1000); %set Verbose=1 to print convergence & debug
            %perhaps try 'OptimizeHyperparameters', 'auto' to set BoxConstraint & KernelScale; initial run gives 'BoxContraint values <1000
            % BoxConstraint is regularization parameter C to prevent overfitting, but trades off with computation time

            % Predict on test set
            XTest = cell2mat(XTest)';
            YPred = predict(SVMModel, XTest);

            % Calculate accuracy
            meanAccuracy(j) = sum(YPred' == YTest) / numel(YTest);
        end
        
        % Display cross-validation results
        mean_accuracy = mean(meanAccuracy);
%         fprintf('Mean Accuracy: %.2f%%\n', mean_accuracy * 100);   
        accuracyOverTime(i) = mean_accuracy;
    end %end loop over timepoints
    accuracyOverTime = accuracyOverTime'; %make column vector for consistent concatenation
%     close(hWait);
end

function plotClassifierAccuracyOverTime(classifiedData, classifiedKinData, timescaleSegClassify, timescaleSegClassifyKin, numSampRepeats)
 %combined ctx or non-ctx units avoid vs react mean responses
 

    probeToPlot = {'Probe1','Probe0','Probe0','Probe2'}; %matched pairs with below
    anatToPlot = {'CTX','CTX','BS','HB'};
%     labels = {'PFC','M1','TH/HY','HB','shuffle'}; 
    labels = {'kinematics','PFC','M1','TH/HY','HB','shuffle'}; 
    lineCmaps = linspecer(size(probeToPlot,2));
    anatCmaps =  [lineCmaps(2,:); lineCmaps(1,:); lineCmaps(3,:); lineCmaps(4,:)]; %reorder to put PFC first

    cmapShuf = [0.2 0.2 0.2];
    cmapKin = [0.4940 0.1840 0.5560]; %same purple as other kin plots
    smoothWindowEphys = [50 0]; %ex. 10 x 20ms bins = 200ms
    smoothWindowKin = [50 0]; % also matched 20ms 'bins'
    smoothOn = 1;
    f = figure; hold on;

    allAccuracyOverTimeShuf = []; %collect shuffles across all areas + kinematics

    % EPHYS data: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:size(probeToPlot,2) %main loop over probe structures
        for rpt = 1:numSampRepeats
            accuracyOverTime(rpt,:) = classifiedData(rpt).(probeToPlot{1,i}).(anatToPlot{1,i}).classifierAccuracyOverTime;
            accuracyOverTimeShuf(rpt,:) = classifiedData(rpt).(probeToPlot{1,i}).(anatToPlot{1,i}).classifierAccuracyOverTimeShuffle;
            if smoothOn
                accuracyOverTime(rpt,:) = smoothdata(accuracyOverTime(rpt,:),2,'gaussian',smoothWindowEphys);
                accuracyOverTimeShuf(rpt,:) = smoothdata(accuracyOverTimeShuf(rpt,:),2,'gaussian',smoothWindowEphys);
            end
%             allAccuracyOverTimeShuf(rpt,:) = [allAccuracyOverTimeShuf(rpt,:); accuracyOverTimeShuf(rpt,:)]; %collect shuffles across all areas
        end
        allAccuracyOverTimeShuf = [allAccuracyOverTimeShuf; mean(accuracyOverTimeShuf,1)]; %collect [mean over repeats] shuffles across all areas

        [meanAccuracy, stdAccuracy] = grpstats(accuracyOverTime,[],{'mean' 'std'}); %#ok<*ASGLU> %mean over repeats
        [h(i),~] = boundedline(timescaleSegClassify, meanAccuracy, stdAccuracy, 'cmap', anatCmaps(i,:), 'LineWidth', 2); 
%         h(i) = plot(timescaleSegClassify, meanAccuracy,'Color',anatCmaps(i,:),'LineWidth',3); %option for only mean if std makes lot too busy

    end %end probe loop

    % KINEMATICS data: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for rpt = 1:numSampRepeats
        kinAccuracyOverTime(rpt,:) = classifiedKinData(rpt).classifierAccuracyOverTime;
        kinAccuracyOverTimeShuf(rpt,:) = classifiedKinData(rpt).classifierAccuracyOverTimeShuffle;
        if smoothOn
            kinAccuracyOverTime(rpt,:) = smoothdata(kinAccuracyOverTime(rpt,:),2,'gaussian',smoothWindowKin);
            kinAccuracyOverTimeShuf(rpt,:) = smoothdata(kinAccuracyOverTimeShuf(rpt,:),2,'gaussian',smoothWindowKin);
        end
    end

    [meanKinAccuracy, stdKinAccuracy] = grpstats(kinAccuracyOverTime,[],{'mean' 'std'}); %#ok<*ASGLU> %mean over repeats
    [meanKinAccuracyShuf, stdKinAccuracyShuf] = grpstats(kinAccuracyOverTimeShuf,[],{'mean' 'std'}); %#ok<*ASGLU> %mean over repeats

    [meanKinAccuracy, stdKinAccuracy] = grpstats(kinAccuracyOverTime,[],{'mean' 'std'}); 
    %[meanKinAccuracyShuf, stdKinAccuracyShuf] = grpstats(kinAccuracyOverTimeShuf,[],{'mean' 'std'});
    [hk,~] = boundedline(timescaleSegClassifyKin, meanKinAccuracy, stdKinAccuracy, 'cmap', cmapKin, 'alpha','y', 'transparency', 0.5, 'LineWidth', 2); %prob too busy plot w/ SEMs as well
    %[hks,~] = boundedline(timescaleSegClassifyKin, meanKinAccuracyShuf, stdKinAccuracyShuf, 'cmap', cmapShuf, 'alpha','y', 'transparency', 0.5, 'LineWidth', 2);
    %hk = plot(timescaleSegKinDec, meanKinAccuracy,'Color',cmapKin,'LineWidth',3);
    %hks = plot(timescaleSegKinDec, meanKinAccuracyShuf,'Color',cmapShuf,'LineWidth',3);

    % SHUFFLE control: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    allAccuracyOverTimeShuf = [allAccuracyOverTimeShuf; meanKinAccuracyShuf]; %collect [mean over repeats] shuffles across all areas
    [meanAccuracyShuf, stdAccuracyShuf] = grpstats(accuracyOverTimeShuf,[],{'mean' 'std'}); %#ok<*ASGLU> %mean over areas + kinematics
    [hs,~] = boundedline(timescaleSegClassify, meanAccuracyShuf, stdAccuracyShuf, 'cmap', cmapShuf, 'LineWidth', 2); 

    xLims = get(gca, 'XLim');
    axis([-1 5 0.3 1]);
    yLims = get(gca, 'YLim');
    zeroXDottedLine; %dotted line 0 at plastformCommand / cue
    plot([3,3],yLims,'Color', 'b', 'LineStyle', '--', 'LineWidth', 1); %dotted line for latency=3 (cold air on - platform inertial delay)
    hold off;
    xlabel('time (sec)');
    ylabel('avoid/react prediction accuracy');
%     legend([h hs], labels, 'Location','northwest');
    legend([hk h hs], labels, 'Location','northwest');
    legend('boxoff')
    set(gcf,'color','w'); %set figure background white
    makepretty;

end





