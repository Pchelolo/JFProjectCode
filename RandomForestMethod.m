classdef RandomForestMethod < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        %Settings. Change before using!!!
        corr_period = 20;   %Period in samples for the local correlation metric. Should be >20, approximetely equal to the time resolution of effect you are looking for
        firstChop = 10000;  %The edge of the training part of the data. in samples.
        secondChop = 20000; %The edge of the cross-validation part of the data, in samples
        cutoff = 1.0;       %The cutoff distance for event clustering
        samplingRate = 20;  %in Hz
        maxlag = 2;         %number of lags to include backwards and forwards.
        smoothPeriod = 5;   %Smooth firind rate with Gaussian, std = smoothPeriod in samples.
        
        paroptions = statset('UseParallel','always'); %Parallel options!
    end
    
    properties
        r;  %firing rate of all neurons. Each row is a neuron
        X;  %Behavioral variables. Each Column!! is a variable.
        Y;  %Smoothed neural firing rate, each Column is a variable
        
        thresholdMean = {}; %each cell has a vector for event significancy threshold
        thresholdSTD = {};  %
        trees = {};         %Stored forests we alreadt calculated
    end
    
    methods
        function this = RandomForestMethod(dataSource, varargin)
           tic;
            %   dataSource - sets the type of data you are loading.
            %   	'extinctionDataset' - need toprovide NS3 file and PopData
            %       'learningDataset' - provide ....
            %       'other' - provide   Y - firing rate for all neurons, each
            %                           row is a neuron
            %                           X - behavioral variables, each row
            %                           is a separate variable
            if matlabpool('size') == 0
                matlabpool local;
            end
            addpath('utils/');
            
            if strcmp(dataSource, 'extDataset')
                if nargin ~= 3
                    error('Wrong number of arguments');
                end 
                
                fprintf('Loading data for extinction dataset');
                
                data = varargin{1}.session;
                NS3 = varargin{2};
                
                minTime = 0;
                maxTime = ceil(size(NS3.Data, 2) / 2);
                diffTime = 1000 / RandomForestMethod.samplingRate;
                fprintf('Using time window %d milliseconds\n', diffTime);
                times = minTime : diffTime : maxTime + diffTime; %create the vector of timebins
                this.r = zeros(length(data.unit), length(times)); %Allocate the memory for the firing rate
                for cell = 1 : length(data.unit)
                    cell_timestamps = data.unit(cell).ts;
                    k = 1;
                    for idx = 1 : length(cell_timestamps)
                        while times(k) < cell_timestamps(idx)
                            k = k + 1;
                        end
                        this.r(cell, k) = this.r(cell, k) + 1;
                    end
                end
                this.r = this.r(:, 1:end-1)*20;
                fprintf('Max number of spikes in a bin %d\n', max(max(this.r)));
                fprintf('Size of the data matrix %d %d\n', size(this.r)); 
    
                this.Y = MakeSmoothWithGaussian(this.r(:,RandomForestMethod.maxlag+1:end-RandomForestMethod.maxlag), 5)';
                
                this.X = [];
                for lag = -RandomForestMethod.maxlag*100:100:RandomForestMethod.maxlag*100
                    this.X = [this.X, NS3.Data(35:40, (RandomForestMethod.maxlag*100+lag+1):100:(end-RandomForestMethod.maxlag*100+lag))'];
                end
            elseif strcmp(dataSource, 'learnDataset')
                if nargin ~= 4
                    error('Wrong number of arguments');
                end
                
                NEV = varargin{1};
                NS3 = varargin{2};
                NEX = varargin{3};
                
                minTime = 0;
                maxTime = ceil(size(NS3.Data, 2) / 2);
                diffTime = 1000 / RandomForestMethod.samplingRate;
                fprintf('Using time window %d milliseconds\n', diffTime);
                times = minTime : diffTime : maxTime + diffTime; %create the vector of timebins in milliseconds!!!
                r_temp = zeros(length(unique(NEV.Data.Spikes.Electrode)) - 2, length(times)-1); %Allocate the memory for the firing rate
                
                hist_edges = zeros(1, length(times)-1);
                timestamps = NEV.Data.Spikes.Timestamps;
                parfor i = 1:length(times)-1
                    hist_edges(i) = find(timestamps > times(i)*30, 1, 'first');
                end
                uniqElec = unique(NEV.Data.Spikes.Electrode);
                electrodes = NEV.Data.Spikes.Electrode;
                parfor n = 1 : length(uniqElec) - 2  
                    r_temp(n, :) = histc(find(electrodes == uniqElec(n)), hist_edges);
                end
                this.r = r_temp;
                fprintf('Max number of spikes in a bin %d\n', max(max(this.r)));
                fprintf('Size of the data matrix %d %d\n', size(this.r)); 
    
                this.Y = MakeSmoothWithGaussian(this.r(:,RandomForestMethod.maxlag+1:end-RandomForestMethod.maxlag), RandomForestMethod.smoothPeriod)';
                
                events = zeros(3, length(times)-1);
                events(1, ceil(NEX.events{1}.timestamps*20)) = 1;
                events(2, ceil(NEX.events{2}.timestamps*20)) = 1;
                events(3, ceil(NEX.events{3}.timestamps*20)) = 1;
                events = MakeSmoothWithGaussian(events, RandomForestMethod.smoothPeriod);
                
                this.X = [];
                for lag = -RandomForestMethod.maxlag:RandomForestMethod.maxlag
                    this.X = [this.X, NS3.Data(35:40, (RandomForestMethod.maxlag*100+lag*100+1):100:(end-RandomForestMethod.maxlag*100+lag*100))'];
                    this.X = [this.X, events(:, (RandomForestMethod.maxlag+lag+1):(end-RandomForestMethod.maxlag+lag))'];
                end
                
            elseif strcmp(dataSource, 'other')
                this.r = varargin{1};
                this.X = varargin{2};
                this.Y = varargin{3};
            else
                error('Unknown datasource');  
            end
            toc;
        end
        
        function result = fullAnalysis(this, neuron, varargin)
            result = struct('corr', 0, 'tree', 0, 'timestamps', 0, 'ts_clust', 0, 'corrSignificancy', 0);
            
            recalculate = 0;
            if nargin > 2
                recalculate = strcmp(varargin{end}, 'recalculate');    
            end
            if recalculate == 1 || length(this.thresholdMean) < neuron || isempty(this.thresholdMean{neuron}) 
                this.doBootstrap(neuron)
            end
            
            if recalculate == 1 || length(this.trees) < neuron || isempty(this.trees{neuron})
                res = TreeBagger(100, this.X(1:RandomForestMethod.firstChop, :), this.Y(1:RandomForestMethod.firstChop), 'method','r','oobpred','on', 'oobvarimp', 'on', 'minleaf',5, 'Options', RandomForestMethod.paroptions);
                this.trees{neuron} = res;
            else
                res = this.trees{neuron};
            end  
            
            Y1 = res.predict(this.X(RandomForestMethod.firstChop:RandomForestMethod.secondChop, :));
            Y2 = this.Y(RandomForestMethod.firstChop:RandomForestMethod.secondChop, neuron);    
            fprintf('Sum square error %f\n', sum(sqrt(sum((Y1-Y2).^2))));
             
            corrCoef = arrayfun(@(t) corr(Y2(t-RandomForestMethod.corr_period:t+RandomForestMethod.corr_period), Y1(t-RandomForestMethod.corr_period:t+RandomForestMethod.corr_period)), RandomForestMethod.corr_period+1:length(Y2)-RandomForestMethod.corr_period);
            result.corr = corrCoef;
            [ts, c] = peakfinder(corrCoef);
            threhold = this.thresholdMean{neuron} + 1.96*this.thresholdSTD{neuron};
            filter = [];
            for i = 1 : length(ts)
                if c(i) > threhold(ts(i))
                    filter = [filter, i];
                end
            end
            ts = ts(filter) + RandomForestMethod.firstChop+RandomForestMethod.corr_period;
            c = c(filter);
            result.timestamps = ts;
    
            cl = classifyNeuralEvents(ts, RandomForestMethod.corr_period, MakeSmoothWithGaussian(this.r(neuron, :), RandomForestMethod.smoothPeriod), RandomForestMethod.cutoff);
            result.ts_clust = cl.C;
    
            
            %Plotting the summary
            figure;
            toseconds = @(x) x/RandomForestMethod.samplingRate;
            color = colormap(hsv(max(cl.C)));
            subplot(3,3,1:3);
            hold on;
            plot(toseconds(RandomForestMethod.firstChop+RandomForestMethod.corr_period+1:RandomForestMethod.secondChop-RandomForestMethod.corr_period), this.Y(RandomForestMethod.firstChop+RandomForestMethod.corr_period+1:RandomForestMethod.secondChop-RandomForestMethod.corr_period), 'k',...
                toseconds((RandomForestMethod.corr_period+1:length(Y1)-RandomForestMethod.corr_period)+RandomForestMethod.firstChop), Y1(RandomForestMethod.corr_period+1:end-RandomForestMethod.corr_period),  'r');
            for idx = 1 : length(ts)
                line(toseconds([ts(idx) ts(idx)]), [min(Y2) max(Y2)], 'Color', color(cl.C(idx),:));
            end
            hold off;
            xlim([min(toseconds((1:length(Y1))+RandomForestMethod.firstChop)), max(toseconds((1:length(Y1))+RandomForestMethod.firstChop))]);
            title('Cross validation results');
            xlabel('Time, seconds');
            ylabel('Firing rate, Hz');
            legend('Real firing rate', 'Model prediction');
    
            subplot(3,3,4:6);
            plot(toseconds((RandomForestMethod.corr_period+1:length(this.Y(RandomForestMethod.firstChop:RandomForestMethod.secondChop))-RandomForestMethod.corr_period)+RandomForestMethod.firstChop), corrCoef, 'k', toseconds(ts), c, 'or', toseconds((RandomForestMethod.corr_period+1:length(this.Y(RandomForestMethod.firstChop:RandomForestMethod.secondChop))-RandomForestMethod.corr_period)+RandomForestMethod.firstChop), threhold, 'g');
            xlim([min(toseconds((1:length(Y1))+RandomForestMethod.firstChop)), max(toseconds((1:length(Y1))+RandomForestMethod.firstChop))]);
            title('Fit estimation');
            xlabel('Time, seconds');
            ylabel('Correlation');
    
            subplot(3,3,7);
            clust_size = arrayfun(@(ind) length(find(cl.C == ind)), 1:max(cl.C));
            h = bar(1:max(cl.C), clust_size);
            h_ch = get(h, 'Children');
            set(h_ch, 'CData', 1:max(cl.C));
            xlabel('Cluster number');
            ylabel('Cluster size');
            title('Set of an event cluster');
    
            best_cluster = find(clust_size == max(clust_size),1);
            lag = -RandomForestMethod.corr_period:RandomForestMethod.corr_period;
            fprintf('Best cluster: %d\n', best_cluster);
            sta = STA(ts(cl.C == best_cluster), MakeSmoothWithGaussian(this.r(neuron,:), 5), lag, 0, 0, 0);
            subplot(3,3,8);
            plot(toseconds(lag), sta.av, 'r', toseconds(lag), sta.up, 'b', toseconds(lag), sta.down, 'b');
            title('Neural ETA');
            xlabel('Time from an event, seconds');
            ylabel('Average firing rate, Hz');
            
            sta = STA(ts(cl.C == best_cluster), this.X(:,4)', lag, 0, 0, 0);
            subplot(3,3,9);
            plot(toseconds(lag), sta.av, 'r', toseconds(lag), sta.up, 'b', toseconds(lag), sta.down, 'b');
            title('Behavioral ETA');
            xlabel('Time from an event, seconds');
            ylabel('Average behavior parameter');
        end
        
        function analyzePredictivePower(this, neuron, parameterSets)
           %parameterSets: matrix of 0s and 1s. each column correspond to one parameter, each row - one parameter set to exemine.
           %0 means the parameter is excluded(wrong alignment), 1 means it is inclusded. 
           
           %returns the mean and std of sqruare error for each parameter set for bootstrap
           %with different subsets of the data. 
           Y_temp = this.Y(:, neuron);
           %powersMean = zeros(1, size(parameterSets, 1));
           %powersSTD = zeros(1, size(parameterSets, 1));
           powersF = {};
           for i = 1 : size(parameterSets, 1)
               powers = [];
               X_temp = zeros(RandomForestMethod.firstChop+1, size(this.X, 2));
               fprintf('Analyzing parameter set number %d', i);
               for chop = 1000 : 1000 : 10000
                   disp(chop)
                   parametersToMove = repmat(parameterSets(i, :), 1, 2*RandomForestMethod.maxlag + 1);
                   x1 = this.X(chop:chop+RandomForestMethod.firstChop, :);
                   x2 = this.X(chop+2000 : chop+2000+RandomForestMethod.firstChop, :);
                   parfor j = 1 : length(parametersToMove)
                      if parametersToMove(j) == 1 
                            X_temp(:, j) = x1(:, j);
                      else
                            X_temp(:, j) = x2(:, j);
                      end
                   end
                   res = TreeBagger(100, X_temp, Y_temp(chop:chop+RandomForestMethod.firstChop), 'method','r','oobpred','on', 'oobvarimp', 'on', 'minleaf',5, 'Options', RandomForestMethod.paroptions);
                   Y1 = res.predict(this.X(RandomForestMethod.firstChop+chop:RandomForestMethod.secondChop+chop, :));
                   Y2 = this.Y(RandomForestMethod.firstChop+chop:RandomForestMethod.secondChop+chop, neuron);    
                   powers = [powers,  sum(sqrt(sum((Y1-Y2).^2)))];
               end
               powersF{i} = powers;
           end 
           save(sprintf('res%d', neuron), 'powersF');
        end
        
    end
    
    methods(Access = private)
        function doBootstrap(this, neuron)
            %Bootstraping over wrong aligned data to obtain a threshold
            %for statistical significant fits.
            fprintf('Making bootstrap. It might take a while');
            boost_val_trees = {};
            for shift = 2000 : 2000 : 4000
                boost_val_trees{end+1} = TreeBagger(100, this.X(shift+1:shift+RandomForestMethod.firstChop, :), this.Y(1:RandomForestMethod.firstChop,neuron), 'method','r','oobpred','on', 'oobvarimp', 'on', 'minleaf',5, 'Options', RandomForestMethod.paroptions);
            end
            corrCoef = [];
            for idx = 1 : length(boost_val_trees)
                Y1 = boost_val_trees{idx}.predict(this.X(RandomForestMethod.firstChop:RandomForestMethod.secondChop, :));
                Y2 = this.Y(RandomForestMethod.firstChop:RandomForestMethod.secondChop, neuron);    
                corrCoef = [corrCoef; arrayfun(@(t) corr(Y2(t-RandomForestMethod.corr_period:t+RandomForestMethod.corr_period), Y1(t-RandomForestMethod.corr_period:t+RandomForestMethod.corr_period)), RandomForestMethod.corr_period+1:length(Y2)-RandomForestMethod.corr_period)];
            end
            this.thresholdMean{neuron} = mean(corrCoef, 1);
            this.thresholdSTD{neuron} = std(corrCoef, [], 1)/sqrt(length(boost_val_trees));
            
        end
    end
end

    

