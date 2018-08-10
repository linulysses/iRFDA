classdef PerformanceMeasure < handle
    properties
        % estimation results
%         muTrue
%         muHat
%         betaTrue
%         betaHat
%         alphaTrue
%         alphaHat
%         phiTrue
%         phiHat
        Kopt
        
        % measures
        muRMSE % rMSE for mean curve
        xRMSE % rMSE for FPCA
        yRMSE % rMSE for FLR
        phiRMSE
        betaRMSE % slope function
        betaRMSEopt
        yRMSEopt;
        lambda
        
        nSimu % # simulations
        nMethod % # of methods
        methodNames
        mfd  % manifold model
        K
    end
    
    methods
        
        % methodNames: cell array of names of methods
        function obj = PerformanceMeasure(nSimu,K,methodNames,mfd,...
                Kopt,muRMSE,xRMSE,yRMSE,phiRMSE,betaRMSE,...
                yRMSEopt,betaRMSEopt,lambda)
            obj.nSimu = nSimu;
            obj.K = K;
            obj.methodNames = methodNames;
            obj.nMethod = length(methodNames);
            obj.mfd = mfd;
            
            if nargin <= 4
                
                obj.Kopt = zeros(obj.nMethod,nSimu);
                obj.muRMSE = zeros(obj.nMethod,nSimu);
                obj.xRMSE = zeros(obj.nMethod,nSimu,K);
                obj.yRMSE = zeros(obj.nMethod,nSimu,K);
                obj.phiRMSE = zeros(2,K,obj.nSimu,obj.nMethod);
                obj.betaRMSE = zeros(obj.nMethod,nSimu,K);
                obj.betaRMSEopt = zeros(obj.nMethod,nSimu);
                obj.yRMSEopt = zeros(obj.nMethod,nSimu);
                obj.lambda = zeros(obj.nMethod,nSimu,K);
                
            else
                obj.Kopt = Kopt;
                obj.muRMSE = muRMSE;
                obj.xRMSE = xRMSE;
                obj.yRMSE = yRMSE;
                obj.phiRMSE = phiRMSE;
                obj.betaRMSE = betaRMSE;
                obj.betaRMSEopt = betaRMSEopt;
                obj.yRMSEopt = yRMSEopt;
                obj.lambda = lambda;

            end
        end
        
        function AddSingle(obj,r,dat,rs,Ks,myRMSEforX)
            if nargin <=4
                Ks = 1:obj.K;
            end
            
            j = 0;
            for jj = 1:obj.nMethod
                if strcmp(rs.Name,obj.methodNames{jj})
                    j = jj;
                    break;
                end
            end
            obj.lambda(j,r,:) = rs.lam;
            obj.muRMSE(j,r) = rMSEforMean(obj,dat.mu,rs.mu);
            if nargin <= 5 || isempty(myRMSEforX)
                obj.xRMSE(j,r,1:length(Ks)) = rMSEforX(obj,dat.test.n,dat.test.X,rs.XXNew,Ks);
            else
                obj.xRMSE(j,r,1:length(Ks)) = myRMSEforX(obj,dat.test.n,dat.test.X,rs.XXNew,Ks);
            end
            obj.yRMSE(j,r,:) = rMSEforY(obj,dat.test.y,rs.FLR.YNew);
            obj.betaRMSE(j,r,:) = rMSEforBeta(obj,dat.beta,rs.FLR.beta);
            obj.yRMSEopt(j,r) = obj.yRMSE(j,r,rs.FLR.optK);
            obj.betaRMSEopt(j,r) = obj.betaRMSE(j,r,rs.FLR.optK);
            obj.Kopt(j,r) = rs.FLR.optK;
            obj.phiRMSE(:,:,r,j) = rMSEforPhi(obj,dat.mu,rs.mu,dat.phi,rs.phiV);
            
        end
        
        function Add(obj,r,dat,rslts)
            for j = 1:length(rslts)
                rs = rslts{j};
                addSingle(obj,r,dat,rs);
            end
        end
        
        function ProduceStatistics(obj,fullinfo)
            if nargin == 1
                fullinfo = false;
            end
            mcMeanMuRMSE = mean(obj.muRMSE,2);
            mcSdMuRMSE = std(obj.muRMSE,1,2);
            
            fprintf('\n======== Result for mean curve estimation:\n');
            for j = 1:obj.nMethod
                fprintf('\t%s -- RMSE = %f,  std = %f \n', ...
                    obj.methodNames{j},mcMeanMuRMSE(j),mcSdMuRMSE(j));
            end
            
            mcMeanXRMSE = squeeze(mean(obj.xRMSE,2));
            mcSdXRMSE = squeeze(std(obj.xRMSE,1,2));
            fprintf('\n======== Result for curve recovery by FPCA:\n');
            for j = 1:obj.nMethod
                fprintf('\t%s',obj.methodNames{j});
            end
            fprintf('\n');
            fprintf('\tRMSE:\n');
                disp(mcMeanXRMSE');
            fprintf('\tStandard Error:\n');
                disp(mcSdXRMSE');
                
            if fullinfo    
                mcMeanBetaRMSE = squeeze(mean(obj.betaRMSE,2));
                mcSdBetaRMSE = squeeze(std(obj.betaRMSE,1,2));
                fprintf('\n======== Result for estimation of beta:\n');
                for j = 1:obj.nMethod
                    fprintf('\t%s',obj.methodNames{j});
                end
                fprintf('\n');
                fprintf('\tRMSE:\n');
                    disp(mcMeanBetaRMSE');
                fprintf('\tStandard Error:\n');
                    disp(mcSdBetaRMSE');


                mcMeanYRMSE = squeeze(mean(obj.yRMSE,2));
                mcSdYRMSE = squeeze(std(obj.yRMSE,1,2));
                fprintf('Result for prediction of Y:\n');
                for j = 1:obj.nMethod
                    fprintf('\t%s',obj.methodNames{j});
                end
                fprintf('\n');
                fprintf('\tRMSE:\n');
                    disp(mcMeanYRMSE');
                fprintf('\tStandard Error:\n');
                    disp(mcSdYRMSE');
            end

            mcMeanBetaRMSE = mean(obj.betaRMSEopt,2);
            mcSdBetaRMSE = std(obj.betaRMSEopt,1,2);
            fprintf('\n======== Result for estimation of beta (CV):\n');
            for j = 1:obj.nMethod
                fprintf('\t%s',obj.methodNames{j});
            end
            fprintf('\n');
            fprintf('\tRMSE:\n');
                disp(mcMeanBetaRMSE');
            fprintf('\tStandard Error:\n');
                disp(mcSdBetaRMSE');
               
            
            mcMeanYRMSE = mean(obj.yRMSEopt,2);
            mcSdYRMSE = std(obj.yRMSEopt,1,2);
            fprintf('\n======== Result for prediction of Y (CV):\n');
            for j = 1:obj.nMethod
                fprintf('\t%s',obj.methodNames{j});
            end
            fprintf('\n');
            fprintf('\tRMSE:\n');
                disp(mcMeanYRMSE');
            fprintf('\tStandard Error:\n');
                disp(mcSdYRMSE');
           
            mcMeanPhiRMSE = squeeze(mean(obj.phiRMSE,3));
            mcSdPhiRMSE = squeeze(std(obj.phiRMSE,1,3));
            
            fprintf('\n======== Result for eigenfunction (method by method):\n');
            for j = 1:obj.nMethod
                fprintf('\t%s: \n',obj.methodNames{j});
            
                fprintf('\n');
                fprintf('\tRMSE:\n');
                    disp(mcMeanPhiRMSE(:,:,j));
                fprintf('\tStandard Error:\n');
                    disp(mcSdPhiRMSE(:,:,j));
            end
            
            fprintf('\n======== Result for eigenfunction (intrinsic only):\n');
            for j = 1:obj.nMethod
                fprintf('\t%s',obj.methodNames{j});
            end
                fprintf('\n');
                fprintf('\tRMSE:\n');
                    disp(squeeze(mcMeanPhiRMSE(1,:,:)));
                fprintf('\tStandard Error:\n');
                    disp(squeeze(mcSdPhiRMSE(1,:,:)));
           
        end
        
        % XX is a cell array
        function [XrMSE] = rMSEforX(obj,n,X,XX,Ks)
            if nargin == 4
                Ks = 1:obj.K;
            end
            XrMSE = zeros(1,length(Ks));
            for j = 1:length(Ks)
                k = Ks(j);
                for i = 1:n;
                    D = obj.mfd.dist(X(:,:,i),XX{k}(:,:,i));
                    XrMSE(j) = XrMSE(j) + mean(D.^2);
                end
                XrMSE(j) = sqrt(XrMSE(j) / n);
            end
        end
        
        function [PrMSE] = rMSEforPhi(obj,muTrue,muHat,phiTrue,phiHat)
            iRMSE = zeros(1,obj.K);
            ExRMSE = zeros(1,obj.K);
            l2vfinprod = @(U,V) mean(sum((U.*V),1));
            l2vfnorm = @(U) sqrt(l2vfinprod(U,U));
            maxK = size(phiTrue,3);
            maxK = min(maxK,obj.K);
            for k = 1:maxK
                phik1 = phiHat(:,:,k);
                phik0 = phiTrue(:,:,k);
                [phik1_p] = obj.mfd.parallel_transport(muHat,muTrue,phik1);
                if obj.mfd.vfinprod_V(phik1_p,phik0,muTrue) < 0
                    phik1_p = - phik1_p;
                end
                iRMSE(k) = sqrt(obj.mfd.vfinprod_V(phik1_p-phik0,phik1_p-phik0,muTrue));

                if l2vfinprod(phik0,phik1) < 0
                    phik1 = -phik1;
                end
                ExRMSE(k) = l2vfnorm(phik1-phik0);
            end
            PrMSE = [iRMSE; ExRMSE];
        end
        
        % XX is a cell array
        function [BrMSE] = rMSEforBeta(obj,betaTrue,betaHat,Ks)
            if nargin == 3
                Ks = 1:obj.K;
            end
            BrMSE = zeros(1,length(Ks));
            for j = 1:length(Ks)
                k = Ks(j);
                BrMSE(j) = mean(obj.mfd.dist(betaTrue,betaHat(:,:,k)).^2);
                %BrMSE(j) = mean(sum((betaTrue-betaHat(:,:,k)).^2,1));
            end
        end
        
        % XX is a cell array
        function [YrMSE] = rMSEforY(obj,yTrue,yPred)
            YrMSE = sqrt(mean((repmat(yTrue,1,obj.K)-yPred).^2,1));
        end
        
        function [rmse] = rMSEforMean(obj,muTrue,muEst)
            rmse = sqrt(mean(obj.mfd.dist(muTrue,muEst).^2));
        end
    end

end