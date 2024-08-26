classdef MOE_LOG_GP_bounds <handle
%class for iterative learning with mixture of experts LoG-GPs
% E: state space dimension, N number of training samples
% Copyright (c) by Armin Lederer (TUM) under BSD License 
% Last modified: Markus Kessler 09/2022
% >> modified predict function: added output of mean gradient

    properties
        loadHyp = true;
    end
    
    properties
        pts = 100; %pts per local GP
        N = 10; %power of 2. Max amount of local GPs
        xSize = 6; %size of the training sample
        fileName = 'name'; %filename with hyperparameters. Named 'sigmaL','sigmaF','sigmaN'
        sigF = 1;
        sigL = 1;
        sigN = 0.01;
        divMethod = 1; %1: median 2: mean 3: mean(max,min)
        wo = 10; %ratio between width and overlapping region
        oEffect = [0 0 0]; %count data in [left, overlapping, right]
        divCount = 0;

        tau = 1e-5 % Grid constant
        delta = 0.2 % probability to violate bound
        Lf = 7 % Lipschitz constant of true f
        Lmean = 7;
    end
    
    properties(Access = protected)
        %properties pts,X,Y,K,invK and alpha must be set to the desired
        %amount of training points
        count; %amount of local GPs
        localCount; %points learned in each GP
        sigmaF;
        sigmaN;
        sigmaL;
        X; %training samples
        Y; %training targets
        K; %covariance matrices
        L; %cholesky factors
        alpha; %L'\(L\y)
        auxAlpha; % L\y
        
        medians; %vector of hyperplanes
        parent;
        children; %line 1: left child, line 2: right child
        overlaps; %line 1:cutting dimension, line 2:size of overlapping region
        
        auxUbic; %map a GP with the position of its data (K,L,alpha,X,Y,auxAlpha)

        Lk;         % Lipschitz constant of the kernel
        Lsig;      % Lipschitz constant for local posterior variance
        Lmu;        % Lipschitz constant for local posterior mean
        Xte;        % Grid to compute Lipschitz constants numerically
        XteMax; 
        XteMin;
        
    end
    
    methods
        function setupData(obj,in_Xte)
            %initialize data
            obj.count = 1;
            
            obj.X = zeros(obj.xSize, obj.pts * 2^obj.N);
            obj.Y = zeros(1, obj.pts * 2^obj.N);
            obj.K = zeros(obj.pts, obj.pts * 2^obj.N);
            obj.alpha = zeros(obj.pts,2^obj.N);
            obj.auxAlpha = zeros(obj.pts,2^obj.N);
            obj.L = zeros(obj.pts, obj.pts * 2^obj.N);
            obj.localCount = zeros(1,2* 2^obj.N -1);
            
            obj.medians =  zeros(obj.xSize, 2*2^obj.N-1);
            
            obj.parent = zeros(1, 2 * 2^obj.N-1);
            obj.children = -1*ones(2, 2 * 2^obj.N-1);
            
            obj.overlaps =  zeros(2, 2 * 2^obj.N-1);
            
            obj.auxUbic = zeros(1, 2 * 2^obj.N-1);
            obj.auxUbic(1,1) = 1;
            
            if obj.loadHyp %load or set hyperparameters
                aux = load(obj.fileName,'sigmaL','sigmaN','sigmaF');
                obj.sigmaL = aux.sigmaL;
                obj.sigmaF = aux.sigmaF;
                obj.sigmaN = aux.sigmaN;
            else
                obj.sigmaL = obj.sigL;
                obj.sigmaF = obj.sigF;
                obj.sigmaN = obj.sigN;
            end

            obj.Lk = norm(obj.sigmaF^2*exp(-0.5)./obj.sigmaL);
            obj.Lsig = zeros(1,obj.pts * 2^obj.N);
            obj.Lmu = obj.Lmean*ones(1,obj.pts * 2^obj.N);
            obj.Xte = in_Xte;
            obj.XteMax = max(in_Xte,[],2);
            obj.XteMin = min(in_Xte,[],2);
        end
        
        function kern = kernel(obj, Xi, Xj)
            kern = (obj.sigmaF^2)*exp(-0.5*sum(((Xi-Xj).^2)./(obj.sigmaL.^2),1))';
        end
        
        function m = mValue(obj, model,cutD)%compute the hyperplane
            if obj.divMethod == 1
                m = median(obj.X(cutD, (obj.auxUbic(model)-1)*obj.pts+1:...
                    obj.auxUbic(model)*obj.pts));
                return
            elseif obj.divMethod == 2
                m = mean(obj.X(cutD, (obj.auxUbic(model)-1)*obj.pts+1:...
                    obj.auxUbic(model)*obj.pts));
                return
            elseif obj.divMethod == 3
                m = (max(obj.X(cutD, (obj.auxUbic(model)-1)*obj.pts+1:...
                    obj.auxUbic(model)*obj.pts))+ min(obj.X(cutD, ...
                    (obj.auxUbic(model)-1)*obj.pts+1:...
                    obj.auxUbic(model)*obj.pts)))/2 ;
                return
            end
        end
        
        function updateParam(obj,x,model)
            pos = obj.auxUbic(model)-1;
            if obj.localCount(model) == 1 %first point in model
                lH = chol(obj.kernel(x, x) + obj.sigmaN.^2);
                obj.K(1,(pos)*obj.pts+1) = obj.kernel(x, x) + obj.sigmaN.^2;
                obj.L(1,(pos)*obj.pts+1) = lH;
                obj.alpha(1,pos+1) = lH'\(lH\obj.Y((pos)*obj.pts+1));
                obj.auxAlpha(1,pos+1) = (lH\obj.Y((pos)*obj.pts+1));
                obj.Lsig(:) = obj.sigmaF*norm(1./obj.sigmaL);
            else
                %set the updated parameters
                %auxX does not consider the new point x
                %auxY does consider the new point y
                auxX =  obj.X(:,(pos)*obj.pts+1:(pos)*obj.pts+obj.localCount(model)-1);
                auxY =  obj.Y((pos)*obj.pts+1:(pos)*obj.pts+obj.localCount(model));
                b = obj.kernel(auxX,x);
                c = obj.kernel(x,x)+obj.sigmaN^2;
                auxL = obj.L(1:obj.localCount(model)-1,...
                    (pos)*obj.pts+1:(pos)*obj.pts+obj.localCount(model)-1);
                newL = [auxL zeros(obj.localCount(model)-1,1); (auxL\b)' sqrt(c - norm(auxL\b)^2)];
                obj.K(1:obj.localCount(model), (pos)*obj.pts+1:...
                    (pos)*obj.pts+obj.localCount(model)) = [obj.K(1:obj.localCount(model)-1,...
                    (pos)*obj.pts+1:(pos)*obj.pts+obj.localCount(model)-1),...
                    b;b',c];
                obj.L(1:obj.localCount(model), (pos)*obj.pts+1:...
                    (pos)*obj.pts+obj.localCount(model)) = newL;
                
                obj.auxAlpha(obj.localCount(model),pos+1) = (auxY(obj.localCount(model))-...
                    newL(end,1:end-1)*obj.auxAlpha(1:obj.localCount(model)-1,pos+1))/...
                    newL(end,end);
                obj.alpha(1:obj.localCount(model),pos+1) = newL'\(obj.auxAlpha(1:obj.localCount(model),pos+1));
                
                obj.Lsig(:) = obj.sigmaF*norm(1./obj.sigmaL);
            end
            
        end
        
        function addPoint(obj, x, y, model)
            
            if obj.localCount(model) < obj.pts %if the model is not full
                obj.X(:,(obj.auxUbic(model)-1)*obj.pts+1+obj.localCount(model)) = x;
                obj.Y((obj.auxUbic(model)-1)*obj.pts+1+obj.localCount(model)) = y;
                obj.localCount(model) = obj.localCount(model) + 1;
                obj.updateParam(x,model)
            end
            if obj.localCount(model) == obj.pts %if full
                obj.divCount = obj.divCount + 1;
                div = 1;
                while div == 1 %divide until no child set has all the data
                    [div,model] = obj.divide(model);
                end
            end
        end
        
        function [div, childModel] =  divide(obj, model)
            if obj.parent(end)~= 0
                div = -1;
                childModel = -1;
                %                 disp('no room for new models')
                return;
            end
            %obtain cutting dimension
            [~,cutD]=max((max(obj.X(:, (obj.auxUbic(model)-1)*obj.pts+1:...
                obj.auxUbic(model)*obj.pts),[],2)-min(obj.X(:, (obj.auxUbic(model)-1)*obj.pts+1:...
                obj.auxUbic(model)*obj.pts),[],2)));
            %obtain hyperplane
            mP = obj.mValue(model,cutD);
            %compute borders
            maxV = max(obj.X(cutD, (obj.auxUbic(model)-1)*obj.pts+1:...
                (obj.auxUbic(model)-1)*obj.pts+obj.pts));
            minV = min(obj.X(cutD, ...
                (obj.auxUbic(model)-1)*obj.pts+1:...
                (obj.auxUbic(model)-1)*obj.pts+obj.pts));
            %compute overlapping region
            o  = (maxV-minV)/(obj.wo*(ceil(log(model))/10+1));
            
            obj.medians(model)=mP;
            obj.overlaps(1,model)=cutD;
            obj.overlaps(2,model)=o;
            
            xL = zeros(obj.xSize,obj.pts); %matrix with x values for the left model
            xR = zeros(obj.xSize,obj.pts); %matrix with x values for the right model
            yL = zeros(1,obj.pts); %vector with y values for the left model
            yR = zeros(1,obj.pts); %vector with y values for the left model
            
            lcount = 0;
            rcount = 0;
            iL = zeros(1,obj.pts); %left index order vector
            iR = zeros(1,obj.pts); %right index order vector
            
            for i=1:obj.pts
                xD = obj.X(cutD,(obj.auxUbic(model)-1)*obj.pts+i);%x value in cut dimension
                if xD<mP-o/2 %if in left set
                    lcount = lcount+1;
                    xL(:,lcount) = obj.X(:,(obj.auxUbic(model)-1)*obj.pts+i);
                    yL(lcount) = obj.Y((obj.auxUbic(model)-1)*obj.pts+i);
                    iL(lcount) = i;
                elseif xD >= mP-o/2 && xD <= mP+o/2 %if in overlapping
                    pL = 0.5 + (xD-mP)/(o);
                    if pL>=rand() %left side
                        lcount = lcount+1;
                        xL(:,lcount) = obj.X(:,(obj.auxUbic(model)-1)*obj.pts+i);
                        yL(lcount) = obj.Y((obj.auxUbic(model)-1)*obj.pts+i);
                        iL(lcount) = i;
                    else
                        rcount = rcount + 1;
                        xR(:,rcount) = obj.X(:,(obj.auxUbic(model)-1)*obj.pts+i);
                        yR(rcount) = obj.Y((obj.auxUbic(model)-1)*obj.pts+i);
                        iR(rcount) = i;
                    end
                elseif xD>mP+o/2 %if in right
                    rcount = rcount + 1;
                    xR(:,rcount) = obj.X(:,(obj.auxUbic(model)-1)*obj.pts+i);
                    yR(rcount) = obj.Y((obj.auxUbic(model)-1)*obj.pts+i);
                    iR(rcount) = i;
                end
            end
            
            obj.localCount(model) = 0; %divided set is now "empty"
            if obj.count == 1
                obj.count = obj.count+1; %update the total number of sets
            else
                obj.count = obj.count+2;
            end
            obj.children(:,model) = [obj.count obj.count+1]'; %assign the children
            obj.parent(obj.count:obj.count+1) = model; %assign the parent
            
            %update parameters of new models
            obj.localCount(obj.count) = lcount;
            obj.auxUbic(obj.count) = obj.auxUbic(model);
            
            obj.localCount(obj.count+1) = rcount;
            obj.auxUbic(obj.count+1) = max(obj.auxUbic)+1;
            
            if lcount == obj.pts  %if left set has all the points now
                div = 1; % output to keep dividing
                childModel = (obj.count);
            elseif rcount == obj.pts %if right set has all the points now
                %move L to the position of right model:
                obj.L(1:end, [(obj.auxUbic(model)-1)*obj.pts+1:...
                    (obj.auxUbic(model)-1)*obj.pts+obj.pts, ...
                    (obj.auxUbic(obj.count+1)-1)*obj.pts+1:...
                    (obj.auxUbic(obj.count+1)-1)*obj.pts+obj.pts]) =...
                    obj.L(1:end, [(obj.auxUbic(obj.count+1)-1)*obj.pts+1:...
                    (obj.auxUbic(obj.count+1)-1)*obj.pts+obj.pts, ...
                    (obj.auxUbic(model)-1)*obj.pts+1:...
                    (obj.auxUbic(model)-1)*obj.pts+obj.pts]);
                %move K to the position of right model:
                obj.K(1:end, [(obj.auxUbic(model)-1)*obj.pts+1:...
                    (obj.auxUbic(model)-1)*obj.pts+obj.pts, ...
                    (obj.auxUbic(obj.count+1)-1)*obj.pts+1:...
                    (obj.auxUbic(obj.count+1)-1)*obj.pts+obj.pts]) =...
                    obj.K(1:end, [(obj.auxUbic(obj.count+1)-1)*obj.pts+1:...
                    (obj.auxUbic(obj.count+1)-1)*obj.pts+obj.pts, ...
                    (obj.auxUbic(model)-1)*obj.pts+1:...
                    (obj.auxUbic(model)-1)*obj.pts+obj.pts]);
                %move alpha to the position of the right model
                obj.alpha(:,[obj.auxUbic(model),obj.auxUbic(obj.count+1)]) = ...
                    obj.alpha(:,[obj.auxUbic(obj.count+1),obj.auxUbic(model)]);
                obj.auxAlpha(:,[obj.auxUbic(model),obj.auxUbic(obj.count+1)]) = ...
                    obj.auxAlpha(:,[obj.auxUbic(obj.count+1),obj.auxUbic(model)]);
                div = 1; % output to keep dividing
                childModel = (obj.count+1);
            else %update alpha and  L values for the new models
                B = (1:obj.pts);
                C = [iL(1:lcount),iR(1:rcount)];
                newK = obj.K(1:end, (obj.auxUbic(model)-1)*obj.pts+1:...
                    (obj.auxUbic(model)-1)*obj.pts+obj.pts);
                %permute K:
                newK(B,:) =  newK(C,:);
                newK(:,B) = newK(:,C);
                
                %comoute child L factors
                lL = chol(newK(1:lcount,1:lcount),'lower');
                rL = chol(newK(lcount+1:end,lcount+1:end),'lower');
                
                obj.K(1:lcount,(obj.auxUbic(obj.count)-1)*obj.pts+1:...
                    (obj.auxUbic(obj.count)-1)*obj.pts+lcount) = newK(1:lcount,1:lcount);
                obj.K(1:rcount,(obj.auxUbic(obj.count+1)-1)*obj.pts+1:...
                    (obj.auxUbic(obj.count+1)-1)*obj.pts+rcount) = newK(lcount+1:end,lcount+1:end);
                
                obj.L(1:lcount,(obj.auxUbic(obj.count)-1)*obj.pts+1:...
                    (obj.auxUbic(obj.count)-1)*obj.pts+lcount) = lL;
                obj.L(1:rcount,(obj.auxUbic(obj.count+1)-1)*obj.pts+1:...
                    (obj.auxUbic(obj.count+1)-1)*obj.pts+rcount) = rL;
                
                obj.auxAlpha(1:lcount, obj.auxUbic(obj.count)) = lL\yL(1:lcount)';
                obj.auxAlpha(1:rcount, obj.auxUbic(obj.count+1)) = rL\yR(1:rcount)';
                
                obj.alpha(1:lcount, obj.auxUbic(obj.count)) = lL'\(obj.auxAlpha(1:...
                    lcount, obj.auxUbic(obj.count)));
                obj.alpha(1:rcount, obj.auxUbic(obj.count+1)) = rL'\(obj.auxAlpha(1:...
                    rcount, obj.auxUbic(obj.count+1)));
                div = -1; %stop the divide loop
                childModel = -1;
            end
            %relocate X Y
            obj.X(:, (obj.auxUbic(obj.count)-1)*obj.pts+1:...
                (obj.auxUbic(obj.count)-1)*obj.pts+obj.pts) = xL;
            obj.X(:, (obj.auxUbic(obj.count+1)-1)*obj.pts+1:...
                (obj.auxUbic(obj.count+1)-1)*obj.pts+obj.pts) = xR;
            obj.Y((obj.auxUbic(obj.count)-1)*obj.pts+1:...
                (obj.auxUbic(obj.count)-1)*obj.pts+obj.pts) = yL;
            obj.Y((obj.auxUbic(obj.count+1)-1)*obj.pts+1:...
                (obj.auxUbic(obj.count+1)-1)*obj.pts+obj.pts) = yR;
            obj.auxUbic(model) = 0; %parent model will not have points
        end
        
        function [pL,pR] = activation(obj, x, model)
            if obj.children(1,model) == -1 %return zeros when model has no children
                pL = 0;
                pR = 0;
                return
            end
            mP = obj.medians(model);
            xD = x(obj.overlaps(1,model)); %x value in cut dimension
            o = obj.overlaps(2,model); %half of the overlapping region
            if xD < mP-o/2
                pL = 1;
            elseif  xD >= mP-o/2 && xD <= mP+o/2 %if in overlapping
                pL = 0.5+(xD-mP)/(o);
                if(pL<=1e-12)%avoid numerical errors
                    pL=0;
                elseif(pL>=1-1e-12)
                    pL=1;
                end
            else
                pL = 0;
            end
            pR = 1-pL;
        end
        
        function update(obj,x,y)
            model = 1;
            while obj.children(1,model)~=-1 %if model is a parent
                %search for the leaf to asign the point
                [pL, ~] = obj.activation(x, model);
                if pL >= rand()
                    model = obj.children(1,model);%left child
                else
                    model = obj.children(2,model);
                end
                
%                 if pL > 0 && pL < 1
%                   disp('overlap')
%                 end
            end
            %add the model to the randomly selected model
            obj.addPoint(x,y,model)
        end
        
        function [out, outVar,outderiv, eta, beta,Lsigmax] = predict(obj,x)
            moP = zeros(2,900);% start from the root
            mCount = 1;
            moP(1,1) = 1;
            moP(2,1) = 1;
            while ~isequal( obj.children(1,moP(1,1:mCount)) , -1*ones(1,mCount) )
                for j=1:mCount
                    [pL, pR] = obj.activation(x,moP(1,j));
                    if pL > 0 && pR == 0
                        moP(1,j) = obj.children(1,moP(1,j));
                        moP(2,j) = moP(2,j)*pL;
                    elseif pR > 0 && pL == 0
                        moP(1,j) = obj.children(2,moP(1,j));
                        moP(2,j) = moP(2,j)*pR;
                    elseif pL>0 && pR>0
                        mCount = mCount + 1;
                        moP(1,mCount) = obj.children(2,moP(1,j));
                        moP(2,mCount) = moP(2,j)*pR;
                        moP(1,j) = obj.children(1,moP(1,j));
                        moP(2,j) = moP(2,j)*pL;
                    end
                end
            end
            out = 0;
            outderiv = 0;
            outVar = 0;
            outSigma = 0;

            gamma = 0;
            eta = 0;
            num = obj.xSize^(obj.xSize/2)*(max(obj.XteMax-obj.XteMin)^obj.xSize)*obj.count;
            logden = log(obj.delta)+obj.xSize*log(2*obj.tau);
            beta = 2*(log(num)- logden);
            Lsigmax = 0;
            for i=1:mCount
                model = moP(1,i);
                loctraindata = obj.X(:,(obj.auxUbic(model)-1)*obj.pts+1:...
                    (obj.auxUbic(model)-1)*obj.pts+obj.localCount(model));
                kstar = obj.kernel(loctraindata, x);
                pred = kstar' * ...
                    obj.alpha(1:obj.localCount(model),obj.auxUbic(model));
                xstardiff = ((loctraindata - x).')./(obj.sigmaL.^2)';
                meanJac = xstardiff.*kstar;
                meanderiv = meanJac.'*obj.alpha(1:obj.localCount(model),obj.auxUbic(model));
                
                v = obj.L(1:obj.localCount(model),(obj.auxUbic(model)-1)*obj.pts+1:...
                    (obj.auxUbic(model)-1)*obj.pts+obj.localCount(model))\...
                    obj.kernel(obj.X(:,(obj.auxUbic(model)-1)*obj.pts+1:...
                    (obj.auxUbic(model)-1)*obj.pts+obj.localCount(model)), x);
                predvar = (obj.kernel(x,x)-v'*v);
                outVar= outVar + (predvar+pred^2)*moP(2,i);
                out = out+pred*moP(2,i);
                outderiv = outderiv + meanderiv *moP(2,i);

                outSigma = outSigma + sqrt(predvar)*moP(2,i);
                gamma = gamma + moP(2,i)*(obj.Lmu(model)*obj.tau + sqrt(beta)*obj.Lsig(model)*obj.tau);
                Lsigmax = max(Lsigmax,obj.Lsig(model));
                
            end
            outVar = outVar - out^2;

            gamma = gamma + obj.Lf*obj.tau;
            eta = sqrt(beta)* outSigma + gamma;
        end

        function [brs] = ebound(obj,x)
            tic
            moP = zeros(2,1000);% line 1: active GPs, line 2: global probability
            mCount = 1; %number of GPs used for predictions
            moP(1,1) = 1; %start in root
            moP(2,1) = 1;
            %while all the GPs found are not leaves
            %get the GPs for prediction and thier global probabilites
            while ~isequal( obj.children(1,moP(1,1:mCount)) , -1*ones(1,mCount) )
                for j=1:mCount
                    [pL, pR] = obj.activation(x,moP(1,j));
                    if pL > 0 && pR == 0
                        moP(1,j) = obj.children(1,moP(1,j));
                        moP(2,j) = moP(2,j)*pL;
                    elseif pR > 0 && pL == 0
                        moP(1,j) = obj.children(2,moP(1,j));
                        moP(2,j) = moP(2,j)*pR;
                    elseif pL>0 && pR>0
                        mCount = mCount + 1;
                        moP(1,mCount) = obj.children(2,moP(1,j));
                        moP(2,mCount) = moP(2,j)*pR;
                        moP(1,j) = obj.children(1,moP(1,j));
                        moP(2,j) = moP(2,j)*pL;
                    end
                end
            end
            
            gamma = 0;
            brs = 0;
            outVar = 0;
            num = obj.xSize^(obj.xSize/2)*(max(obj.XteMax-obj.XteMin)^obj.xSize)*obj.count;
            logden = log(obj.delta)+obj.xSize*log(2*obj.tau);
            beta = 2*(log(num)- logden);
            
            %prediction: weigthing prediction and variance with proabilities
            for i=1:mCount
                model = moP(1,i);
                sig2m = obj.local_predict_var(x,model);
                gamma = gamma + moP(2,i)*(obj.Lmu(model)*obj.tau + sqrt(beta)*obj.Lsig(model)*obj.tau);
                outVar= outVar + moP(2,i)/sig2m;
                brs = brs + moP(2,i)/sqrt(sig2m);
            end
            gamma = gamma + obj.Lf*obj.tau;
            outVar = 1/outVar;
            brs = sqrt(beta)* brs * outVar + gamma;
        end
    end
end
