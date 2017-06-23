function [Sigma results] = sample_conditional_correlation(X,options)
%SAMPLE_CONDITIONAL_CORRELATION
% 
% Sigma = Sigma_obs - Sigma_nuisance
% 
% Obtained by 
% Schur complement using covariance decomposition
% Residuals after regressing out confounds and modified cholesky decomposition 
% 	
% USAGE: 
% 
% INPUT
% 	- X : A single n x p Data matrix
% 	- options.scaled: Columns of X scaled to unit diagonals
% 	- options.verbose: Plot and print outputs
% 	- options.outputdir: Defaults to 'tmp/'
% 	- options.filename: Defaults to 'conditionalCov'
% 
% OUTPUT
% - results.Sigma - Denoised Correlation Matrix
% - results.nCov - Nuisance Covariance Matrix
% - results.NSR - Nuisance to Observed Signal Ratio (trace(nCov)/trace(Sigma_obs))
% - results.X_perpY - Projection of Data onto orthogonal complement of nuisance Y
% - results.nCorr - Rescaled Nuisance Correlation Matrix	
	
	if(~exist('options','var'))
		options = []
	end
	if(isempty(options))
		options = []
	end
	if(~isfield(options,'verbose'))
		options.verbose = false;
	end
	if(~isfield(options,'outputdir'))
		options.outputdir = fullfile('tmp',datestr(now,'dd-mmm-yyyy-HHMM')); 
	end	
	if(~isfield(options,'filename'))
		options.filename = mfilename;
	end
	
	% TBD: Add process_options or processArgs
	% ~isfield(options,'verbose')
	
	if(~isfield(options,'corrfun'))
        %options.corrfun = @(X)(cov(X));
        options.corrfun = @(X)(covariance.mle_sample_covariance( ...
                                    X, ...
                                    struct('standardize',false) ...
                          ));
	end
	
	function Sig = my_xcorr(X,y,options)
	% Returns p x q cross correlation between n x p and n x q	
	
		[n1 p] = size(X); 
		[n2 q] = size(y); 			
		assert(n1==n2,'X and y need to have same row length')
		assert(n2>q,'Too many nuisance variables'); 
	
		if(~exist('options','var'))
			options = [];
		end
		if(~isfield(options,'standardize'))	
			options.standardize = false;
		end
		
		if(options.standardize)
            % Xt = standardize.successive_normalize(X');
            % X  = Xt';
            X = standardize.standardize_cols(X);	                 	
            if(q>1)
                %y = standardize.successive_normalize(y);
    			[~, mu_y, sig_y] = standardize.standardize_cols(y); 
            else
    			[~, mu_y, sig_y] = standardize.standardize_cols(y); 
            end
			Sig = X'*y/n1;
            sig_y(sig_y<=0) = 1;
            Sig = Sig*diag(1./sig_y); 
					
		else
			[~, mu_x, sig_x] = standardize.standardize_cols(X); 
			[~, ~, sig_y] = standardize.standardize_cols(y); 
			%y = zscore(y);
			
			if(q>1)
				Sig = bsxfun(@minus,X,mu_x)'*y/n1; 
			elseif(q==1)
				Sig = bsxfun(@minus,X,mu_x)'*y/n1; 
			end
			
            sig_x(sig_x<=0) = 1;
			Sig = diag(1./sig_x)*Sig;
			%Sig = diag(1./sig_x)*Sig*diag(1./sig_y);			
            			
		end
		
	end
	
    proj = @(X)(X*pinv(X'*X)*X');
	if(~isfield(options,'xcorrfun'))
		%options.xcorrfun = @(X,y)(my_xcorr(X,y,struct('standardize',true)));
        options.xcorrfun = @(X,y)(proj(y)*standardize.standardize_cols(X)); 
	end
	
	if(~isfield(options,'nuisance'))
		options.nuisance = mean(X,2); 
	end
	
	y = options.nuisance;    
	Sig_xx = options.corrfun(X); 		
	Sig_xy = options.xcorrfun(X,y); 
	Sig_yy = options.corrfun(y);
    if(size(Sig_xy,2)==size(Sig_yy,1))
        tmpSigmaL = Sig_xy*pinv(Sig_yy)*Sig_xy';
        [LV LD] = eig(tmpSigmaL);
        [sortLD, LDidx] = sort(diag(LD),'descend'); LDidx = LDidx(real(sortLD)>0);    
        SigmaL = LV(:,LDidx)*real(LD(LDidx,LDidx))*LV(:,LDidx)';
        SigmaL = real(SigmaL);
    else
        tmpSigmaL = cov(Sig_xy); 
        SigmaL = tmpSigmaL;
    end    

    % Remove low rank component
    Sigma = Sig_xx - SigmaL;
    varcorr = @(D,Sig)(diag(1./D)*Sig*diag(1./D));
    D = real(diag(Sigma)); 
    D(D<=0) = 1; 
	
	results.SigXX = Sig_xx;
	results.SigXY = Sig_xy;
	results.SigYY = Sig_yy;
	results.nCov  = SigmaL;
    results.Sigma = varcorr(sqrt(D),Sigma);
	results.SigmaCov = Sigma; %
    results.nCorr = varcorr(sqrt(diag(SigmaL)),SigmaL);
	results.NSR = trace(SigmaL)/trace(Sig_xx);
    results.SNR = trace(Sigma)/trace(SigmaL);
    results.Y = y;
	
	function Xy =  Xy_orthogonalize(X,Y)
	% Nuisance signal regression
		
		% p = size(X,2);
		% Xy = zeros(size(X));
		%
		% for ii=1:p
		% 	gsProj = gsr'*X(:,ii)/(gsr'*gsr)*gsr;
		% 	Xy = X(:,ii) - gsProj;
		% end
		
		mu = mean(X); 
		Xcen = bsxfun(@minus,X,mu); 
        
		proj = @(X)(X*pinv(X'*X)*X');
		Xy = Xcen - proj(Y)*Xcen;
		Xy = bsxfun(@plus,Xy,mu);
		 
	end

	results.X_perpY = Xy_orthogonalize(standardize.standardize_cols(X),options.nuisance);
	
	
	if(options.verbose)
		% Only perform for comparison purposes;
		nuisance = options.nuisance; 				
		Xgsr = Xy_orthogonalize(standardize.standardize_cols(X),nuisance); 
		
    	if(exist('brewermap'))	
    		colormapfun = @()(brewermap(length(colormap),'RdYlBu'));
    		close all;
    	else
    		colormapfun = @winter;
    	end
        
		figure; 
		subplot(2,2,1); 
		imagesc(Sig_xx); colorbar; axis equal image; 
		title('Correlation'); 
		subplot(2,2,2); 
		imagesc(results.Sigma); colorbar; axis equal image; 
		title('Conditional Corr.'); 
		subplot(2,2,3);  
		imagesc(results.nCov); colorbar; axis equal image; 
		title('Nuisance Factor'); 
		subplot(2,2,4); 
		imagesc(options.corrfun(Xgsr)); colorbar; axis equal image; 
		title('Correlation (NSR)');
        fullfile(options.outputdir,options.filename)
		print('-dpdf','-r150',fullfile(options.outputdir,options.filename)); 		pause(5)	
		close all;
	end	
	
	
	
end