function results =  demo_successive_norm(X,varargin)
	
	opts = struct();
	opts.exportfig = false;
    opts.verbose = false;
	opts.exportfun = @(fname)(print('-dpng','-r300',fname));	
	exportfig = opts.exportfig;
	exportfun = opts.exportfun;
	
	results = {};
	
	results{1}.method = 'Column Standardize';
	results{1}.output = standard_correlation(X);
	
	results{2}.method = 'Row-Column Standardize';
	results{2}.output = standard_correlation_sn(X);

    if(opts.verbose)		
        disp(sprintf('Frob. MSE: Sigma_rc - Sigma_c = %.3f',  ...
            norm(abs(results{2}.output.corr-results{1}.output.corr),'fro')));

        disp(sprintf('Frob. SSE: Sigma_rc - Sigma_c = %.3f', ...	
            sum(sum((results{2}.output.corr-results{1}.output.corr).^2)) ));
	end
    
    if(exist('brewermap'))	
		colormapfun = @()(brewermap(length(colormap),'RdYlBu'));
		close all;
	else
		colormapfun = @winter;
	end
    
    histogram = @(x,varargin)(hist(x,linspace(-1.0,1.0,100),1));
    
	figure('Position',[50 100 800 650]);
    upper_idx = find(reshape(triu(ones(size(results{1}.output.corr)),1),  ...
                            [1 numel(results{1}.output.corr)])); 	
	subplot(1,2,1); 
	imagesc(results{1}.output.corr); axis equal image;
	colormap(colormapfun()); 
	title(results{1}.method);
    set(gca,'fontsize',16);
	subplot(1,2,2); 
	imagesc(results{2}.output.corr); axis equal image;
	colormap(colormapfun()); 
	title(results{2}.method);
    set(gca,'fontsize',16);

	figure('Position',[50 100 800 650]);
	subplot(1,2,1); 	
	histogram(results{1}.output.corr(upper_idx),'Normalization','pdf'); 
    ylim([-1.2 1.2]);axis tight; set(gca,'fontsize',16);

	%axis equal image;
	title(results{1}.method); xlabel('correlation'); ylabel('pdf')
	subplot(1,2,2); 
	histogram(results{2}.output.corr(upper_idx),'Normalization','pdf'); 
    ylim([-1.2 1.2]);axis tight; set(gca,'fontsize',16);

	%axis equal image;
	title(results{2}.method);xlabel('correlation'); ylabel('pdf')
	
	fname = ['tmp' filesep datestr(now,'dd-mmm-yyyy-HHMMSS')]; 
	if(~exist(fname,'dir'))
		mkdir(fname)
	end
    if(exportfig)
	    exportfun(fullfile(fname,mfilename));	
    end
	
end


function output = standard_correlation(X)
	% Usual standard correlation matrix
	
	output = struct();
	
	[Sigma results] = covariance.mle_sample_covariance(X, ...
												struct('standardize',false));
	
	output.corr = Sigma;
	
end


function output = standard_correlation_sn(X)
	% Automatically applies row-first successive norm
	%standardize.successive_normalize(X');
	
	output = struct();
	
	[Sigma results] = covariance.mle_sample_covariance(X, ...
												struct('standardize',true));
	
	output.corr = Sigma;
	
end

