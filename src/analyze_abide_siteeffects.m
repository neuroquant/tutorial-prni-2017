function results = analyze_abide_siteeffects(varargin)
% Compute and return conditional correlations, nuisance correlations, and NSR on all subjects x sites.
% INPUT
% (optional)
%   - nuisance_type: 'within_subject' or 'within_site' nuisance
% 
    datadir = '../data';
	studydata = load(fullfile(datadir,'ABIDE_controlData_7Sites.mat')); 
    studydata.yeo_order = load(fullfile(datadir,'cc200Map2Communities'));
    [~,reorder_communities] = sort(studydata.yeo_order.mapToCommunities,'ascend'); 
    
    % Save filename
    exportfun = @(fname)(print('-dpng','-r300',fname));
	fname = ['tmp' filesep datestr(now,'dd-mmm-yyyy-HHMM')]; 
	if(~exist(fname,'dir'))
		mkdir(fname)
	end
    
	results = {};
    results.observed = [];
    results.tw_observed = [];
    results.denoised = [];
    results.nuisance = [];
    results.nsr = [];
    sitelabels = [];
    denoised = [];
    nuisance = [];
    observed = [];
    tw_observed = [];
    nsr = [];
    
    % Site Group Options
    sitegroup = 2;
    switch sitegroup
    case 1
        sites = [3 5];
    case 2
        sites = [1 2 4 6 7];
    case 3
        sites = [1:7];
    end
    nsites = length(sites);

    if(nargin==0)
        nuisance_type = 'within_subject'; 
    else
        nuisance_type = varargin{1}; 
    end

    for studyno=sites;  
        sitelabels = cat(2,sitelabels, ...
                ones(1,length(studydata.data{studyno}.subIDs))*studyno);
        
        if(strcmp(nuisance_type,'within_site'))         
            Y = squeeze(mean( ...
                studydata.data{studyno}.signals(:,reorder_communities,:),2))'; 
        end
                 
    	for ii=1:length(studydata.data{studyno}.subIDs)
            
            if(mod(ii,5)==0)
                disp(sprintf('Site %s, Processing Subject %d ...', ...
                        studydata.data{studyno}.dataName,ii));
            end
    		X = squeeze(studydata.data{studyno}.signals(ii,reorder_communities,:))';
            
            if(strcmp(nuisance_type,'within_subject'))
                Y = mean(X,2); 
            end
            
%             %Check subject level conditional correlation
%             if(mod(ii,10))
%                 demo_conditional_correlation(X,Y);
%             end
            
            try 
                % Standard observed correlation               
                results(studyno).observed(ii,:,:) = standard_correlation(X);
                observed = cat(1,observed,results(studyno).observed(ii,:,:)); 
                results(studyno).tw_observed(ii,:,:) = standard_correlation_sn(X);  
                tw_observed = cat(1,observed,results(studyno).tw_observed(ii,:,:)); 
                % Denoised Correlation             
                output = conditional_correlation(X,Y);
                results(studyno).denoised(ii,:,:) = output.corr;
                denoised = cat(1,denoised,results(studyno).denoised(ii,:,:)); 
                % Nuisance Correlation
                results(studyno).nuisance(ii,:,:) = output.nuisance;
                results(studyno).nsr(ii) = output.NSR; 
                nuisance = cat(1,nuisance,results(studyno).nuisance(ii,:,:));                 
                nsr = cat(1,nsr,output.NSR); 
                
            catch me
                disp(me)
                disp(me.stack)
            end
    	end
    end
	 

	if(exist('brewermap'))	
		colormapfun = @()(flipud(brewermap(length(colormap),'RdYlBu')));
		close all;
	else
		colormapfun = @winter;
	end
	
    figure(1);
	set(gcf,'Position',[10 450  1350  375],'PaperPosition', [.5 1.5 12.0 7.0]);
    subplot(1,3,1); 
    site_effect{1} = detect_site_effect(observed,sitelabels); 
    imagesc(real(site_effect{1}.similarity)); 
    colormap(colormapfun()); colorbar; axis image; 
    title('Similarity Matrix (Observed)','fontsize',24);
    xlabel(sprintf('(wit,bet,rat,test) = (%.2f, %.2f, %.2f)', ...
                    site_effect{1}.within,  ...
                    site_effect{1}.between, ...
                    site_effect{1}.ratio));
    set(gca,'fontsize',16);
    
    subplot(1,3,2); 
    site_effect{2} = detect_site_effect(nuisance,sitelabels); 
    imagesc(real(site_effect{2}.similarity)); 
    colormap(colormapfun()); colorbar; axis image; 
    title('Similarity Matrix (Nuisance)','fontsize',24);    
    xlabel(sprintf('(wit,bet,rat,test) = (%.2f, %.2f, %.2f)', ...
                    site_effect{2}.within,  ...
                    site_effect{2}.between, ...
                    site_effect{2}.ratio));
    set(gca,'fontsize',16);
    
    
    subplot(1,3,3); 
    site_effect{3} = detect_site_effect(denoised,sitelabels); 
    imagesc(site_effect{3}.similarity); 
    colormap(colormapfun()); colorbar; axis image; 
    title('Similarity Matrix (Denoised)','fontsize',24);    
    xlabel(sprintf('(wit,bet,rat,test) = (%.2f, %.2f, %.2f)', ...
                    site_effect{3}.within,  ...
                    site_effect{3}.between, ...
                    site_effect{3}.ratio));
    set(gca,'fontsize',16);
    
	exportfun(fullfile(fname,[mfilename '1']));
    

    function Sigma = standard_correlation(X)
        % Usual standard correlation matrix


        [Sigma results] = covariance.mle_sample_covariance(X, ...
                                                    struct('standardize',false));


    end


    function Sigma = standard_correlation_sn(X)
        % Automatically applies row-first successive norm
        %standardize.successive_normalize(X');


        [Sigma results] = covariance.mle_sample_covariance(X, ...
                                                    struct('standardize',true));

    end


    function site_effect = detect_site_effect(X,sitelabels)

        corrfun = @corr;
        upper_idx = find(reshape(triu(ones(size(X,2), size(X,3)),1), [1 size(X,2)^2])); 
        X = reshape(X,[size(X,1) size(X,2)*size(X,3)]);
        tmpSimilarity = corrfun(X(:,upper_idx)'); 
        [hom sep mw] = compareWithinAndBetweenGroupsSim(tmpSimilarity,sitelabels);

        site_effect.similarity = tmpSimilarity;           
        site_effect.within = hom; 
        site_effect.between = sep; 
        site_effect.stat = mw;  
        site_effect.ratio = (hom-sep)/(hom+sep); 


    end

    function [nuisance] = get_shared_nuisance(Y); 

       % Y is time-series x subjects
       Yz = zscore(Y')'; 
       [U S V] = svd(Y);  
       nfactors = min(size(Y,2),3);   
       nuisance = U(1:nfactors,:); 
       nuisance = reshape(nuisance,[size(Y,1) nfactors]); 

    end


    function output =  conditional_correlation(X,Y)
        % Only uses usual column standardize (i.e. correlation)
        output = struct();

        [Sigma results] = covariance.sample_conditional_correlation(X, ...
                                        struct('verbose',false,...
                                                'nuisance',Y) ...	
                                                );

        output.corr = results.SigmaCov;
        output.nuisance = results.nCov; 
        output.corr2 = covariance.mle_sample_covariance(results.X_perpY, ...
                                                    struct('standardize', false));
        output.NSR = results.NSR;

    end

end