function [X output] = successive_normalize(X,varargin)
	%SUCCESSIVE_NORMALIZE implements normalization strategy inspired 
	% by Olshen and Rajaratnam (2010); Allen and Tibhsirani (2011); 
	% The symmetric variant is new. 
	% 
	% USAGE: [X output] = successive_normalize(X)
	% INPUT
	% 	- X : a n x p transposable data matrix
	% 	- options: (optional)
	% 	- options.method: 'sym' or 'original'. 
	% 		Original method converges to different solutions depending on row or column first. 'Sym' enforces symmetry at every iteration. 
	% OUTPUT: 
	% 	- X : normalized n x p data matrix
	% 	- output : Contains original row,column parameters
	% 
	%  Copyright 2017, Manjari Narayan
	%  BSD-2 Clause License
	% 

	stop_err = Inf;
	tol = 1e-3;
	n_iter = 5; 
	k = 0;
	
	narginchk(1,2);
	if(nargin==2)
		options = varargin{1}; 
	else
		options = [];
	end
	
	if(isfield(options,'method'))
		method = options.method;
	else
		method = [];
	end
	if(isempty(method))
		method = 'original';
	end


	while ((stop_err > tol) && (k<=n_iter))
	
		if(k==0)
			
			switch method 
				
			case 'sym'
				
				[Xc_colpolish mu_c sig_c] = standardize.standardize_cols(X); 
				Xc_rowpolish = standardize(Xc_colpolish');	
			
				[Xr_rowpolish mu_r sig_r] = standardize.standardize_cols(X');	
				Xr_colpolish = standardize(Xr_rowpolish'); 
			
				X_colpolish = (Xr_colpolish + Xc_colpolish)/2;
				X_rowpolish = (Xr_rowpolish + Xc_rowpolish)/2;
				
				sym_err(k+1) = norm(Xr_colpolish-Xc_rowpolish','fro');
				
				
			otherwise
				
				[X_colpolish mu_c sig_c] = standardize.standardize_cols(X); 
				[~, mu_r, sig_r] = standardize.standardize_cols(X'); 
				
				X_rowpolish = standardize.standardize_cols(X_colpolish');	
				
				sym_err = [];				
				
			end
			
			oldXc = X;
			oldXr = X';
		else

			switch method 
				
			case 'sym'
				
				Xc_colpolish = standardize.standardize_cols(X_rowpolish'); 
				Xc_rowpolish = standardize.standardize_cols(Xc_colpolish');	
			
				Xr_rowpolish = standardize.standardize_cols(X_colpolish');	
				Xr_colpolish = standardize.standardize_cols(Xr_rowpolish'); 
			
				X_colpolish = (Xr_colpolish + Xc_colpolish)/2;
				X_rowpolish = (Xr_rowpolish + Xc_rowpolish)/2;
				
				sym_err(k+1) = norm(Xr_colpolish-Xc_rowpolish','fro');
				
				
			otherwise
				
				[X_colpolish] = standardize.standardize_cols(X_rowpolish'); 
				[X_rowpolish] = standardize.standardize_cols(X_colpolish');	
											
			end
			
		end
	
		stop_err = (.5*norm(X_colpolish-oldXc,'fro') + ...
			 		.5*norm(X_rowpolish-oldXr,'fro'))/numel(X);
	
		if(k==n_iter)
			if(options.verbose)				
				disp(sprintf(...
				'Iter:%d, stop_err:%.8f. Exiting loop.',k,stop_err));
			end
		end			
	
	
		k = k+1; 
		err(k) = stop_err;
		
		oldXr = X_rowpolish;
		oldXc = X_colpolish;

	end

	X = (X_colpolish + X_rowpolish')/2;
	output = struct('mu_r', mu_r, ...
					'sig_r', sig_r, ...
					'mu_c', mu_c, ...
					'sig_c', sig_c, ...
					'iter', k, ...
					'err', err, ...
					'sym_err', sym_err ...
					);	

end

% function [X mu sig]  = standardize(X)
%     %STANDARDIZE is a helper function for SUCCESSIVE_NORMALIZE

%     [X mu sig] = standardize.standardize_cols(X);

% end