function [R results] = var_corr(Sigma,options)
%VAR_CORR - returns variance correlation decomposition of any symmetric, positive semi-definite matrix (usually covariance matrix) 
% 				Sigma = Gamma^{1/2} R Gamma^{1/2}
% USAGE
% 
% INPUTS
%
% OUTPUTS

	if(~exist('options','var'))
		options = [];
	end
	
	if(isfield(options,'verbose'))
		verbose=options.verbose;
	else
		options.verbose = true;
		verbose = options.verbose;
	end

	import covariance.check_symposdef
	[issym isposdef results.symposdef] = check_symposdef(Sigma);
	try
		assert(issym>=1,'Not symmetric')
	catch me
		if(verbose)
			disp(me)
			results.symposdef.sym_err
		end		
	end
	try	
		assert(isposdef>=1,'Negative Eigenvalues')
	catch me
		if(verbose)
			disp(me)
			min(results.symposdef.eigs)
		end
	end
		
	results.input.Sigma = Sigma;
	results.input.options = options;
	results.issym = issym; 
	results.isposdef = isposdef;
	
	if(issym==1)
		if(verbose)
			disp('Frob. Error: A-At')
			results.symposdef.sym_err
			disp('Symmetricizing input');
		end
		Sigma = (Sigma + Sigma')/2;
	end
	
	p 	= length(Sigma); 

	gam = sqrt(diag(Sigma));
	d 	= 1./gam;

	R 	= diag(d)*Sigma*diag(d);
	
	results.var = d;
	results.corr = R; 
	
end