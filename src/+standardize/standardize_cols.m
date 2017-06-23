function [X mu sig]  = standardize_cols(X,varargin)
	%STANDARDIZE_COLS z-scores the columns of X
	% Used as a is a helper function for SUCCESSIVE_NORMALIZE
	% 
	% USAGE: [X mu sig] = standardize_cols(X)
	% 		 [X mu sig] = standardize_cols(X,'center')
	% 		 [X mu sig] = standardize_cols(X,'scale')
	% 
	% 
	% SEE ALSO SUCCESSIVE_NORMALIZE
	
	narginchk(1,2);
	
	if(nargin==1)
		method = 'all';
	elseif(nargin>=1)
		method = varargin{1};
	end
	
	mu 	= mean(X); 		 
	sig = std(X,1); 
	assert(any(sig<1e-5)==0,'Check data matrix for near constant or 0 rows or colums'); 	
	
	switch method
		
	case 'center'
		X 	= bsxfun(@minus,X,mu); 
	case 'scale'
		X 	= bsxfun(@minus,X,mu); 		
		X 	= bsxfun(@rdivide,X,sig);
		X 	= bsxfun(@plus,X,mu./sig); 				
	otherwise		
		X 	= bsxfun(@minus,X,mu); 
		X 	= bsxfun(@rdivide,X,sig); 		
	end 	
	
end