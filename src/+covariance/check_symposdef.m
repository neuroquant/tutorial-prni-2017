function [issym isposdef results] = check_symposdef(A)
	
	[p1 p2] = size(A); 	
	assert(p1==p2,'Not square matrix'); 
	
	vecA = A(find(triu(ones(p1,p2),1))); 
	vecB = A(find(tril(ones(p1,p2),-1)));
	sym_err = sum((vecA-vecB).^2)/length(vecA);
	issym  	= sym_err<.5;
	if(issym & sym_err> 1e-2)
		warning('Numerical Asymmetry. Symmetricize by A + A^T'); 
		issym = 1;
	else
		issym = 2;
	end
	
	[V D] = eig(A);
	tol = 1e-5;
	if(min(D)>tol)
		isposdef = 2; % Positive Definite
	elseif(abs(min(D))<=tol)
		isposdef = 1; % Positive Semi-Definite
	else
		isposdef = 0; % Negative Definite
	end
	
	results.sym_err = sym_err;
	results.eigs = D;
	results.eigv = V;
end