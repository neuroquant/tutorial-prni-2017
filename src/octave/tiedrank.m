## Copyright (C) 2017 
## 
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*- 
## @deftypefn {Function File} {@var{retval} =} tiedrank (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author:  Manjari Narayan
## Created: 2017-05-24

function [xranks tieadj] = tiedrank (x)

    if(ndims(x)==2)
        [dim1 dim2] = size(x);
        if(dim1>1 && dim2>1)
            disp('WARNING. This tiedrank function only takes vectors.');
        end
    else
        disp('WARNING. This tiedrank function only takes vectors as input');
    end

    xranks = ranks(reshape(x,[length(x) 1]),1);
    %% To normalize tiedcorrection for length. Not recommended. Not typical use case. 
    % tieadj = tiedcorrection(xranks)/(length(x)*(length(x)^2-1));
    tieadj = tiedcorrection(xranks);

endfunction



function correction =  tiedcorrection(ranks)
	
	u_ranks = unique(ranks); 
	n_unique = length(u_ranks);
	correction = 0;
	
	if(n_unique<length(ranks)) 
		for grp_no=1:n_unique
			tmp_idx = find(ranks==u_ranks(grp_no)); 
			n_ties = length(tmp_idx); 
			correction = correction + (n_ties^3-n_ties);
		end
	end
end
