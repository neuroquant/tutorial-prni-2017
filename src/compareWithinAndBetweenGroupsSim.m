function [ hom sep mw wbRat] = compareWithinAndBetweenGroupsSim( simMat, grouping )
%compareWithinAndBetweenGroupsSim compares similaritiy between elements
%from same group to that between elements from different groups (can be
%used to examnie correspondence with a known organization, site
%distribution etc.
%   Detailed explanation goes here
%Input: simMat - a symmetric similarity matrix in which entry i,j and j,i
%                both contain the similarity between element i and element j
%       grouping - a known partition of the elements which is to be
%                  examined (0 is treated as non group)
%Output: stats of wheather within group similarity is greater than between group similarity :
%  hom - homogeneity (average pairwise similarities) within groups 
%  sep - separation between groups (average pairwise similarities for pairs that are from different groups)
%  mw - stats from a 1 tailed Mann�Whitney (Ranksum) test
% 
% 2017, Adi Maron Katz
withinSimVals=[];
betweenSimVals=[];
numElems=size(simMat,1);
for i=1:(numElems-1)
    for j=(i+1):numElems
        if(grouping(i)>0 && grouping(i)==grouping(j))
            withinSimVals=[withinSimVals,simMat(i,j)];
        else
            betweenSimVals=[betweenSimVals,simMat(i,j)];
        end
    end
end
hom=mean(withinSimVals);
sep=mean(betweenSimVals);
wbRat=mean(withinSimVals)/mean(betweenSimVals);

 % [mw.p,mw.h,mw.stats] =
 % ranksum(withinSimVals,betweenSimVals,'tail','right');%This function is
 % The above two lines have been replaced with the code below for Octave compatibility:
 mwRes=mwwtest(betweenSimVals,withinSimVals);
 direction=sign(mean(withinSimVals)-mean(betweenSimVals));
 mw.stats.zval=direction*mwRes.Z;
 if(direction>0)
     mw.p=mwRes.p(1);
 else
     mw.p=1-mwRes.p(1);
 end
 if(mw.p<0.05)
     mw.h=1;
 else
     mw.h=0;
 end
 mw.stats.ranksum=mwRes.W(1);
end
 %[mw.p, mw.stats.zval] = u_test (withinSimVals,betweenSimVals, '<');
 %if(mw.p<0.05)
 %    mw.h=1;
 %else
 %    mw.h=0;
 %end


