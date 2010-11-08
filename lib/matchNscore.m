function [scores, index]=matchNscore(W1,W2)

%This function matches column from W1 and W2. 
%scores gives the absolute dot product between columns of W1 and their
%matched columns in W2. index gives the identity of the columns in W2 that
%were matched to the columns of W1. For example, index(1)=3 means that the
%first column of W1 was matched with the 3rd column of W2. 

s=size(W1,2);
scores=zeros(1,s);
index=zeros(1,s);
WW1=normc(W1);
WW2=normc(W2);
Scores=abs(WW1'*WW2);
for i=1:s
    [m,irow]=max(Scores);
    [m,icol]=max(m);
    irow=irow(icol);
    Scores(irow,:)=-Inf;
    Scores(:,icol)=-Inf;
    index(irow)=icol;
    scores(irow)=m;
end