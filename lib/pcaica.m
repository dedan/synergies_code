function [Wpi,loadings, resid]=pcaica(mat,Ns)

%mat = observations X channels


m=mean(mat);
lm=size(mat,1);

%compute principle components
%COEFF - channels X channels
%SCORE - observations X channels  
[COEFF, SCORE] = princomp(mat);

%recons=(COEFF(:,1:Ns)*SCORE(:,1:Ns)')'+ones(lm,1)*m;

if Ns>1
    [W,SP]=runica(SCORE(:,1:Ns)','extended',1, 'verbose', 'off');
    Wpi=COEFF(:,1:Ns)*inv(W*SP);
    D=diag(sign(sum(Wpi)));
    Wpi=Wpi*D;   
    loadings=inv(D)*(W*SP)*SCORE(:,1:Ns)';
else
    Wpi=COEFF(:,1);
    loadings=SCORE(:,1)';
end

SS = sum(sum(mat));
dif = mat - ones(lm,1)*m-(COEFF(:,1:Ns)*SCORE(:,1:Ns)')';
resid = (sum(sum(abs(dif)))/SS) * 100;
