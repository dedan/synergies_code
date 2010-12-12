function Wpi = pcaica(mat,Ns)

%mat = observations X channels

builtin = false;


if builtin
    [Wpi, ~] = runica(mat', 'pca', 3, 'verbose', 'off');
    Wpi      = Wpi';
else
        
    %compute principle components
    %COEFF - channels X channels
    %SCORE - observations X channels
    [COEFF, SCORE] = princomp(mat);
    
    if Ns>1
        [W,SP]=runica(SCORE(:,1:Ns)', 'verbose', 'off');
        Wpi=COEFF(:,1:Ns)/(W*SP);
        D=diag(sign(sum(Wpi)));
        Wpi=Wpi*D;
    else
        Wpi=COEFF(:,1);
    end
end