function res = circ_std(vec_rad)

% circular standard deviation

% implementation from the book: Statistical Analysis of Circular Data of
% N. I. Fisher



n = length(vec_rad);
C = sum(cos(vec_rad));
S = sum(sin(vec_rad));

% resultant length
R_square = C^2 + S^2;

% mean resultant length
R_mean = sqrt(R_square)/n;


res = sqrt(-2 * log(R_mean));



end