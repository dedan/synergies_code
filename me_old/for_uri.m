

% What I want to do with this projections, is to find out wether the
% synergies that I found are somehow related to the synergies which have
% been found in a experiment without stimulation. Therefore I project the
% results of my computation, the activation patterns in a 11x148 matrix on
% the subspace which contains the synergies previously found (11x3).

% to get a measurment for the result I experimented with projections and
% random data, and found out that the ratio of the sums of the first norm
% of all the vectors in the range of the projections is the same as the
% ration between the dimensions of the subspaces. This I want to use as a
% measurement.

% I did it for different basises from several distributions. Now, what you
% see in the following code is what I really need. Orthonormal basis chosen
% from a modal distribution (similar to the 11x3 matrix I have) and my data
% will also be nonnegative.

% When I use the first norm, I come up with a nice ratio, but somehow only
% for normal distributed data. With nonnegative data I come to this offset
% and I don't know from where it comes and how to get rid of it.

% Another questions, which is not necessary for the project, but just out
% of curiosity. When I chose my basis from the normal distribution, I get
% this strange noise in the function for nonnegative data, do you know
% where it comes from ? (just uncomment the line where i create the space
% to see the difference)

n     = 1000;
dim   = 100;

rel_res = zeros(3,dim/2);

for sub_dim = 1:dim/2

   % construct space and complement and projections
  
   space = orth(rand(dim));

   part  = space(:,1:sub_dim);
   comp  = null(part');
   
   p = part  * part';
   p2 = comp  * comp';


   % create uniform and normal distributed data 
   pos = rand(dim,n);
   nor = randn(dim,n);

   % apply projections and back to unity space
   pos_on_b = space \ (p * pos);
   pos_on_o = space \ (p2 * pos);
   nor_on_b = space \ (p * nor);
   nor_on_o = space \ (p2 * nor);

   % make distribution
   num = 4;
   res = zeros(num,n);
   for i = 1:n
      res(1,i) = sum(abs(pos_on_b(:,i)));
      res(2,i) = sum(abs(pos_on_o(:,i)));
      res(3,i) = sum(abs(nor_on_b(:,i)));
      res(4,i) = sum(abs(nor_on_o(:,i)));
   end

%    for i = 1:n
%       res(1,i) = norm(abs(pos_on_b(:,i)));
%       res(2,i) = norm(abs(pos_on_o(:,i)));
%       res(3,i) = norm(abs(nor_on_b(:,i)));
%       res(4,i) = norm(abs(nor_on_o(:,i)));
%    end


   % ratio between the sum of the first norm of the vectors on one space
   % and sum of first nor of the vectors on its complement
   rel_res(1,sub_dim) = sum(res(1,:)) / sum(res(2,:));
   rel_res(2,sub_dim) = sum(res(3,:)) / sum(res(4,:));
   rel_res(3,sub_dim) = sub_dim / (dim - sub_dim);

end

figure(5);
plot(rel_res')

hold on
plot(rel_res(1,:) - rel_res(3,:));