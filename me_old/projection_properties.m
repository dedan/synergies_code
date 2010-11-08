
% in this file I do some experiments to find the statistical properties of
% projections on subspaces. Dependent on the distributions from which I
% chose the basis and the data

% still open questions are the offset that appears when using a basis and
% data from the unimodal distribution
% also the noise in using a basis from the normal distribution with data
% from the unimodal is an open question


%% 3 d example

n = 100;
dim = 3;

u_space = eye(3);
u_2dim  = u_space(:,1:2);
u_1dim  = u_space(:,3);
pu = u_2dim * u_2dim';
pu2 = u_1dim * u_1dim';


r_space = orth(rand(3));
r_2dim  = r_space(:,1:2);
r_1dim  = null(r_2dim');
pr = r_2dim * r_2dim';
pr2 = r_1dim * r_1dim';


rr_space = orth(randn(3));
rr_2dim  = rr_space(:,1:2);
rr_1dim  = null(rr_2dim');
prr = rr_2dim * rr_2dim';
prr2 = rr_1dim * rr_1dim';


pos = rand(dim,n);
nor = randn(dim,n);
figure(1)
hold off
scatter3 (pos(1,:), pos(2,:),pos(3,:))
hold on

% on_u_b = pu * pos;
% on_u_o = pu2 * pos;
% scatter3 (on_u_b(1,:), on_u_b(2,:),on_u_b(3,:))
% scatter3 (on_u_o(1,:), on_u_o(2,:),on_u_o(3,:))

% on_r_b = pr * pos;
% on_r_o = pr2 * pos;
% scatter3 (on_r_b(1,:), on_r_b(2,:),on_r_b(3,:))
% scatter3 (on_r_o(1,:), on_r_o(2,:),on_r_o(3,:))
% 
% figure(2)
% hold off
% scatter3 (pos(1,:), pos(2,:),pos(3,:))
% hold on
% on_r_b = inv(r_space) * pr * pos;
% on_r_o = inv(r_space) * pr2 * pos;
% scatter3 (on_r_b(1,:), on_r_b(2,:),on_r_b(3,:))
% scatter3 (on_r_o(1,:), on_r_o(2,:),on_r_o(3,:))

% rr_space = rr_space(:,[3 1 2]);

on_rr_b = prr * pos;
on_rr_o = prr2 * pos;
scatter3 (on_rr_b(1,:), on_rr_b(2,:),on_rr_b(3,:))
scatter3 (on_rr_o(1,:), on_rr_o(2,:),on_rr_o(3,:))

figure(2)
hold off
scatter3 (pos(1,:), pos(2,:),pos(3,:))
hold on
on_rr_b = rr_space \ (prr * pos);
on_rr_o = rr_space \ (prr2 * pos);
scatter3 (on_rr_b(1,:), on_rr_b(2,:),on_rr_b(3,:))
scatter3 (on_rr_o(1,:), on_rr_o(2,:),on_rr_o(3,:))


%% relation between distributions and dimensions


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

figure(3);
plot(rel_res')

hold on
plot(rel_res(1,:) - rel_res(3,:));


%% check the upper limit

% after I found out somethin about the relation between subspace dimensions
% and the distributions of first norms of projection ranges, I have to find
% a upper limit for the quantification of my results.
% therefore I want to construct data, which projects only onto the basis of
% the nonevoked synergists, while still having the same statistical
% properties as my data.
% in previous experiments I've seen, that I don't have to normalize the
% data what makes it possibly easier.


% ich suche mir ein basis die dat aufspannen kann, schaue dann auf die
% eigenschaften (verteilungen) der faktoren, die ich benötige um aus dieser
% basis den datenraum aufzuspannen.
% dann erzeuge ich faktoren von genau diesen verteilungen, mit welchen ich
% aus meiner basis (nonevoked) einen datenraum schaffe, der komplett auf
% die basis projeziert werden kann. dies sollte mir ein upper limit geben
% das klappt so nicht, ich habe ja jetzt elf faktoren mit unterschiedlichen
% eigenschaften, welche soll ich denn jetzt für linearkombinationen der
% basis benutzen. 
% also alter ansatz, ich versuche die faktoren zu finden, die mir aus
% linearkombinationen der basis möglichst genau meine daten beschreiben.


load ~/Documents/uni/yifat_lab/results/30122008__213/data.mat
dat = fin_res(1).dat';

grouped  = group(res, 'Wopt');
c        = grouped.center;
basis    = orth(c');
dim      = size(basis,2);
mean_dat = mean(dat,2);
std_dat  = std(dat,0,2);
x        = linsolve(basis,dat);
mean_x   = mean(x,2);
std_x    = std(x,0,2);
x_basis  = NaN(size(x));

dist_length = 100000;
% verteilungen erzeugen
for i = 1:size(x,1)
   distribution = (randn(1,dist_length) * std_x(i)) + mean_x(i);
   perms        = randperm(dist_length);
   x_basis(i,:) = distribution(perms(1:size(x,2)));
end

test_dat = basis * x_basis;
noise_10p_gaus = randn(size(dat_test)) * mean(test_dat(:)) * 0.5;
test_dat       = test_dat + noise_10p_gaus;


disp(['mean of dat vs test_dat: ' num2str(mean(dat(:))) ' -> ' num2str(mean(test_dat(:)))]);
disp(['std of dat vs test_dat: ' num2str(std(dat(:))) ' -> ' num2str(std(test_dat(:)))]);
disp(['var of dat vs test_dat: ' num2str(var(dat(:))) ' -> ' num2str(var(test_dat(:)))]);






