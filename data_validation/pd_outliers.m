
% pd outlier detection


%% load
load('~/Documents/uni/yifat_lab/results/data/all_data_chalva.mat');

hand = 1;
ses2take = find([sessions.hands] == hand);

for i = ses2take
    data(i,:) = sessions(i).pd(hand,:); %#ok<SAGROW>
end
data(data > pi)  = data(data > pi) - 2* pi;


% find sessions which have in more than 2 channels a pd that is more then a
% std away from the mean

all_mean = repmat(circ_mean(data), size(data,1), 1);
all_std  = repmat(circ_std(data) * 2, size(data,1), 1);
out      = sum(abs(data - all_mean) > all_std, 2);


% look at it
figure(1)
for i = 1:12
    subplot(4,4,i);
    [x y] = pol2cart(data(:,i), ones(size(data(:,i))));
    scatter(x, y, 50, out, '.');
    title(circ_std(data(:,i)));
    axis([-1 1 -1 1])
end
subplot(4,4, [13 14])
hist(out, 0:16)

data = data(out < 4, :);
out = out(out < 4);
out(1) = 8;

figure(2)
for i = 1:12
    subplot(4,4,i);
    [x y] = pol2cart(data(:,i), ones(size(data(:,i))));

    scatter(x, y, 50, out, '.');
    title(circ_std(data(:,i)));
    axis([-1 1 -1 1])
end

