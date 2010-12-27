
% simple test to check functionality of the group function

randa = rand(3, 5);
for i = 1:10
    data(i).dat = randa; %#ok<SAGROW>
end

figure(1)
subplot 211
imagesc(vertcat(data.dat))
subplot 212
g_res = group(data, 'dat');
imagesc(vertcat(g_res.dat));



