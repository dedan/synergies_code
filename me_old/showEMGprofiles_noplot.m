function e2take = showEMGprofiles_noplot( emg )

pre = -500;
post = 1000;

if isempty(emg),
    return;
end
Lng = length(emg(1).hand(1).target(1,:));
dt = (post-pre)/Lng;
ax = pre:dt:post;
ax=ax(1:Lng);
  
zeroindx = max(find(ax<= 0));
baseline = find(ax>= pre & ax <= -300);
for i=1:length(emg),
    Nh = length(emg(i).hand);
    for k=1:Nh
        mresp = mean(emg(i).hand(k).target(:,zeroindx:end)');
        mbase = (mean(emg(i).hand(k).target(:,baseline)'));
        sbase = (std(emg(i).hand(k).target(:,baseline)'));
        tmp  = mresp-(mbase+2*sbase);
        if ~isempty(find(tmp> 0)),
            e2take(i,k) = 1;
        else
            e2take(i,k) = 0;
        end
     end;
end
    