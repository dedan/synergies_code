

test_patterns = [randn(20,2) * 10; randn(20,2) * 2;10+ randn(20,2) * 2;-15+ randn(20,2) * 3; -45 + randn(20,2) * 5;60+ randn(20,2) * 2;];

%
hold off;
figure(1);
plot(test_patterns(:,1), test_patterns(:,2),'*');
hold on;
prots = ng_protos(test_patterns, 5, 1,10, 100, 0.999, 0.9);
plot(prots(:,1),prots(:,2),'*g');


