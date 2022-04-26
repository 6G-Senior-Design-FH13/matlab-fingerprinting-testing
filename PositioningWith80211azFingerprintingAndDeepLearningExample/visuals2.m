figure()
array=features(:,1,2,1);
plot(array);
hold on
r = randi(480,10,1);
for i=1:10
    array=features(:,1,2,r(i));
    plot(array);
end
hold off
title('CIRs of 10 Random Receivers')
set(gca, 'YScale', 'log')
xlabel('Time Slice "n"') 
ylabel('log(Channel Impulse Response)') 