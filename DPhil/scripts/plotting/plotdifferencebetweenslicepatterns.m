coords = [24 22]
phantom_y = squeeze(FPimages(coords(1),coords(2),1,1:10))
phantom_y2 = squeeze(FPimages(coords(1),coords(2),2,1:10))

phantom_y(:) = phantom_y(:)/phantom_y(1)
phantom_y2(:) = phantom_y2(:)/phantom_y2(1)
figure; plot(phantom_y,'*'), hold on, plot(phantom_y2,'o')

simSignal(2,:) = simSignal(2,:)/simSignal(2,1)
plot(1:5, simSignal(2,:),'+')

% the slices have slightly different patterns...


