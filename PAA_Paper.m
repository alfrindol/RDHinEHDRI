seeds = [792103, 981244, 570819, 443221, 102938];
input = 'rosette_oC92.txt';
numMapsList = [1, 2, 4];
useLZCList = [true, false];
useHorizList = [true, false];
%Test(input, numMaps, useLZC, useHoriz, seeds);
Test(input, 4, true, false, seeds);

