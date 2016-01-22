offsetList = generate_offset_list([900, 500, 90, 170],[1100 600 100 190],24,1);

TC = generateTestFPdata(48,9,0,1,workingdir);
[dictionaryParams, paramList] = setDictionaryParams('sphereD170',3);
[similarity, matchedT1, matchedT2, matchedFAdev, M0fit_grad, bestMatch, match_time] = calcSimilarity(TC,48, dictionaryParams, paramList, savingdir, 'sphereD170','test',9);
