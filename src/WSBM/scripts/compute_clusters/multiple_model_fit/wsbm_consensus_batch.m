
%cd("/scratch/floraliu/worm-connectomics/WSBM");
%cd('/home/witvliet/Dropbox/ZhenLab/worm-connectomics/WSBM');
cd('/n/home03/dwitvliet/WSBM');
run('setup.m')
numThread = str2num(getenv('NUM_THREAD'));
graphIdx = str2num(getenv('DS_IDX'));
k = str2num(getenv('NUM_CLUSTERS'));
currentFolderName = getenv('FOLDER_NAME');

fprintf("matlab prepare to run ds%d k %d using %d cores\n", graphIdx, k, numThread);
repeats = 100; % first set it small
muIter = 100;
mainIter = 100;

outputDir = sprintf("./results/%s/ds%d", currentFolderName, graphIdx);
mkdir(outputDir);
E = load(sprintf('%s_%d.mat', currentFolderName, graphIdx), 'arr');
E = double(E.arr);
E(E == 0) = NaN;
fprintf("results is gonna to save in %s", outputDir);
tic
run_wsbm_large_repeats(E,k, repeats, outputDir, graphIdx, muIter, mainIter, numThread)
toc
