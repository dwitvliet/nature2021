%cd(getenv('workdir'));
run('setup.m');
folderName=getenv('FOLDER_NAME');
graphIdx=str2num(getenv('DS_IDX'));
k=str2num(getenv('NUM_CLUSTERS'));
numThread = str2num(getenv('NUM_THREAD'));

folderName='B-all';
graphIdx=0;
k=6;
numThread = 16;

resultDir = sprintf("results/%s", folderName);

for graphIdx=1:3
    for k=6:8
        
       if graphIdx == 1 && k ~= 6
           continue
       end
       if graphIdx == 3 && k ~= 6
           continue
       end
       
       disp([graphIdx k])

tic
fitConsensusModel(resultDir, graphIdx, k, numThread);
fprintf("\n fitted the consensus model for ds %d num clusters %d\n", graphIdx, k);
toc

    end
end