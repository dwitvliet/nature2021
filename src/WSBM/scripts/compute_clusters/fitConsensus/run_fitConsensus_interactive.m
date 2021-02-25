function run_fitConsensus_interactive()
  dirNames = ["C-no-noise"  "D-no-variable"  "E-stable"  "F-common-exclude1"  "G-common"];
  DS = 0:7
  ks = [2     2     2     2     4     4     5     4; %C
        2     2     2     2     2     3     5     4; %D
        2     2     2     2     2     2     2     2; %E
        2     2     2     2     2     2     2     2; %F
        2     2     2     2     2     2     2     2];%G
  for d = 1:numel(dirNames)
    resultDir=sprintf('./results/%s', dirNames(d))
    ks(d,:)
    for i = 1:numel(DS)
      ds = DS(i)
      k = ks(d,i)
      fitConsensusModel(resultDir, ds, k);
    end
  end
end
