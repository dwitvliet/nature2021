function [model, label] = load_cns_model(inputdir, ds, k)
    model = load(sprintf("%s/consensus/consensus_model_ds%d_k%d.mat", inputdir, ds, k)); 
    label = load(sprintf("%s/consensus/consensus_label_ds%d_k%d.mat", inputdir, ds, k));    
    try
        model = model.cnsModel;
        label = label.cnsLabel;
        disp('this is a new model');
    catch
        model = model.templateModel;
        label = label.consensus_ca;
        warning('this is an old model');
    end
end