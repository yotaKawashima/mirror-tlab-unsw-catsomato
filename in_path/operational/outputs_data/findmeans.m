function data_out = findmeans(data)

fprintf('Entered findmeans.\n')

fprintf('\t%s\n', data.custom.filename)

data_out = data;
data_out.trial = mean(data.trial, 3);

data_out.custom.ntrials = 1;
data_out.custom.datatype.ismeantrl = true;
data_out.custom.filename = [data_out.custom.filename '_meantrl'];


fprintf('Exited findmeans.\n\n')