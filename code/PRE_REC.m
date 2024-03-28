function [precision,recall,F1] =  PRE_REC(re_cts,pre_cts,n)
data_index = (1:n)';
re_n_cts = setdiff(data_index,re_cts);
pre_n_cts = setdiff(data_index,pre_cts);
TN = length(intersect(re_n_cts,pre_n_cts));
FN = length(pre_n_cts)-TN;
TP = length(intersect(re_cts,pre_cts));
FP = length(pre_cts)-TP;
precision = TP/(TP+FP);
recall = TP/(TP+FN);
F1=2*(precision*recall)/(precision+recall);
end