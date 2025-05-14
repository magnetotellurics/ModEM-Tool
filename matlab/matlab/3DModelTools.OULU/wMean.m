function mean = wMean(values,weights)
% Calculates the weighted mean
% 

W = sum(weights);


mean = weights * values';
if length(mean) > 1 % robust 
    mean = weights' * values;
    if length(mean) > 1
        disp('error')
    end
end

mean = mean/W;