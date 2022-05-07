function [a, b] = make_logistic(w0, w1)
    a = log((w1-w0*w1)/(w0-w0*w1));
    b = log((1-w0)/w0);
end