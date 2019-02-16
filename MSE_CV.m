function [ mse_cv ] = MSE_CV( diff_cv )

m1=length(diff_cv);
mse_cv=(1/(2*m1))*sum(diff_cv.^2);



end

