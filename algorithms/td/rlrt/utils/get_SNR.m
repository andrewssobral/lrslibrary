function [ rel_err, snr ] = get_SNR( X, X0 )

noise = X - X0;
trueNorm = norm(X0);
rel_err = norm( noise ) / trueNorm;
snr = -20*log10( rel_err );

end