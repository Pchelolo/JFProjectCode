function Rsmth = MakeSmoothWithGaussian( r, sigma ) 
   %Convolves r with Gaussian.
   x = -10000:10000;
   w = exp(-(x/sigma).^2/2)/sqrt(2*pi)/sigma;
   for i = 1 : size(r, 1)
        z = conv(r(i,:), w);
        z=z(ceil(length(w)/2):end-floor(length(w)/2));
        r(i,:) = z;
   end
   Rsmth = r;  
end
