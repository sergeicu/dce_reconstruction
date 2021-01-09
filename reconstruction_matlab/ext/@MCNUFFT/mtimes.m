function ress = mtimes(a,bb)


%% b1 is the coil profile. Its size defines the image size
%% bb is the image to be transformed (from or to) fourier space

 if a.adjoint,
     % Multicoil non-Cartesian k-space to Cartesian image domain
     % nufft for each coil and time point
     res = zeros( a.imSize(1), a.imSize(2), size(bb,3),  size(bb,4) );
     for tt=1:size(bb,4),   %% JCF: per each slice
     for ch=1:size(bb,3), %% JCF: for each chanel
         b = bb(:,:,ch,tt).*a.w(:,:,tt);
         res(:,:,ch,tt) = reshape( nufft_adj(b(:),a.st{tt})/sqrt(prod(a.imSize)), [a.imSize(1),a.imSize(2)] );
     end
     end
     % compensate for undersampling factor
     res=res*size(a.b1,1)*pi/2/size(a.w,2);
     ress = zeros(size(res,1),size(res,2),size(res,4));   
     % coil combination for each time point
     ress = squeeze( sum(res.*conj(a.b1),3)./sum(abs((a.b1)).^2,3) ); 
    %  for tt = 1:size(bb,4),   %% JCF: per each slice
    %     ress(:,:,tt) = squeeze( sum(res(:,:,:,tt).*conj(a.b1),3)./sum(abs((a.b1)).^2,3) ); 
    %  end
 else
    ress = zeros(a.dataSize(1),a.dataSize(2), size(a.b1,3),  size(bb,3) ); %[  N, NSpk, NCha, NSli  ]
     % Cartesian image to multicoil non-Cartesian k-space 
     for tt=1:size(bb,3),   %% JCF: per each slice
     for ch=1:size(a.b1,3), %% JCF: for each chanel
        res=bb(:,:,tt).*a.b1(:,:,ch); %#ok<AGROW>
        ress(:,:,ch,tt) = reshape( nufft(res,a.st{tt})/sqrt(prod(a.imSize))  , [a.dataSize(1),a.dataSize(2)]  ).*a. w(:,:,tt);
     end
     end
 end
