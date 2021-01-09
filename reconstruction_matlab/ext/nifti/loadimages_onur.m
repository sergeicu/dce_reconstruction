clear all
im1=zeros(64,64,50,169);
for ii=7:215
    ii
    if ii<10
    aaa=strcat('e6411s3a3i_00',num2str(ii),'.nii');
    else if ii<100
        aaa=strcat('e6411s3a3i_0',num2str(ii),'.nii');
        else
        aaa=strcat('e6411s3a3i_',num2str(ii),'.nii');
        end
    end
    ax=load_untouch_nii(aaa);
    im1(:,:,:,ii-6)=ax.img;
end