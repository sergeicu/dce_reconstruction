clear all
im2=zeros(128,96,20,100);
for ii=1:100
    ii
    if ii<10
    aaa=strcat('rrarest_00',num2str(ii),'.nii');
    else if ii<100
        aaa=strcat('rrarest_0',num2str(ii),'.nii');
        else
        aaa=strcat('rrarest_',num2str(ii),'.nii');
        end
    end
    ax=load_untouch_nii(aaa);
    im2(:,:,:,ii)=ax.img;
end