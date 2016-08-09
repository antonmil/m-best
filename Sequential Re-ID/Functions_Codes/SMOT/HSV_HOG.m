function   [Feats,Feats_loc] = HSV_HOG(I,idl0,k)

D_si=size(idl0(k).rect,1);
Feats=zeros(D_si,32+4123);
Feats_loc=zeros(D_si,2);
for j=1:D_si
    xr=max(round(idl0(k).rect(j,1)),1);
    yr=max(round(idl0(k).rect(j,2)),1);
    wr=round(idl0(k).rect(j,3));
    hr=round(idl0(k).rect(j,4));
    
    Feats_loc(j,:)=[idl0(k).rect(j,1)+idl0(k).rect(j,3)/2 idl0(k).rect(j,2)+idl0(k).rect(j,4)/2 ];
    Feats(j,1:32) = hsvHistogram(I(yr:yr+hr-1,xr:xr+wr-1,:));
    Ihog=imresize(I(yr:yr+hr-1,xr:xr+wr-1,:),[150,56]);
    HOGMat = vl_hog(single(Ihog), 8) ;
    Feats(j,33:end)=reshape(HOGMat ,1,numel(HOGMat));
end

end
