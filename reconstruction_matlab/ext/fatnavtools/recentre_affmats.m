function newAffMats = recentre_affmats(affMats,iCentreVol)

nV = size(affMats,3);
newAffMats = zeros(size(affMats));
if length(iCentreVol) == 1
    affMatRef = affMats(:,:,iCentreVol);
else
    affMatRef = mean(affMats(:,:,iCentreVol),3);
end

for iV = 1:nV
    newAffMats(:,:,iV) = affMats(:,:,iV) / affMatRef;    
end