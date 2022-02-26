count=0;
featsC=feats;
labelsC=labels;
for k=1:1:100
    for i=1:1:size(featsC,2)
        curr=featsC(:,i);
        zcount=0;
        for j=1:1:48
            samp=curr(j,1);
            if(samp==0)
                zcount=zcount+1;
            end
        end
        if(zcount==48)
            count=count+1;
            featsC(:,i)=[];
            labelsC(:,i)=[];
        end
    end
end
fCt=0;
for i=1:1:size(featsC,2)
        curr=featsC(:,i);
        zcount=0;
        for j=1:1:48
            samp=curr(j,1);
            if(samp==0)
                zcount=zcount+1;
            end
        end
        if(zcount==48)
            fCt=fCt+1;
            featsC(:,i)=[];
            labelsC(:,i)=[];
        end
end
    