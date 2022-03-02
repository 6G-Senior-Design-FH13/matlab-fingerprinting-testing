ft1=zeros(48,480);
ft2=zeros(48,480);
ft3=zeros(48,480);
ft4=zeros(48,480);
for k=1:4
    for l=1:480
        for j=1:4
            temp=zeros(48, 1);
            for i=1:48
                temp(i,1)=features(i,j,k,l);
            end
            if(j==1)
                ft1(:,l)=temp;
            end
            if(j==2)
                ft2(:,l)=temp;
            end
            if(j==3)
                ft3(:,l)=temp;
            end
            if(j==4)
                ft4(:,l)=temp;
            end
        end
    end
end

