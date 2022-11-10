function h=plotEpoch(data,chNum,label,fs,StartEndT,Chan_Name,fileName)

n=35;
nRow=5;
figNum=ceil(size(data,3)/n);
for i =1:figNum
    
    h(i)=figure;
    for j=1:n
        ind=(i-1)*n+j;
        if ind <= size(data,3)
        subplot(nRow,n/nRow,j);
       
        t=1:length(data);
        t=t/fs;
        plot(t,data(:,chNum,ind));
        xlim([0,5])
        ylim([-4,4])
        
        a=[string(ind),string(label(ind))];
        a=join(a);
        title(a);
        end
        
        
    end
    sgtitle(join(['ch:',string(Chan_Name(chNum)),fileName])) 
    set(gcf,'Position',[100,100,1300,900]);
end




end