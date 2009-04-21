function [gridsmooth]=smoothwinds_1(wind_temp_nam,wind_temp_gfs)
gridsmooth=wind_temp_nam;
maxwind=max(max(wind_temp_gfs))+10;%+std(std(wind_temp_gfs));
[eta xi]=size(gridsmooth);
smcells=100; %number or length of cells to smooth over
inc=((100/smcells)/100);%percentage to smooth by

    
for thisrow=1:eta
    [cfstx cfsty] = find(gridsmooth(thisrow,:)<maxwind,1,'first');%logic matrix, what to do...
    [cendx cendy] = find(gridsmooth(thisrow,:)<maxwind,1,'last');
    %replace nam with gfs for data outside interp zone
        if cfsty>1 && cendy<xi
            gridsmooth(thisrow,1:cfsty+4)=wind_temp_gfs(thisrow,1:cfsty+4);
            gridsmooth(thisrow,cendy-4:end)=wind_temp_gfs(thisrow,cendy-4:end);
        elseif cfsty==1 && cendy<xi
            gridsmooth(thisrow,cendy-4:end)=wind_temp_gfs(thisrow,cendy-4:end);
        end
    %create a loop to smooth first and last five rows and columns using a
    if cfsty>1
        cfsty=cfsty+4;
    end
    if cendy<xi
        cendy=cendy-4;
    end
    if (cendy-cfsty)>=(smcells*2)
        for count=0:(smcells-1)
            gridsmooth(thisrow,cfsty+count)=((wind_temp_nam(thisrow,cfsty+count)*((count+1)*inc))+(wind_temp_gfs(thisrow,cfsty+count)*(1-((count+1)*inc))));
            gridsmooth(thisrow,cendy-count)=((wind_temp_nam(thisrow,cendy-count)*((count+1)*inc))+(wind_temp_gfs(thisrow,cendy-count)*(1-((count+1)*inc))));
        end
    elseif (cendy-cfsty)<(smcells*2)
        dy=(cendy-cfsty);
        inc2=((100/(dy/2))/100);
        for count=0:((dy/2)-1)
            gridsmooth(thisrow,cfsty+count)=((wind_temp_nam(thisrow,cfsty+count)*((count+1)*inc2))+(wind_temp_gfs(thisrow,cfsty+count)*(1-((count+1)*inc2))));
            gridsmooth(thisrow,cendy-count)=((wind_temp_nam(thisrow,cendy-count)*((count+1)*inc2))+(wind_temp_gfs(thisrow,cendy-count)*(1-((count+1)*inc2))));
        end
    end
end
    cleft = find(gridsmooth>maxwind);%logic matrix, what to do...
    gridsmooth(cleft)=wind_temp_gfs(cleft);

