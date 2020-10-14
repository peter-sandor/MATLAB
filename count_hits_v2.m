function varargout = count_hits_v2(pic,thrs,Ndivide)

% area_sngl=1100;
% pic2=pic;
% pic2(pic2<=thrs)=0;
% pic2(pic2>thrs)=1;
% nrhits1=round(sum(sum(pic2))/area_sngl);
% Ndivide=ceil(sqrt(nrhits1/5));

% Ndivide=[5 5];  % [y x]
overlap=40;
A=size(pic);
B=floor(A./Ndivide);
nrhits=0;
coords=[];
% dfigure;
fprintf(1,['zone progress (' num2str(Ndivide(1)*Ndivide(2)) ' total): ']);
for ind1=1:Ndivide(1)
    for ind2=1:Ndivide(2)
        zone_x0=(ind2-1)*B(2)+1:ind2*B(2);
        zone_y0=(ind1-1)*B(1)+1:ind1*B(1);
        if ind1<Ndivide(1) && ind1>1 && ind2<Ndivide(2) && ind2>1
            zone_x=(ind2-1)*B(2)+(1-overlap):(ind2*B(2)+overlap);
            zone_y=(ind1-1)*B(1)+(1-overlap):(ind1*B(1)+overlap);
        elseif ind1==Ndivide(1) && ind2<Ndivide(2) && ind2>1
            zone_x=(ind2-1)*B(2)+(1-overlap):(ind2*B(2)+overlap);
            zone_y=(ind1-1)*B(1)+(1-overlap):ind1*B(1);
        elseif ind1<Ndivide(1) && ind1>1 && ind2==Ndivide(2)
            zone_x=(ind2-1)*B(2)+(1-overlap):ind2*B(2);
            zone_y=(ind1-1)*B(1)+(1-overlap):(ind1*B(1)+overlap);
        elseif ind1==1 && ind2<Ndivide(2) && ind2>1
            zone_x=(ind2-1)*B(2)+(1-overlap):(ind2*B(2)+overlap);
            zone_y=(ind1-1)*B(1)+1:(ind1*B(1)+overlap);
        elseif ind1<Ndivide(1) && ind1>1 && ind2==1
            zone_x=(ind2-1)*B(2)+1:(ind2*B(2)+overlap);
            zone_y=(ind1-1)*B(1)+(1-overlap):(ind1*B(1)+overlap);
        elseif ind1==1 && ind2==1
            zone_x=(ind2-1)*B(2)+1:(ind2*B(2)+overlap);
            zone_y=(ind1-1)*B(1)+1:(ind1*B(1)+overlap);
        elseif ind1==1 && ind2==Ndivide(2)
            zone_x=(ind2-1)*B(2)+(1-overlap):ind2*B(2);
            zone_y=(ind1-1)*B(1)+1:(ind1*B(1)+overlap);
        elseif ind1==Ndivide(1) && ind2==1
            zone_x=(ind2-1)*B(2)+1:(ind2*B(2)+overlap);
            zone_y=(ind1-1)*B(1)+(1-overlap):ind1*B(1);
        elseif ind1==Ndivide(1) && ind2==Ndivide(2)
            zone_x=(ind2-1)*B(2)+(1-overlap):ind2*B(2);
            zone_y=(ind1-1)*B(1)+(1-overlap):ind1*B(1);         
        end
        pic2=pic(zone_y,zone_x);
        offset=[min(zone_x)-1 min(zone_y)-1];
        temp2=centroid_v2(pic2,thrs);
        if size(temp2,1)==1
            temp4=temp2+offset;
        elseif size(temp2,1)>1
            temp4=temp2+permute(extend(offset,size(temp2,1)),[2 1]);
        else
            temp4=[];
        end
        index_delete=[];
        for ind3=1:size(temp4,1)
            index_delete(ind3,1)=0;
            if temp4(ind3,1)<min(zone_x0) || temp4(ind3,1)>max(zone_x0) || temp4(ind3,2)<min(zone_y0) || temp4(ind3,2)>max(zone_y0)
                index_delete(ind3,1)=1;
            end
        end
        
        hitimg=plot_FragImgP(temp2,size(pic2));        
        
%         colormap bone
%         colormap(flipud(colormap));
%         imagesc(zone_x,zone_y,hitimg)
%         line(zone_x0,min(zone_y0),'color','r');
%         line(zone_x0,max(zone_y0),'color','r');
%         line(min(zone_x0),zone_y0,'color','r');
%         line(max(zone_x0),zone_y0,'color','r');
        
%         nrhits=nrhits+temp1;
        temp4(logical(index_delete),:)=[];
        coords=cat(1,coords,temp4);
        fprintf(1, '.');
    end
end
fprintf(1, '\n');
% varargout{2}=nrhits;
varargout{1}=coords;
end


