clear;
M2=12;
M=sqrt(M2);
beam_width=0.9/2; % from second moment measurement [mm]
input.w=beam_width/M; % [mm]
input.R=1e6; % [mm]
input.lambda=532e-6; % [mm]
system(1).element='fsp';
system(1).param=640;
system(2).element='lns';
system(2).param=+200;
system(3).element='fsp';
system(3).param=190;
system(4).element='lns';
system(4).param=-100;
system(5).element='fsp';
system(5).param=745;
system(6).element='lns';
system(6).param=+100;
system(7).element='fsp';
system(7).param=505;
system(8).element='lns';
system(8).param=+500;
system(9).element='fsp';
system(9).param=+390;
system(10).element='lns';
system(10).param=+500;
system(11).element='fsp';
system(11).param=255;
system(12).element='rff';
system(12).param=[1,1.5];
system(13).element='fsp';
system(13).param=10;
system(14).element='rfc';
system(14).param=[1.5,1,1000];
system(15).element='fsp';
system(15).param=750;

% system.s={640;+200;100;-100;835;+100;700;+600;390;+600;255;[1,1.5];10;[1.5,1,1000];1000};
% system.param={640;+200;190;-100;745;+100;505;+500;390;+500;255;[1,1.5];10;[1.5,1,1000];750}; % actual setup for pump imaging [mm]
% system.s={640;+200;190;-100;745;+100;700;+600;390;+600;255;[1,1.5];10;[1.5,1,1000];1000}; % ideal setup for pump imaging [mm]

N=200;
ind_fsp=[];
ind_lns=[];
ind_rff=[];
ind_rfc=[];
for ind1=1:length(system)
    if strcmp(system(ind1).element,'fsp')
        ind_fsp=cat(2,ind_fsp,ind1);
    elseif strcmp(system(ind1).element,'lns')
        ind_lns=cat(2,ind_lns,ind1);
    elseif strcmp(system(ind1).element,'rff')
        ind_rff=cat(2,ind_rff,ind1);
    elseif strcmp(system(ind1).element,'rfc')
        ind_rfc=cat(2,ind_rfc,ind1);  
    end
end
zplot=[];
ind3=1;
for ind2=ind_fsp
%     z_start=sum(system.s(ind_fsp(ind_fsp<ind2)));
    z_start=0;
    z_end=system(ind2).param;
    z=z_start:(z_end-z_start)/(N-1):z_end;
    if isempty(ind_fsp(ind_fsp<ind2))
        zplot=cat(2,zplot,z);
    else
        zplot=cat(2,zplot,sum([system(ind_fsp(ind_fsp<ind2)).param])+z);
    end
%     input.element=system.element(1:ind2,:);
    subsystem=system(1:ind2);
    for ind1=1:N
        subsystem(ind2).param=z(ind1);
        output=beam_tracer(input,subsystem);
        R(ind3)=output.R;
        w(ind3)=output.w*M;
        ind3=ind3+1;
    end
end
%%
% figure;
hndl1=plot(zplot,w,'k');
xlabel('z [mm]')
ylabel('M\cdotw [mm]')
zrange=max(zplot)-min(zplot);
wmax=1.3*max(w);
for ind1=ind_lns
    dist=sum([system(ind_fsp(ind_fsp<ind1)).param]);
    temp=1.3*w(zplot==dist);
    line([dist dist],[0 wmax],'color','r');
    if system(ind1).param>0
        line([dist-0.03*zrange dist],[0.95*wmax wmax],'color','r');
        line([dist dist+0.03*zrange],[wmax 0.95*wmax],'color','r');
    elseif system(ind1).param<0
        line([dist-0.03*zrange dist],[1.05*wmax wmax],'color','r');
        line([dist dist+0.03*zrange],[wmax 1.05*wmax],'color','r');
    end
end

for ind1=ind_rff
    dist=sum([system(ind_fsp(ind_fsp<ind1)).param]);
    temp=1.3*w(zplot==dist);
    line([dist dist],[0 temp(1)],'color','b');
    if system(ind1).param(2)/system(ind1).param(1)<=1
        line([dist-0.03*zrange dist],[wmax wmax],'color','b');
    else
        line([dist dist+0.03*zrange],[wmax wmax],'color','b');
    end
end

for ind1=ind_rfc
    dist=sum([system(ind_fsp(ind_fsp<ind1)).param]);
    temp=1.3*w(zplot==dist);
    line([dist dist],[0 wmax],'color','g');
    if system(ind1).param(2)/system(ind1).param(1)<=1
        if system(ind1).param(3)>0
        line([dist-0.03*zrange dist],[1.05*wmax wmax],'color','g');
        elseif system(ind1).param(3)<=0
        line([dist-0.03*zrange dist],[0.95*wmax wmax],'color','g');
        end
	else
        if system(ind1).param(3)>0
        line([dist dist+0.03*zrange],[wmax 0.95*wmax],'color','g');
        elseif system(ind1).param(3)<=0
        line([dist dist+0.03*zrange],[wmax 1.05*wmax],'color','g');
        end
    end
end

%%
% figure;plot(zplot,R,'k')
% xlabel('z [mm]')
% ylabel('R [mm]')
% Rmax=0.1*max(R);
% 
% for ind1=ind_lns
%     dist=sum(cell2mat(system.s(ind_fsp(ind_fsp<ind1))));
%     line([dist dist],[0 Rmax],'color','r');
%     if cell2mat(system.s(ind1))>0
%         line([dist-0.03*zrange dist],[0.95*Rmax Rmax],'color','r');
%         line([dist dist+0.03*zrange],[Rmax 0.95*Rmax],'color','r');
%     elseif cell2mat(system.s(ind1))<0
%         line([dist-0.03*zrange dist],[1.05*Rmax Rmax],'color','r');
%         line([dist dist+0.03*zrange],[Rmax 1.05*Rmax],'color','r');
%     end
% end
% 
% for ind1=ind_rff
%     dist=sum(cell2mat(system.s(ind_fsp(ind_fsp<ind1))));
%     line([dist dist],[0 temp(1)],'color','b');
%     if system.s{ind1}(2)/system.s{ind1}(1)<=1
%         line([dist-0.03*zrange dist],[Rmax Rmax],'color','b');
%     else
%         line([dist dist+0.03*zrange],[Rmax Rmax],'color','b');
%     end
% end
% 
% for ind1=ind_rfc
%     dist=sum(cell2mat(system.s(ind_fsp(ind_fsp<ind1))));
%     line([dist dist],[0 Rmax],'color','g');
%     if system.s{ind1}(2)/system.s{ind1}(1)<=1
%         if system.s{ind1}(3)>0
%         line([dist-0.03*zrange dist],[1.05*Rmax Rmax],'color','g');
%         elseif system.s{ind1}(3)<=0
%         line([dist-0.03*zrange dist],[0.95*Rmax Rmax],'color','g');
%         end
% 	else
%         if system.s{ind1}(3)>0
%         line([dist dist+0.03*zrange],[Rmax 0.95*Rmax],'color','g');
%         elseif system.s{ind1}(3)<=0
%         line([dist dist+0.03*zrange],[Rmax 1.05*Rmax],'color','g');
%         end
%     end
% end
