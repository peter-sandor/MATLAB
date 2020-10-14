function out = rotational_distribution(simresult)
h_SI=6.626e-34; % [J*s]
c_SI=2.99792e8; % [m/s]
MaxM=size(simresult.amps0,3);
Jvec=map2colvec(0:simresult.input.maxJ);
E_SI=h_SI*c_SI*100*(simresult.input.rot_const_B*Jvec.*(Jvec+1)+simresult.input.centrif*Jvec.^2.*(Jvec+1).^2);
E_au=map2rowvec(E_SI/(27.211385*1.60218e-19));

Jpop=zeros([length(simresult.delay0) simresult.input.maxJ+1 1]);
Jpop0=zeros([simresult.input.maxJ+1 1]);
Jpop1=zeros([simresult.input.maxJ+1 1]);
amps0=simresult.amps0;
if length(size(amps0))==4
    for indM=1:MaxM;
        if indM>1
            amps1(:,:,indM,:)=sqrt(2)*squeeze(amps0(:,:,indM,:)).*permute(extend(extend(sqrt(simresult.probs_JM(1:MaxM)),length(simresult.delay0)),size(amps0,2)),[2 3 1]);
        else
            amps1(:,:,indM,:)=squeeze(amps0(:,:,indM,:)).*permute(extend(extend(sqrt(simresult.probs_JM(1:MaxM)),length(simresult.delay0)),size(amps0,2)),[2 3 1]);
        end
        ReAmp1(indM,:)=real(squeeze(sum(amps1(end,:,indM,:),4)));
        phase1(indM,:)=unwrap(angle(squeeze(sum(amps1(end,:,indM,:),4))));
        Jpop=Jpop+map2colvec(squeeze(abs(sum(amps1(:,:,indM,:),4))).^2);
        Jpop0=Jpop0+map2colvec(abs(sum(amps1(1,:,indM,:),4)).^2);
        Jpop1=Jpop1+map2colvec(abs(sum(amps1(end,:,indM,:),4)).^2);
    end
elseif length(size(amps0))==5
    MaxK=size(simresult.amps0,4);
    for indK=1:MaxK
        for indM=1:MaxM;
            amps1(:,:,indM,indK,:)=squeeze(amps0(:,:,indM,indK,:)).*permute(extend(extend(sqrt(squeeze(simresult.probs_JMK(1:size(amps0,5),indK))),length(simresult.delay0)),size(amps0,2)),[2 3 1]);
            ReAmp1(indM,indK,:)=real(squeeze(sum(amps1(end,:,indM,indK,:),5)));
            phase1(indM,indK,:)=unwrap(angle(squeeze(sum(amps1(end,:,indM,indK,:),5))));
            Jpop1=Jpop1+map2colvec(abs(squeeze(sum(amps1(end,:,indM,indK,:),5))).^2);
            Jpop0=Jpop0+map2colvec(abs(sum(amps1(1,:,indM,indK,:),5)).^2);
            Jpop=Jpop+map2colvec(abs(sum(amps1(:,:,indM,indK,:),5)).^2);
        end
    end
end

out=[map2colvec(Jvec) map2colvec(Jpop0) map2colvec(Jpop1)];

figure;
hold on;
hax=bar(Jvec,[map2colvec(Jpop0) map2colvec(Jpop1)]);
% bar(Jvec,Jpop1,'r')
xlim([-0.5 min(vec2ind(Jpop1<1e-4))])
xlabel('J')
legend('before pulse','after pulse');
ylabel('state population');
set(hax(1),'facecolor','k','edgecolor','k')
set(hax(2),'facecolor','r','edgecolor','r')
setfigP;
end