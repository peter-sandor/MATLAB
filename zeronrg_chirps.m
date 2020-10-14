function varargout = zeronrg_chirps(chirps,peaks,pulsedur)
% function for evaluating the energies of various peaks in the
% photoelectron spectrum at zero field intensity. Works on the assumption
% that the peak positions shift linearly with intensity.
% input:
% 'chirps'      - cell array containing the chirps set with the pulse shaper, which correspond to the entries in the cell array 'peaks'
% 'peaks'       - cell array, each cell contains an array of peak position
%               energies for a specific peak, each array element within a given cell
%               corresponds to a specific AOM amplitude contained in the corresponding entry of the cell array amps
% 'pulsedur'    - cell array with the same size as 'chirps', containing the
%               values for pulse duration in [fs]
% output:
% 'nrgout'      - determined peak postitions at zero intensity

M=size(peaks,2);
COLORS=colormap(lines(M));
figure;
xlabel('Intensity [normalized units]','Fontsize',14)
ylabel('Photoelectron peak positions [meV]','Fontsize',14)
hold on
if iscell(chirps)==1
    intensities=cellfun(@chirp2int,pulsedur,chirps,'UniformOutput',0);
    xmax=max(cell2mat(cellfun(@max,cellfun(@max,intensities,'UniformOutput',0),'UniformOutput',0)));
    for ind1=1:M
        handle(ind1)=plot(intensities{ind1},peaks{ind1},'Color',COLORS(ind1,:),'Marker','+','Linewidth',2);
        FIT=ezfit(intensities{ind1},peaks{ind1},'poly1');
        plot((0:intensities{ind1}(1)/9:intensities{ind1}(1)),FIT.m(1)+FIT.m(2)*((0:intensities{ind1}(1)/9:intensities{ind1}(1))),'r--')
        varargout{1}(ind1)=FIT.m(1);
    end
elseif isnumeric(chirps)==1
    intensities=chirp2int(pulsedur,chirps);
    xmax=max(max(intensities));
    for ind1=1:M
        handle(ind1)=plot(intensities,peaks(:,ind1),'Color',COLORS(ind1,:),'Marker','+','Linewidth',2);
        FIT=ezfit(intensities,peaks(:,ind1),'poly1');
        plot((0:intensities(1)/9:intensities(1)),FIT.m(1)+FIT.m(2)*((0:intensities(1)/9:intensities(1))),'r--')
        varargout{1}(ind1)=FIT.m(1);
    end
    % ylim([0,max(max(peaks)*1.1)]);
end
xlim([0 xmax]);
hold off
varargout{2}=handle;
end