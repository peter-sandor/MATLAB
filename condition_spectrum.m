function spectrum_out = condition_spectrum(spectrum_in)

% This function takes a raw spectrum, subtracts background from it, normalizes it, and
% makes spurious pixel values zero. Then the function 'SpectrumConvert' can
% be used to calculate electric field strenght vs omega.
% spectrum_in(:,1) = wavelengths in [nm]
% spectrum_in(:,2) = spectral intensity (unnormalized)

data = spectrum_in;
N_sp = size(data);
data_min = data(:,2);
data_range0 = max(data(:,2)) - min(data(:,2));
data_range = data_range0/10;

Nbin0 = round(data_range/100);
binwidth0 = data_range/(Nbin0-1);
bin_edges0 = (data_min-binwidth0/2:binwidth0:data_min+data_range+binwidth0/2);
bin_centers = (data_min:binwidth0:data_min+data_range);
[counts0,bin_edges] = histcounts(data(:,2),bin_edges0);
bkg_level0 = bin_centers(counts0 == max(counts0)); % determine bkg level roughly.

data_range1 = 30;
binwidth1 = 1;
Nbin1 = data_range1;
bin_edges1 = (bkg_level0-data_range1/2-binwidth1/2:binwidth1:bkg_level0+data_range1/2+binwidth1/2);
bin_centers1 = (bkg_level0-data_range1/2-binwidth1/2:binwidth1:bkg_level0+data_range1/2);
[counts1,bin_edges2] = histcounts(data(data(:,2)<bkg_level0+data_range1/2,2),bin_edges1);
% bkg_level1 = bin_centers1(counts1 == max(counts1)); % determine bkg level with higher accuracy
bkg_level1 = sum(counts1.*bin_centers1)/sum(counts1);

if 0
    figure;
    subplot(211);hold on;
    plot(data(:,2),'ko');
    title('background scan')
    subplot(212);bar(bin_centers1,counts1)
    xlabel('signal level')
    ylabel('counts')
    title('Histogram of background trace')
end

spectrum = (data(:,2)-bkg_level1)/max(data(:,2)-bkg_level1);
wl = data(:,1);
vec0 = 1:N_sp(1);

thrs = 10*std(data(data(:,2)<2*bkg_level1-data_min,2))/max(data(:,2)-bkg_level1); % determine threshold for zeroing as 5*sigma of background noise amplitude
ind_to_zero = (spectrum<thrs);
ind_comp = zeros([N_sp(1),1]);
ind_comp(min(vec2ind(~ind_to_zero)):max(vec2ind(~ind_to_zero))) = 1;

% while sum(ind_to_zero.*ind_comp)
%     thrs = thrs + 5e-4;
%     ind_to_zero = (spectrum<thrs);
%     ind_comp = zeros([N_sp(1),1]);
%     ind_comp(min(vec2ind(~ind_to_zero)):max(vec2ind(~ind_to_zero))) = 1;
% end
spectrum(ind_to_zero) = 0;
spectrum_out = [map2colvec(wl) map2colvec(spectrum)];

% [mean_wl, sigma_wl, FWHM_wl, ind_FWHM, ind_crop] = peak_props([map2colvec(spectrum_out(:,1)) map2colvec(spectrum_out(:,2))]);
% disp(['mean wavelength = ' num2str(mean_wl) ' nm']);
% disp(['FWHM = ' num2str(FWHM_wl) ' nm']);
end