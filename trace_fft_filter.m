function varargout = trace_fft_filter(trace_in,filter_range)
% This function applies Fourier-filtering on the input 'trace_in'
% 'filter_range' specifies the range of frequencies for which the content
% should be set to zero. The units for these frequencies is the natural
% frequency associated with the trace.
% (e.g. stepsize of 'trace_in' is in nanoseconds, then the natural frequency is 2*pi*GHz)
% The DC and the highest frequency component is always set to zero.

Ntrace=length(trace_in);
trace_fft=fft(trace_in);
faxis=FourierAxis(1:length(trace_in));
index1(1)=max(vec2ind(faxis(1:floor(Ntrace/2))<=filter_range(1)));
index1(2)=max(vec2ind(faxis(1:floor(Ntrace/2))<=filter_range(2)));
figure;plot(faxis,abs(trace_fft),'k')
trace_fft_filtered=trace_fft;
trace_fft_filtered(index1(1):index1(2))=0;
trace_fft_filtered(Ntrace-index1(2)+2:Ntrace-index1(1)+2)=0;
trace_fft_filtered(1)=0;
trace_fft_filtered(floor(Ntrace/2)+1)=0;
trace_out=ifft(trace_fft_filtered);

figure;hold on;
plot(trace_in,'k')
plot(abs(trace_out),'r')

varargout{1}=trace_out;
varargout{2}=trace_fft_filtered;
varargout{3}=faxis;
end
