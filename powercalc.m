%Calculates rms voltage, power, and converts if using an attenuator

voltage=input('What was the peak-to-peak voltage? ');
attenuate=input('How much attenuation do you have(dB)? ');
rmsvoltage=voltage/(2*sqrt(2))
power=(rmsvoltage)^2/50
totalpowerout=power*10^(attenuate/10)
dbm=10*log10(totalpowerout*1000)
totalvoltage=rmsvoltage*10^(attenuate/20)



