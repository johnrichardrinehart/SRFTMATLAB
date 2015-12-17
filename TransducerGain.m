function TG = TransducerGain(EqualizerZ,LoadZ)
TG = 4*real(EqualizerZ).*real(LoadZ)./abs(EqualizerZ+LoadZ).^2;