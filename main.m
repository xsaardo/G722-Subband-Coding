i = 1;
for ber = [0 1e-4 1e-3 1e-2 1e-1]
    subplot(5,1,i);
    [x,y] = g722('woman1_wb.wav',1,ber);
    plot(y);
    xlim([1 16000]);
    title(['BER = ' num2str(ber)]);
    soundsc(y,16000);
    pause;
    i = i + 1;
end

%%
i = 1;
for mode = [1 2 3]
    subplot(3,1,i);
    [x,y] = g722('woman1_wb.wav',mode,0);
    plot(y);
    xlim([1 16000]);
    title(['Mode = ' num2str(mode)]);
    soundsc(y,16000);
    pause;
    i = i + 1;
end