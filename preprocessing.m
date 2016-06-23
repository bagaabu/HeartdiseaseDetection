
%% put everything into the wav cell
clear wav;
wav(1,:) = {'index','location','name','raw','normalised'};
d = dir(PATH);
for i = 3:size(d,1)
   if(d(i).isdir ~= 0)
        wavset = wavSet([PATH,'/',d(i).name]);
        for j = 1:wavset.Count
            wave = audioread(wavset.WavLocation{j});           
            wave = decimate(wave,decimatedRate);
            ori = getSignal(wave,Fs);
            % normalise the original signal
            wave = (wave-min(wave))/(max(wave)-min(wave));
            wave = wave - mean(wave);
            ww = getSignal(wave,Fs);
            wav(size(wav,1)+1,:) = {i-3,wavset.WavLocation{j},d(i).name,ori,ww};
        end
   end
end
string = {'normal','artifact','murmur','extrahls','normal_noisy','murmur_noisy'};
for i  = 1:size(string,2)
    temp = find(strcmp(string{i},wav(:,3)) == 1);
    if (temp ~= 0)
        eval([string{i},' = ','[',num2str(temp(1)),',',num2str(temp(length(temp))),']']);
    end
end
all = [2,size(wav,1)];
%% saving rawFFT into wav cell
% NO.5 column is the FFT data, 2X2 cell include x axis & y axis
wav(1,RAW+1) = {'rawFFT'};
for i = 1:size(wav,1)-1
    w = wav{i+1,RAW};
    FFT = getFFT(w{2,2},Fs);
    wav(i+1,RAW+1) = {FFT};
end
%% Using wavelet to denoise signal (DEN), and store in wav col 6&7
% wav(1,DEN) = {'waveletDEN'};
% wav(1,DEN+1) = {'DENFFT'};
% 
% for i = 1:size(wav,1)-1
%     Fs = wav{i+1,3};
%     w = wav{i+1,RAW};
%     % [c,l] = wavedec(w,5,'db10');
%     sigDEN = wden(w,'minimaxi','s','mln',level,'db10');
%     wav(i+1,DEN) = {sigDEN};
%     DFFT = getFFT(sigDEN,Fs);
%     wav(i+1,DEN+1) = {DFFT};
% end
%% Using wavelet decomposition to denoise signal (DEN)
wav(1,DEN2) = {'waveletDEN2'};
wav(1,DEN2+1) = {'DEN2FFT'};
for i = 1:size(wav,1)-1
    w = wav{i+1,RAW}{2,2};
    sigDEN = wden(w,'minimaxi','s','mln',level,'db10');
    thr = thselect(w,'minimaxi');
    thrSettings =  [...
        thr ; ...
        thr ; ...
        0 ; ...
        0 ; ...
        0 ; ...
        ];
    sigDEN2 = cmddenoise(sigDEN,'db10',level,'s',NaN,thrSettings);
    sigDEN2 = (sigDEN2-min(sigDEN2))/(max(sigDEN2)-min(sigDEN2));
    sigDEN2 = sigDEN2 - mean(sigDEN2);
    s = getSignal(sigDEN2',Fs);
    wav(i+1,DEN2) = {s};
    DFFT2 = getFFT(s{2,2},Fs);
%     DFFT2{2,2} = envelope(DFFT2{2,2},10,'peak');
    wav(i+1,DEN2+1) = {DFFT2};
end
%% Shannon energy 1, Liang's method
wav(1,ENERGY) = {'Shannon energy 1'};
for i = 1:size(wav,1)-1
    w = wav{i+1,DEN2}{2,2};
    % normalise first time
    w = w./max(abs(w));
    % average shannon energy method
    e = zeros(0);
    for j = 22:22:size(w,1)-22
        t1 = w(j-21:j+22,1);
        t2 = -1/44*sum(t1.^2.*log10(t1.^2));
        e(length(e)+1,1) = t2;
    end
    % normalise again
    e = (e-mean(e))/std(e);
    e = (e-min(e))/(max(e)-min(e)); 
    
    T = getSignal(e,Fs/22);
    wav(i+1,ENERGY) = {T};
end
%% Shanoon energy 2?Our methods
wav(1,ENERGY+1) = {'Shannon energy 2'};
for i = 1:size(wav,1)-1
    w = wav{i+1,DEN2}{2,2};
    
    % normalise
    w = w./max(abs(w));
    % shanoon energy
    e = -(w.^2).*log10(w.^2);
    
    e = (e - mean(e))/std(e);
    e = (e-min(e))/(max(e)-min(e));
    
    % envelope by 44
    e = envelope(e,44,'peak');
    T = getSignal(e,Fs);
    wav(i+1,ENERGY+1) = {T};
end
%% Get the cardiac cycle time
wav(1,CYCLE) = {'Cycle'};
for i = 1:size(wav,1)-1
    t = wav{i+1,ENERGY}{2,1};
    w = wav{i+1,ENERGY}{2,2};
    c = getCycle(t,w);
    wav(i+1,CYCLE) = {c};
end
%% Get S1 S2 location
wav(1,S1S2) = {'S1S2 locs'};
for i = 1:size(wav,1)-1
    t = wav{i+1,ENERGY}{2,1};
    w = wav{i+1,ENERGY}{2,2};
    interval = wav{i+1,CYCLE}{2,3};
    interval50 = wav{i+1,CYCLE}{2,4};
    t_st = wav{i+1,CYCLE}{2,1}(1);
    s = getS1S2(t,w,interval,interval50,t_st);
    wav(i+1,S1S2) = {s};
end
%% Get Sysole and Diastole cycle time
wav(1,CYCLE2) = {'2 cycle times'};
clear h1;
clear label;
for i = 1:size(wav,1)-1
    t = wav{i+1,S1S2}{2,1};
    d = diff(cell2mat(t));
    [h,~] = hist(d,linspace(0.05,1,20));
    h1(i,:) = [h,wav{i+1,CYCLE}{2,2},wav{i+1,CYCLE}{2,2}];
    label(i,1) = wav{i+1,1};
end