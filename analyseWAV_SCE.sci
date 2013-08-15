
// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);

clc
clear
xdel(winsid())

stacksize('max')

//[y,Fs,bits] = wavread("/opt/work2/1.wav");Fs,bits
//[y,Fs,bits] = wavread("/opt/work2/da.wav");Fs,bits
[y,Fs,bits] = wavread("/opt/work2/r.wav");Fs,bits
//[y,Fs,bits] = wavread("/opt/work2/20130801115515_16000.wav");Fs,bits
//[y,Fs,bits] = wavread("/opt/work2/20130801115553_16000.wav");Fs,bits
//[y,Fs,bits] = wavread("/opt/work2/20130801115415_16000.wav");Fs,bits
//[y,Fs,bits] = wavread("/opt/work2/wavrecorder/20130808152550_16000.wav");Fs,bits // Pattern gut erkennbar
//[y,Fs,bits] = wavread("/opt/work2/wavrecorder/20130808152632_16000.wav");Fs,bits // kein Pattern erkennbar
//[y,Fs,bits] = wavread("/opt/work2/wavrecorder/20130808152656_16000.wav");Fs,bits // Pattern gut erkennbar
//[y,Fs,bits] = wavread("/opt/work2/wavrecorder/20130808152804_16000.wav");Fs,bits // kein Pattern erkennbar
//[y,Fs,bits] = wavread("/opt/work2/wavrecorder/20130812163819_16000.wav");Fs,bits 



PeriodeProBit = 12;
lengthOfFrame = 8;
tf = 2400;
fcutl = 1900;
fcuth = 2900;
lengthOfSignal = size(y,2)
indFcutl = ceil((lengthOfSignal/Fs)*fcutl);
indFcuth = ceil((lengthOfSignal/Fs)*fcuth);

t = 0:1/Fs:(size(y,2)-1)*1/Fs;
subplot(3,1,1)
plot(t,y)
//mtlb_axis([0,6, -2, 2]);

tmuster = mtlb_imp(0,1/mtlb_double(Fs),2/tf);
muster = sin(((2*%pi)*tf)*tmuster);
subplot(3,1,2)

use_svd = 0;
    fLsignal = fft(y);
    fFilter = zeros(max(size(fLsignal)),1);
    fFilter(indFcutl:indFcuth)=1;
    fFilter(lengthOfSignal-indFcuth:lengthOfSignal-indFcutl)=1;
    fsignal = fLsignal.* fFilter';
    Signal = mtlb_ifft(fsignal);
//plot(t,abs(Signal))
b1 = (tf/mtlb_double(Fs))*ones(1,ceil(mtlb_double(Fs)/tf));
Signal2 = filter(b1,1,abs(Signal));
sigma = max(Signal2)/3;

threshold = sigma*ones(length(Signal2),1);
plot(t,[Signal2' threshold]);
//plot(t,threshold);

Signal3 = 0.5*sign(Signal2'-threshold)+0.5;
subplot(3,1,3)
plot(t,abs(Signal3));
FrameOnEdge=sum(Signal3(1:Fs/tf*lengthOfFrame*PeriodeProBit));

if 0~=FrameOnEdge
  ind=Fs/tf*lengthOfFrame*PeriodeProBit+1+...
      find(Signal3(Fs/tf*lengthOfFrame*PeriodeProBit+1:),1);
else
  ind=find(Signal3,1);
end
result=zeros(1,lengthOfFrame);
for i=1:lengthOfFrame
  pos=ind+ceil(PeriodeProBit*Fs/tf/2)+(i-1)*ceil(PeriodeProBit*Fs/tf);
  result(i)=Signal3(pos);
end
result
