function CarrierEst=dualCostasLoop(r,fs,fc)

t = 0:(1/fs):((length(r)/fs)-(1/fs));                  % Initialize Time Vector
fl=50; ff=[0 .01 .02 1]; fa=[1 1 0 0];
h=firpm(fl,ff,fa);                        % LPF design
fl2=200;
h2=firpm(fl2,ff,fa);
mu=.13;                                  % algorithm stepsize
mu2=0.00018;
CarrierEst=zeros(1,length(t));                             % assumed freq. at receiver
theta=zeros(1,length(t)); theta(1)=0;     % initialize estimate vector
theta2=zeros(1,length(t)); theta2(1)=0;
zs=zeros(1,fl+1); zc=zeros(1,fl+1);       % initialize buffers for LPFs
zs2=zeros(1,fl2+1); zc2=zeros(1,fl2+1);
for k=1:length(t)-1                       % z's contain past fl+1 inputs
  zs=[zs(2:fl+1), 2*r(k)*sin(2*pi*fc*t(k)+theta(k))];
  zc=[zc(2:fl+1), 2*r(k)*cos(2*pi*fc*t(k)+theta(k))];
  lpfs=fliplr(h)*zs'; lpfc=fliplr(h)*zc'; % new output of filters
  theta(k+1)=theta(k)-mu*lpfs*lpfc;       % algorithm update

  zs2=[zs2(2:fl2+1), 2*r(k)*sin(2*pi*fc*t(k)+theta2(k)+ theta(k))];
  zc2=[zc2(2:fl2+1), 2*r(k)*cos(2*pi*fc*t(k)+theta2(k)+ theta(k))];
  lpfs2=fliplr(h2)*zs2'; lpfc2=fliplr(h2)*zc2'; % new output of filters
  theta2(k+1)=theta2(k)-mu2*lpfs2*lpfc2;       % algorithm update
  CarrierEst(k)=cos(2*pi*fc*t(k)+theta2(k)+theta(k));
end
figure,
subplot(2,1,1), plot(t,theta)              % plot first theta
title('Output of first CL')
subplot(2,1,2), plot(t,theta2)              % plot second theta
title('Output of second CL')


