    num=input('Choose the signal number from 1 to 3\n');
    if(num==1)
     load('mysteryA.mat');
              SRRCLength=     5;
              SRRCrolloff=    0.25;
              T_t=            9e-6;
              f_if=           1.88e6;
              f_s=            820e3;
              fileID='Mystery A.txt';
    elseif(num==2)
     load('mysteryB.mat');
              SRRCLength=     3;
              SRRCrolloff=    0.35;
              T_t=          7.3e-6;
              f_if=           1.92e6;
              f_s=            680e3;
              fileID='Mystery B.txt';
    elseif(num==3)
     load('mysteryC.mat');
              SRRCLength =    4;
              SRRCrolloff =   0.3;
              T_t   =        8.2e-6;
              f_if =          2.4e6;
              f_s=            760e3;
              fileID='Mystery C.txt';
    end
    
    r=r';
    k= round(f_if/f_s);
    fc=abs(f_if-k*f_s);                     
    m= f_s*T_t;
    Ts=1/f_s;

    %% filtering 
    lcut = 0.8; 
    b=firpm(500,[0 lcut lcut+0.01 1],[1 1 0 0]);
    fr=filtfilt(b,1,r);
    figure, plotspec(fr,1/f_s);
    title('Filtered signal');

    %% Carrier Recovery
    CarrierEst=dualCostasLoop(fr,f_s,fc);

    %% Demodulation
    r1=fr.*CarrierEst;                      % demod received signal

    figure, plotspec(r1,Ts);
    title('Demodulated signal');

    %% Matched Filtering

    recfilt=srrc(SRRCLength,SRRCrolloff,m);               % receive filter H sub R
    r2=filter(fliplr(recfilt),1,r1);           % matched filter with data
    figure, plotspec(r2, Ts);
    title('Matched filter output');

   %% Time Recovery

    l= SRRCLength;  n= round(length(r2)/m);
    tnow=l*m+1; tau=0; tau1=0; xs=zeros(1,n); xs1=zeros(1,n); r3=zeros(1,n);         % initialize variables
    tausave=zeros(1,n); tausave(1)=tau; i=0;
    tausave1=zeros(1,n); tausave1(1)=tau1; 
    mu=0.1; mu1=0.0035;                                   % algorithm stepsize
    delta=0.1; delta1=0.07;                               % time for derivative

    while tnow<length(r2)-2*l*m                  % run iteration
     i=i+1;
     xs(i)=interpsinc(r2,tnow+tau,l);           % interpolated value at tnow+tau
     x_deltap=interpsinc(r2,tnow+tau+delta,l);  % get value to the right
     x_deltam=interpsinc(r2,tnow+tau-delta,l);  % get value to the left
     dx=x_deltap-x_deltam;                      % calculate numerical derivative
     qx=quantalph(xs(i),[-3,-1,1,3]);           % quantize xs to nearest 4-PAM symbol
     tau=tau+mu*dx*(qx-xs(i));                  % alg update: DD
     tausave(i)=tau;                            % save for plotting

     xs1(i)=interpsinc(r2,tnow+tau+tau1,l);             % interpolated value at tnow+tau
     x_deltap1=interpsinc(r2,tnow+tau+tau1+delta1,l);   % get value to the right
     x_deltam1=interpsinc(r2,tnow+tau+tau1-delta1,l);   % get value to the left
     dx1=x_deltap1-x_deltam1;                           % calculate numerical derivative
     qx1=quantalph(xs1(i),[-3,-1,1,3]);                 % quantize xs to nearest 4-PAM symbol
     tau1=tau1+mu1*dx1*(qx1-xs1(i));                    % alg update: DD
     tausave1(i)=tau1;                                  % save for plotting
     r3(i)=interpsinc(r2,tnow+tau+tau1,l);
     tnow=tnow+m;
     end
    r3= 1.2.*r3(1:i-2);

figure,                                                 % plot results
    subplot(3,1,1), plot(r3(1:i-2),'b.')                % plot constellation diagram
    title('Constellation diagram');
    ylabel('Estimated symbol values');
    subplot(3,1,2), plot(tausave(1:i-2))                % plot trajectory of tau
    title('Time Recovery 1st loop'); 
    ylabel('Offset estimates'), xlabel('Iterations');
    subplot(3,1,3), plot(tausave1(1:i-2))               % plot trajectory of tau
    title('Time Recovery 2nd loop');
    ylabel('Offset estimates'), xlabel('Iterations');

%% Correlation

head = letters2pam('A0Oh well whatever Nevermind');
c=xcorr(head, r3);                                      % do cross correlation
[m,ind]=max(abs(c));                                    % location of largest correlation
headstart=length(r3)-ind+1;                             % place where header starts
headstart = mod(headstart, (112+400));
frameLength= 512;

r4 =r3(headstart+512-8:end);                            % we discared the first frame of the signal but we take the last 2 letters from it "8"
start = 1; ende=512;

figure,
subplot(3,1,1), stem(head)                       % plot header
title('Header')
subplot(3,1,2), stem(r3)                         % plot data sequence
title('Data with embedded header')
subplot(3,1,3), stem(c)                          % plot correlation
title('Correlation of header with data');


ref = mod(length(r4), frameLength)-1;                           
r5 = r4(1:length(r4)-ref);                       %removing extra noise

num3 = mod(length(r4), frameLength);
r4 = r4(1:length(r4)-num3);
num4 = length(r4);
num5 = mod (num4, frameLength);
noFrames = num4/frameLength;


while ende<=length(r5)
    pos = start:1:ende;
    frame = r5(start:ende);                  % assigning new frame of header+data
    c=xcorr(head, frame);                    % do cross correlation
    [m,ind]=max(c);                          % location of largest correlation
    [m1,ind1]=min(c);                        % location of largest correlation
    if (abs(m1) > m )
        r5(pos) = -r5(pos);
    end
    start = start + 512;                     % iterate               
    ende = ende + 512;                       % iterate
end

r6 = zeros(1,noFrames*400);
rh = zeros(1,noFrames*112);
c=xcorr(head, r5);                     % do cross correlation
figure,
subplot(3,1,1), stem(head)             % plot header
title('Header')
subplot(3,1,2), stem(r5)               % plot data sequence 
title('Data with embedded header')
subplot(3,1,3), stem(c)                % plot correlation
title('Correlation of header with data')

%% Equalizer

headPos = 1; frameEnd = 512;
no=9; f=[0 1 0 0 0 0 0 0 0]';          % initialize equalizer
mu=.008;                               % stepsize
index = 1;
eqPos =1; eqEnd=400;
hstart=1;

for i=no+1:length(r5)                        % iterate
  rr=r5(i:-1:i-no+1)';                       % vector of received signal
  e=quantalph(f'*rr,[-3 -1 1 3])-f'*rr;      % calculate error
  f=f+mu*e*rr;                               % update equalizer coefficients
  r5(index) = rr'*f;
  index = index +1;                          %iterate   
end

figure, plot([1:length(r5)],r5,'.')
title('EQ Constellation diagram');

%% Frame sync
while frameEnd <= length(r5) 
        locData = eqPos:1:eqEnd;
        frame = r5(headPos:frameEnd);           % assigning new frame of header+data
        r6(locData) = frame(113:512);           % saving the data part in r6 
        headPos = headPos +512;                 % iterate
        frameEnd = frameEnd + 512;              % iterate
        eqPos =eqPos + 400;                     % iterate
        eqEnd=eqEnd+400;                        % iterate

end

%% Decoding

mprime=quantalph(r6,[-3,-1,1,3])';              % quantize to +/-1 and +/-3 alphabet

% decode decision device output to text string
reconstructed_message=pam2letters(mprime);      % reconstruct message
disp(reconstructed_message);
