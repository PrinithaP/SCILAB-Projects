
data=128 ; 
M = 2;
nblock= 1;
b = round(rand(1, data*log(M)/log(2)*nblock));//bits
fs= 20e6; //sampling frequency
rs= fs/data; //symbol rate
rb= rs*log(M)/log(2); //bit rate
tb= 1/rb //bit duration
ts= 1/fs //sampling rate
nb= 1 //bit index    
subband= 8;
Nfft= 128;
subcarriers= 16;
k_data= [[-64 :-1],[0:63]];
k_data(k_data<0)= k_data(k_data<0)+Nfft;
k_data= fftshift(k_data) 

k_center= [-54 -38 -22 -6 10 26 42 58];
  
k_center(k_center<0) = k_center(k_center<0)+Nfft;
k_center=  fftshift(k_center);

taps= round(3/4*Nfft); 
L= Nfft+taps-1;  // filter length
nor_cutoff= 10/128;  
n= 0:taps-1;
[h_c,H_k,fre]=wfir('lp',taps,nor_cutoff,'hm',[0,0])
Co=1.5;
[h_c,H_k,fre]=wfir('lp',taps,nor_cutoff,'hm',[0,0])
h_n=h_c.*exp(2i*%pi*(-Co)/Nfft*n) // center moves to -1.5
H_K=fft([h_n zeros(1,Nfft-taps) ])

figure
subplot(3,1,1);title('IMPULSE RESPONSE OF PROTOTYPE LOW PASS FILTER(h_c)');plot(h_c);
xlabel('Sample Index (n)');ylabel('Amplitude')
subplot(3,1,2);title('FREQUENCY RESPONSE OF PROTOTYPE LOW PASS FILTER (H_k)');plot(fre*Nfft,abs((H_k)))
xlabel('Frequency Bin (k)');ylabel('|H(k)|')
subplot(3,1,3);title('PHASE RESPONSE OF PROTOTYPE LOW PASS FILTER (H_k)');plot(unwrap(angle(H_k)));
xlabel('Frequency Bin (k)');ylabel('Phase ')
figure
subplot(4,1,1);title('IMPULSE RESPONSE OF FREQUENCY SHIFTED PROTOTYPE FILTER(h_n)');plot(h_n);xlabel('Sample Index (n)');ylabel('Amplitude')
subplot(4,1,2);title('FREQUENCY RESPONSE OF FREQUENCY SHIFTED PROTOTYPE FILTER (H_K)');plot(abs(H_K))
xlabel('Frequency Bin (k)');ylabel('|H(k)|')
subplot(4,1,3);title('PHASE  RESPONSE OF FREQUENCY SHIFTED PROTOTYPE FILTER(H_K)');plot(unwrap(angle(H_K))) 
xlabel('Frequency Bin (k)');ylabel('Phase ')

l=0:subband-1;
pilot_intervalLp=16
k_data_pilot=l*pilot_intervalLp
pilindex=zeros(1,subband)
pilindex([k_data_pilot]+1)=ones(1,subband);
 

x_sum=zeros(1,Nfft+taps-1) //filter output    
s_sum=zeros(1,Nfft) //ifft 
sub_S=zeros(1,Nfft);
X_sum=zeros(1,256) //filter output 
for ii=1:nblock //block
    for j=1:8
        S=zeros(1,Nfft) //subband
        S([k_data_pilot(j)]+1)=pilindex([k_data_pilot(j)]+1)
        for kk=k_center(j)-10+1:k_center(j)+5
          
          if b(nb)==1
              S(kk+1)=1
          else
              S(kk+1)=-1  
          end 
           nb=nb+1
       end 
       disp(S)
       sub_S=sub_S+S;
       s=ifft(S);
       s_sum=s_sum+ s ;
       h_shift=h_n.*exp(2i*%pi*k_center(j)/Nfft*n)
       H_Shift=fft([h_shift,zeros(1,Nfft-taps)])
       subplot(4,1,4);title('FREQUENCY RESPONSE OF H_Shift');plot(abs((H_Shift)))
       xlabel('Frequency Bin (k)');ylabel('|H(k)|')
       x=convol(s,h_shift) 
       x_sum=x_sum+x
       X=fft([x zeros(1,Nfft-taps+1)]);
       //delay compensation
       k=0:255
       X_comp=X.*exp(2i*%pi*((k-k_center(j)+Co)/(2*Nfft))*((taps-1)/2))
       X_sum=X_sum+X_comp
       
    end
    
end
figure

subplot(2,1,1);
title('TIME DOMAIN BEFORE FILTERING');
plot(s_sum);
xlabel('Sample Index (n)');ylabel('Amplitude')

subplot(2,1,2);
title('FREQUENCY DOMAIN BEFORE FILTERING');
s_fft=fft(s_sum)
plot(s_fft);
xlabel('Frequency Bin (k)');ylabel('Magnitude')

figure
subplot(4,1,1);
title('TIME DOMAIN BEFORE DELAY COMPENSATION');
plot(x_sum((taps-1)/2:(taps-1)/2+Nfft-1))
xlabel('Sample Index (n)');ylabel('Amplitude')

subplot(4,1,2);
title('FREQUENCY DOMAIN BEFORE DELAY COMPENSATION');
x_sum_fft=fft(x_sum)
plot(x_sum_fft);
xlabel('Frequency Bin (k)');ylabel('Magnitude')

xsum=ifft([X_sum ])
xsum_fft=X_sum
xsum_fft_down=xsum_fft(1:2:end)

subplot(4,1,3);
title('TIME DOMAIN AFTER DELAY COMPENSATION');
plot(xsum(1:128))
xlabel('Sample Index (n)');ylabel('Amplitude')

subplot(4,1,4);
title('FREQUENCY DOMAIN AFTER DELAY COMPENSATION');
plot(xsum_fft(1:2:end));
xlabel('Frequency Bin (k)');ylabel('Magnitude')

//channel estimation and equivalisation 


kp = zeros(1, subband);
Hp = zeros(1, subband);
Xeq_pilot_all=[]
for p=1:subband
    p_idx = k_data_pilot(p)+1;
    kp(p) = p_idx;
    Hp(p) = xsum_fft_down(p_idx)/sub_S([k_data_pilot(p)]+1);  //pilot estimation
    Xeq_pilot=xsum_fft_down(p_idx)/Hp(p)       //equivalsing the pilot
    Xeq_pilot_all=[Xeq_pilot_all Xeq_pilot]
end

k = 1:Nfft;
Hdata = interp1(kp, Hp, k, 'linear', 'extrap'); //data estimation
Xeq_data=xsum_fft_down./Hdata   //data equivalising
figure

subplot(2,1,2)
title('equalised signal') 
plot(sign(Xeq_data))
xlabel('Frequency Bin (k)');ylabel('Magnitude')

subplot(2,1,1);
title('FREQUENCY DOMAIN BEFORE FILTERING');
s_fft=fft(s_sum)
plot(s_fft);
xlabel('Frequency Bin (k)');ylabel('Magnitude')

real([sub_S; round(s_fft*100)/100; round(xsum_fft_down*100)/100; round(abs(xsum_fft_down)*100)/100; round(unwrap(angle(xsum_fft_down))*100)/100;sign(Xeq_data)])




























