clc;
clear all;
close all;
data=[0 0 1 0 1 1 1 0 0 1 0 1 0 0 0 1 1];%given message signal
a=[];
b=[];
L=length(data);
for i=1:L
    if data(i)==1
        a=ones(1,100);
    else
        a=-ones(1,100);
    end
    b=[b a];
end

odd=data(1:2:L);
x=[];
o=[];
m1=length(odd);
for i=1:m1
    if odd(i)==1;
        x=ones(1,2*100);
    else 
        x=-ones(1,2*100);
    end
       o=[o x]
end
even=data(2:2:L);
z=[];
e=[];
m2=length(even);
for i=1:m2
    if even(i)==1;
        z=ones(1,2*100);
    else 
        z=-ones(1,2*100);
    end
       e=[e z]
end

t1=linspace(0,m1,2*m1*100);
t2=linspace(0,m2,2*m2*100)
f=2;
c1=cos(2*pi*f.*t1);
c2=sin(2*pi*f.*t2);

s=o.*c1;
p=e.*c2;

n1=1:1800;
y1=s;
n2=1:1600;
y2=p;
u=min(min(n1),min(n2));
t=max(max(n1),max(n2));
r=u:1:t;
z1=[];
temp=1;
for i=1:length(r)
    if(r(i)<min(n1) || r(i)>max(n1))
      z1=[z1 0];
    else
        z1=[z1 y1(temp)];
        temp=temp+1;
    end
end
z2=[];
temp=1;
for i=1:length(r)
    if(r(i)<min(n2) || r(i)>max(n2))
      z2=[z2 0];
    else
        z2=[z2 y2(temp)];
        temp=temp+1;
    end
end
z=z1+z2;

figure(1)
subplot(311)
plot(b)
xlabel('n')
title('Message signal')
axis([1 100*L -1.5 1.5])
subplot(312)
plot(o);
xlabel('n');
title('odd signal');
axis([0 100*L -1.5 1.5]);

subplot(313)
plot(e);
xlabel('n');
title('even signal');
axis([0 100*L -1.5 1.5]);
figure(2)
subplot(311)
plot(s);
xlabel('n');
title('odd*cos');

subplot(312)
plot(p);
xlabel('n');
title('even*cos');

subplot(313)
plot(z);
xlabel('n');
title('Quadrature Phase shift keying signal');

y = fft(z);
figure(3)
plot(fftshift(abs(y)));
title('spectrum of modulated signal');