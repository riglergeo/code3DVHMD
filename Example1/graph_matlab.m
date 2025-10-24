# simple Matlab code for simple visualization of
# normalized components by H0

clear all;
hy=load('hy.dat');
hz=load('hz.dat');

s=3.66 #spacing between coil
H0 = (-4*pi*(s^3));
hy(:,2)=hy(:,2).*H0;
hy(:,3)=hy(:,3).*H0;
hz(:,2)=hz(:,2).*H0;
hz(:,3)=hz(:,3).*H0;

figure(1),plot(hy(:,1),hy(:,2),'b','linewidth', 2);
figure(2),plot(hz(:,1),hz(:,2),'k','linewidth', 2);
figure(3),plot(hy(:,1),hy(:,3),'b','linewidth', 2);
figure(4),plot(hz(:,1),hz(:,3),'k','linewidth', 2);


