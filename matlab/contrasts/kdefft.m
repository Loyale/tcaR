function [xa,y]=kdfft1(x,inc,h);

% KDFFT1      KDE plot for 1D data using FFTs.
%
%             Example: load kddata
%                      [x,y]=kdfft1(suicide,'knorm',128,40);
%
%             KDFFT1(A,K,INC,H,G,XR) plots the Kernel Density Estimate with 
%             Normal kernel.  
%             A: The data is stored in the rows of A.  For 1D data A 
%                will have one column.
%             K: The kernel function, one of 'KNORM', 'KEPAN', etc..
%             INC: The number of increments to use in plotting the 
%                  density function.  THIS MUST BE A POWER OF 2, for example 
%                  32, 64 and 128 give good results.  
%             H: Window width.  If omitted, an 'optimal' value (calculated using
%                the STE rule) is used.
%             G: Optional, G=0 turns graphics off.
%             XR: Optional, XR=[XMIN, XMAX] specifies the x range - XR must 
%                 include the range of the data.
%
%             Christian C. Beardah 1994

P=2*inc;

r=4;

n=length(x);      
xmin=min(x);      % Find the minimum value of x.
xmax=max(x);      % Find the maximum value of x.
xrange=xmax-xmin; % Find the range of x.

% xa holds the x 'axis' vector, defining a grid of x values where 
% the k.d. function will be evaluated and plotted.

    ax=xmin-xrange/2;   
    bx=xmax+xrange/2;

xa=linspace(ax,bx,inc); 

c=zeros(inc,1);

deltax=(bx-ax)/(inc-1);
binx=floor((x-ax)/deltax)+1;

% Obtain the grid counts.

for i=1:n, % Loop over data points in x direction.
  c(binx(i))=c(binx(i))+(xa(binx(i)+1)-x(i))/deltax;
  c(binx(i)+1)=c(binx(i)+1)+(x(i)-xa(binx(i)))/deltax;  
end;

% Obtain the kernel weights.

L=min(floor(r*h*(inc-1)/(bx-ax)),inc-1);

kw=zeros(1,2*inc);

s= (bx-ax)*[0:L]/((inc-1)*h) ;

temp=(2*pi)^(-.5)*exp(-0.5*s.*s)/(n*h);           

kw(inc+1:inc+L+1)=temp;
kw(inc-L+1:inc)=temp([L+1:-1:2]);

% Apply 'fftshift' to kw.

kw=fftshift(kw)';

% Perform the convolution.

y=real(ifft(fft(c,P).*fft(kw)));

y=y(1:inc).*(y(1:inc)>0);
