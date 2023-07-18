function [hh,m,sem]=plot_sem(arg1,arg2,arg3,varargin)
%plot patch  mean plus minus sed, and array (Ntr x time) has to be
%provided, average is calcuated over first dimension (trials)
%can be called
%plot(y), just an array 
%plot(time, y), x axis is given
%plot(time, y, 'r'), colour is given as well
%plot(y,varargin), varargin has arguments for matlab function patch
%plot(time,y,varargin), varargin has arguments for matlab function patch

if nargin<1
    fprintf('Array input has to be provided, (Ntr x time)\n');
    return;
end
%if only one argument is given assume it is y
if nargin<2
    y=arg1; x=1:size(y,2);c='r'; 
elseif nargin<3 
    %if two arguments are given check they type of second argument
    x=arg1; y=arg2; c='r';
else
    x=arg1; y=arg2; c=arg3;
%     t=varargin;
%     varargin{1}=arg2; varargin{2}=arg3;
%     for i=1:length(t), varargin{i+2}=t{i}; end
end

%check varargin
n=length(varargin);
if n/2~=round(n/2) 
    fprintf('Check arguments\n');
    return;
end
b=0;
for i=1:2:n, if ~ischar(varargin{i}), b=1; break; end; end
if b
    fprintf('Check arguments\n');
    return;
end
args='';
smo=0;
i=1;
while i<=length(varargin)
    t=varargin{i};
    if strcmp('smooth',t)
        smo=varargin{i+1};
        i=i+2;
    else        
        if ischar(t), args=sprintf('%s,''%s''',args,t);
        else
            args=sprintf('%s,%d',args,t);
        end
        i=i+1;
    end
end

nn=size(y,2);
[NN]=size(y,1);
m=nanmean(y);
sem=NaN*ones(1,nn);
for ii=1:nn
    sem(ii)=nanstd(y(:,ii))/sqrt(sum(~isnan(y(:,ii))));
end
sem=nanstd(y)/sqrt(NN);
%just add NaN in the beginning and in the end
x=[NaN x NaN];
m=[NaN m NaN];
sem=[NaN sem NaN];

t=find(isnan(m)); tt1=t(diff(t)~=1)+1;
t=find(~isnan([m 0])); tt2=t(diff(t)~=1);
h0=NaN*ones(1,length(tt1));
for i=1:length(tt1)
    ind=tt1(i):tt2(i);
    xx=x(ind); mm=m(ind); semm=sem(ind);
    if smo, mm=smooth(mm,smo)'; semm=smooth(semm,smo)'; end
    xxx=[xx(1:end) xx(end:-1:1)];
    yyy=[mm+semm mm(end:-1:1)-semm(end:-1:1)];
    h0(i)=patch(xxx,yyy,c,'edgecolor','none','facealpha',0.5);
end

for i=1:length(tt1)
    ind=tt1(i):tt2(i);
    xx=x(ind); mm=m(ind); semm=sem(ind);
    if smo, mm=smooth(mm,smo)'; semm=smooth(semm,smo)'; end
    line(xx,mm,'color',c);
end

if ~isempty(args), eval(sprintf('set(h0%s);',args)); end
if nargout>=1, hh=h0; end

end