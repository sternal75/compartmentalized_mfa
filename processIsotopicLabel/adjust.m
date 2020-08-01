function [R,QunlityControl]=adjust(data,N,prob_natural_C13, prob_labeled_carbon_is_C12)  

% data= measuring data with 1st column: indicator of labeling number ,2nd
% column: signal
% N: number of carbon in the compound

%q13=0.008;   %C12 inpurity in labeled compound
% could use Ac-coA labeling as inpurity



M=data;

O=zeros(N+1,N+1);

for i=0:N   %i- number of labeled C

    X=0:N-i;
    C12=binopdf(X,N-i,prob_natural_C13)';
    Y=0:i;
    C13=binopdf(Y,i,1-prob_labeled_carbon_is_C12)'; 
    
    V=condense(C12,C13);
    O(:,i+1)=V;

end 

%O(N+2,:)=ones(1,N+1);
%R=inv(O)*M;
%R=(M\O)';


%R=mldivide([O;ones(1,N+1)],[M;1]);
%QunlityControl=sum(abs(O*R-M));

%[R,QunlityControl]=lsqcurvefit(@join,ones(N+1,1)/N+1+1,[O;10*ones(1,N+1)],[M;10],zeros(N+1,1),inf*ones(N+1,1));

xdata = [O;10*ones(1,N+1)]; 
ydata = [M;10];
[R,QunlityControl] = fmincon(@(x) my_func(x, xdata, ydata),  ones(N+1,1)/N+1+1,[], [], [], [], zeros(N+1,1),inf*ones(N+1,1));

R=R/sum(R);


%  mldivide \, is better

function err = my_func(x, xdata, ydata)
    y = xdata*x;
    err = sum( (y-ydata).^2 );



function MDV=condense(part1,part2)

l1=length(part1);
l2=length(part2);

MDV=zeros(l1+l2-1,1);
v=zeros(l1+l2-1,1);
v(1:l2)=part2;
O=zeros(l1+l2-1,l1+l2-1);
for i=1:l1   % i-1 is the shift
   for j=i:l1+l2-1
       O(j,j-i+1)=part1(i);
   end
end

MDV=O*v;



%test = [         0    0.4891
 %   1.0000    0.0188
  %  2.0000    0.4921]