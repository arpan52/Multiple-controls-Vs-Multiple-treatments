function L = EQ_132324252627(y)

x=y';

n= sqrt(57);
del= 6.6;
sigma = 5.9;

 e1=[1,0,0,0,0,0,0]';
 e2=[0,1,0,0,0,0,0]';
 e3=[0,0,1,0,0,0,0]';
 e4=[0,0,0,1,0,0,0]';
 e5=[0,0,0,0,1,0,0]';
 e6=[0,0,0,0,0,1,0]';
 e7=[0,0,0,0,0,0,1]';

 e13=e1+e3;
 e23=e2+e3;
 e24=e2+e4;
 e25=e2+e5;
 e26=e2+e6;
 e27=e2+e7;


 A13=e1*e3';
 A23=e2*e3';
 A24=e2*e4';
 A25=e2*e5';
 A26=e2*e6';
 A27=e2*e7';


 mu13 = sqrt((x'*A13*x)/(x'*e13));
 mu23 = sqrt((x'*A23*x)/(x'*e23));
 mu24 = sqrt((x'*A24*x)/(x'*e24));
 mu25 = sqrt((x'*A25*x)/(x'*e25));
 mu26 = sqrt((x'*A26*x)/(x'*e26));
 mu27 = sqrt((x'*A27*x)/(x'*e27));

p13=normpdf(norminv(0.95)- (n*del/sigma)*mu13)-normpdf(-norminv(0.95)- (n*del/sigma)*mu13);
p23=normpdf(norminv(0.95)- (n*del/sigma)*mu23)-normpdf(-norminv(0.95)- (n*del/sigma)*mu23);
p24=normpdf(norminv(0.95)- (n*del/sigma)*mu24)-normpdf(-norminv(0.95)- (n*del/sigma)*mu24);
p25=normpdf(norminv(0.95)- (n*del/sigma)*mu25)-normpdf(-norminv(0.95)- (n*del/sigma)*mu25);
p26=normpdf(norminv(0.95)- (n*del/sigma)*mu26)-normpdf(-norminv(0.95)- (n*del/sigma)*mu26);
p27=normpdf(norminv(0.95)- (n*del/sigma)*mu27)-normpdf(-norminv(0.95)- (n*del/sigma)*mu27);

L13_1 = e13'*x*x'*A13';
L13_2 = x'*A13*(e13'*x*eye(7)-x*e13');
L13_3 = 2*sqrt(x'*A13*x)*(x'*e13)^(3/2);
L13 = (L13_1+L13_2)/L13_3;

L23_1 = e23'*x*x'*A23';
L23_2 = x'*A23*(e23'*x*eye(7)-x*e23');
L23_3 = 2*sqrt(x'*A23*x)*(x'*e23)^(3/2);
L23 = (L23_1+L23_2)/L23_3;

L24_1 = e24'*x*x'*A24';
L24_2 = x'*A24*(e24'*x*eye(7)-x*e24');
L24_3 = 2*sqrt(x'*A24*x)*(x'*e24)^(3/2);
L24 = (L24_1+L24_2)/L24_3;

L25_1 = e25'*x*x'*A25';
L25_2 = x'*A25*(e25'*x*eye(7)-x*e25');
L25_3 = 2*sqrt(x'*A25*x)*(x'*e25)^(3/2);
L25 = (L25_1+L25_2)/L25_3;

L26_1 = e26'*x*x'*A26';
L26_2 = x'*A26*(e26'*x*eye(7)-x*e26');
L26_3 = 2*sqrt(x'*A26*x)*(x'*e26)^(3/2);
L26 = (L26_1+L26_2)/L26_3;

L27_1 = e27'*x*x'*A27';
L27_2 = x'*A27*(e27'*x*eye(7)-x*e27');
L27_3 = 2*sqrt(x'*A27*x)*(x'*e27)^(3/2);
L27 = (L27_1+L27_2)/L27_3;

L1 = (p13*L13)+(p23*L23)+(p24*L24)+(p25*L25)+(p26*L26)+(p27*L27);
L = abs(L1(1))+abs(L1(2))+abs(L1(3))+abs(L1(4))+abs(L1(5))+abs(L1(6))+abs(L1(7));

%L = abs([L1(1),L1(2),L1(3),L1(4),x'*ones(4,1)-1]);
%L =abs(L1(1))+abs(L1(2))+abs(L1(3))+abs(L1(4))+abs(x'*[1,0,-1,0]')+abs(x'*[0,1,0,-1]')
end
