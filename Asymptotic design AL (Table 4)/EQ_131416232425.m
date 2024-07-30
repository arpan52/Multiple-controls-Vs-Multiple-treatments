function L = EQ_131416232425(y)

x=y';

n= sqrt(57);
del= 6.6;
sigma = 5.9;

 e1=[1,0,0,0,0,0]';
 e2=[0,1,0,0,0,0]';
 e3=[0,0,1,0,0,0]';
 e4=[0,0,0,1,0,0]';
 e5=[0,0,0,0,1,0]';
 e6=[0,0,0,0,0,1]';

 e13=e1+e3;
 e14=e1+e4;
 e16=e1+e6;
 e23=e2+e3;
 e24=e2+e4;
 e25=e2+e5;

 A13=e1*e3';
 A14=e1*e4';
 A16=e1*e6';
 A23=e2*e3';
 A24=e2*e4';
 A25=e2*e5';



 mu13 = sqrt((x'*A13*x)/(x'*e13));
 mu14 = sqrt((x'*A14*x)/(x'*e14));
 mu16 = sqrt((x'*A16*x)/(x'*e16));
 mu23 = sqrt((x'*A23*x)/(x'*e23));
 mu24 = sqrt((x'*A24*x)/(x'*e24));
 mu25 = sqrt((x'*A25*x)/(x'*e25));


p13=normpdf(norminv(0.95)- (n*del/sigma)*mu13)-normpdf(-norminv(0.95)- (n*del/sigma)*mu13);
p14=normpdf(norminv(0.95)- (n*del/sigma)*mu14)-normpdf(-norminv(0.95)- (n*del/sigma)*mu14);
p16=normpdf(norminv(0.95)- (n*del/sigma)*mu16)-normpdf(-norminv(0.95)- (n*del/sigma)*mu16);
p23=normpdf(norminv(0.95)- (n*del/sigma)*mu23)-normpdf(-norminv(0.95)- (n*del/sigma)*mu23);
p24=normpdf(norminv(0.95)- (n*del/sigma)*mu24)-normpdf(-norminv(0.95)- (n*del/sigma)*mu24);
p25=normpdf(norminv(0.95)- (n*del/sigma)*mu25)-normpdf(-norminv(0.95)- (n*del/sigma)*mu25);

L13_1 = e13'*x*x'*A13';
L13_2 = x'*A13*(e13'*x*eye(6)-x*e13');
L13_3 = 2*sqrt(x'*A13*x)*(x'*e13)^(3/2);
L13 = (L13_1+L13_2)/L13_3;

L14_1 = e14'*x*x'*A14';
L14_2 = x'*A14*(e14'*x*eye(6)-x*e14');
L14_3 = 2*sqrt(x'*A14*x)*(x'*e14)^(3/2);
L14 = (L14_1+L14_2)/L14_3;

L16_1 = e16'*x*x'*A16';
L16_2 = x'*A16*(e16'*x*eye(6)-x*e16');
L16_3 = 2*sqrt(x'*A16*x)*(x'*e16)^(3/2);
L16 = (L16_1+L16_2)/L16_3;

L23_1 = e23'*x*x'*A23';
L23_2 = x'*A23*(e23'*x*eye(6)-x*e23');
L23_3 = 2*sqrt(x'*A23*x)*(x'*e23)^(3/2);
L23 = (L23_1+L23_2)/L23_3;

L24_1 = e24'*x*x'*A24';
L24_2 = x'*A24*(e24'*x*eye(6)-x*e24');
L24_3 = 2*sqrt(x'*A24*x)*(x'*e24)^(3/2);
L24 = (L24_1+L24_2)/L24_3;

L25_1 = e25'*x*x'*A25';
L25_2 = x'*A25*(e25'*x*eye(6)-x*e25');
L25_3 = 2*sqrt(x'*A25*x)*(x'*e25)^(3/2);
L25 = (L25_1+L25_2)/L25_3;

L1 = (p13*L13)+(p14*L14)+(p16*L16)+(p23*L23)+(p24*L24)+(p25*L25);
L = abs(L1(1))+abs(L1(2))+abs(L1(3))+abs(L1(4))+abs(L1(5))+abs(L1(6));

%L = abs([L1(1),L1(2),L1(3),L1(4),x'*ones(4,1)-1]);
%L =abs(L1(1))+abs(L1(2))+abs(L1(3))+abs(L1(4))+abs(x'*[1,0,-1,0]')+abs(x'*[0,1,0,-1]')
end
