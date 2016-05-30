function h = ellips(x0,y0,P,col)
D=det(P);
teta=0.5*atan2(2*P(1,2),P(1,1)-P(2,2));
%rho = P(1,2)/sqrt(P(1,1)*P(2,2));
A = sqrt(D/(P(2,2)*(cos(teta))^2 + P(1,1)*(sin(teta))^2 - P(1,2)*sin(2*teta)));
B = sqrt(D/(P(2,2)*(sin(teta))^2 + P(1,1)*(cos(teta))^2 + P(1,2)*sin(2*teta)));
%sqrt(9.21) = 99%  sqrt(5.99) = 95%, sqrt(4.61) = 90%

a = sqrt(4.61)*A;
b = sqrt(4.61)*B;

% a = sqrt(9.21)*A;
% b = sqrt(9.21)*B;

i=1;
for t=0:0.3:2*pi+0.2;
	x(i)=a*cos(t)*cos(teta) - b*sin(t)*sin(teta);
	y(i)=a*cos(t)*sin(teta) + b*sin(t)*cos(teta);
	i = i+1;
end
h = plot((x+x0),(y+y0),'color',col);

end