% carga las concentraciones iniciales de cada componente
%Conc factor eres 1:n-1 ejemplo 2 factores, 3 vectores blanco (se dan n-1
%vectores)

nes=input('para cuales espectros tiene blancos? ');

for ii=1:n-1
    esp=[];
for jj=1:length(nes);

   eval(['esp= [ esp ' Mdatos '(:,nes(jj))];']);
    aa=num2str(ii);
    a=num2str(nes(jj));
b(ii,jj)=input(['ingrese los valores blancos del factor ' aa ' para el espectro ' a '  ']) ;
end
b=b';
bb=1-sum(b,2);

b=[b bb];
b=b';
end
b=pinv(b);
xx=esp*b;
T=Mdiagonal*R'*xx;

    
    
% calcula la matriz de factores reales
X = R*T;

% calcula la matriz de concentraciones
Y = T^(-1)*v';


f = figure('Units','characters','Position',[  10   38.0  112.4000   25]);
plot(nu,X,'-')
ylabel('Componentes principales reales')
xlabel('Numero de onda [cm^{-1}]')


YY=Y';

f = figure('Units','characters','Position',[  10   4.0  112.4000   25]);
eval(['plot(' tiempo ' ,YY, ''.'')'])
ylabel('Concentraciones principales reales')
xlabel('Tiempo')
RRR=X*Y;
f = figure('Units','characters','Position',[  162   38.0  112.4000   25]);
plot(nu,RRR,'-')
ylabel('Componentes principales reales')
xlabel('Tiempo')
RMS=norm(dr-RRR)/norm(dr)
figure(1)
