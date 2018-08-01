close all;clear all;
load('MO1_117_03.mat');

[svody,delka] = size(x); %pocet svodu EKG
if svody >= 6 %pocet sloupcu pro subplot
  sloupce = 2;
  radky = svody/2;
else
  sloupce = 1;
  radky = svody;
end

for i = 1:svody %vykresleni svodu puvodniho signalu EKG
    subplot(radky,sloupce,i)
    plot(x(i,:));
    title(['pùvodní signál (svod ' num2str(i) ')'])
    xlim([1 delka])
end

%% ______1. exrakce znakového signálu____
%______filtrace pasmovou propusti______
h = fir1(27,[0.072 0.14],'bandpass'); %impulzni char. PP
a = 1;
xf = zeros(sloupce,delka);
figure(2)
for i = 1:svody
    xf(i,:) = filtfilt(h,a,x(i,:)); %filtrovany signal PP
    subplot(radky,sloupce,i)
    plot(xf(i,:));
    title(['signál po filtraci PP(svod ' num2str(i) ')'])
    xlim([1 delka])
end

%______nelineární transformace_____
y = zeros(sloupce,delka);
figure(3)
for i = 1:svody
    for m = 1:length(xf(i,:))
        y(i,m) = sign(xf(i,m))*(xf(i,m))^2; %upraveny signal
    end
  subplot(radky,sloupce,i)
  plot(y(i,:));
  title(['signál po nelineární filtraci (svod ' num2str(i) ')'])
  xlim([1 delka])
end

%______amplituda sekvence_____
lamK = 0.99; %nejlip
c = 2.6;     %nejlip
K = zeros(sloupce,delka);
for i = 1:svody
    K(i,1) = abs(y(i,1))*c; 
    for n = 2:length(y(i,:))
        K(i,n) = lamK*K(i,n-1)+(1-lamK)*abs(y(i,n))*c; %vektor zhodnoceni amplitud
    end
end

%______pøidání sekvence_____
b = zeros(sloupce,delka);
z = zeros(sloupce,delka);
figure(4)
for i = 1:svody
    for n = 1:length(K(i,:))
    b(i,n)= ((-1)^n)*K(i,n); %pridavana sekvence
    end
    z(i,:)= y(i,:)+b(i,:); %signal s pridanou sekvenci
    subplot(radky,sloupce,i)
    plot(z(i,:))
    title(['signál s pøidanou sekvencí (svod ' num2str(i) ')'])
    xlim([1 delka])
end

%______detekce a poèítání prùchodù nulou____
lamD = 0.6;
d = zeros(sloupce,delka-1);
D = zeros(sloupce,delka-1); 
for i = 1:svody
    for r = 2:length(z(i,1:end))
        d(i,r-1)= abs((sign(z(i,r))-sign(z(i,r-1)))/2);
    end
    D(i,1) = d(i,1);
    for s = 2:length(d(i,:))
        D(i,s) = lamD*D(i,s-1)+(1-lamD)*d(i,s); %znakovy signal - pocet pruchodu nulou na segment
    end
end

%% ______2. detekce událostí__________
%______výpoèet adaptivního prahu____
lamFi = 0.99; %cim mensi, tim vic kopiruje prubeh znakoveho signalu D
Fi = zeros(sloupce,delka-1);
figure(5)
for i = 1:svody
    Fi(i,1) = (lamFi)*D(i,1); 
    for j = 2:length(D(i,:))
        Fi(i,j) = lamFi*Fi(i,j-1)+(1-lamFi)*D(i,j); %adaptivni prah
    end
   subplot(radky,sloupce,i)
   plot(D(i,:)) 
   hold on;
   plot(Fi(i,:),'Color','black')
   title(['Poèty prùchodù nulou a adaptivní práh (svod ' num2str(i) ')'])
   xlim([1 delka-1])
end

%______detekce událostí_____
zacatek = cell(1,svody);
konec = cell(1,svody);
udalost = cell(1,svody);
for i = 1:svody
    zacatek{i}=0;
    konec{i}=0;
    for k = 2:length(D(i,:))
    if D(i,k) < Fi(i,k)
        if D(i,k-1) >= Fi(i,k-1)
            zacatek{i}(end+1)=k;
        end
    else
        if D(i,k-1) < Fi(i,k-1)
            konec{i}(end+1)=k;
        end
    end
    end
    zacatek{i}=zacatek{i}(2:end);
    konec{i}=konec{i}(2:end);
    if length(konec{i})<length(zacatek{i})
        zacatek{i} = zacatek{i}(1:end-1);
    end
    udalost{i}(:,1)=zacatek{i};
    udalost{i}(:,2)=konec{i};
    
    %____casova prodleva mezi udalostmi___
    j = 2;
    while j <= length(udalost{i})
       if (udalost{i}(j,1) - udalost{i}(j-1,2)) < 75 %kratsi prodleva nez 150 ms odp. 75vz - rezerva
           udalost{i}(j-1,2) = udalost{i}(j,2); %slouèení událostí
           udalost{i}(j:end-1,:)=udalost{i}(j+1:end,:); %zbavení se øádku po slouèení
           udalost{i}=udalost{i}(1:end-1,:); %posunuti
           j=j-1;   %vraceni se o radek kvuli posunuti
       end
       j = j+1;
    end
end

%% ________3. lokalizace R vln________

for i = 1:svody
    for k = 1:length(udalost{i})
        zac = udalost{i}(k,1);
        kon = udalost{i}(k,2);
        if  max(abs(y(i,zac:kon))) > 1.2*max(y(i,zac:kon))%pokud je minimum vic nez 1,2 krat maximum
        [maximum(i,k), pozice(i,k)] = max(abs(y(i,zac:kon)));%R vlna pak bude zakmit dolu 
        else
        [maximum(i,k), pozice(i,k)] = max(y(i,zac:kon));
        end
        pozice(i,k)= pozice(i,k) + zac;
    end
end

%________zisk 1 vektoru pozic QRS ze vsech svodu_______
del=length(pozice);
vektor = zeros(1,svody*del); %vytvoreni vektoru z matice pozic R vln
for i = 1:svody
    vektor((i*del)-del+1:(i*del)) = pozice(i,:);
end

mez = 75;
serazeny = sort(vektor);
ind = find(serazeny);       %indexy nenulovych pozic v serazeny
serazeny = serazeny(ind);   %pro odebrani nul, ktere "zarovnavaji" matici pozic
rozdil = find((diff(serazeny) > mez)); %pozice, kde je vetsi rozdil nez mez
rozdil2 = [1, rozdil+1, length(serazeny)+1]; %pozice od ktere je dalsi skupina
n = diff(rozdil2);  %pocet cisel ve skupinach
skupiny = mat2cell(serazeny,1,n);
poziceR=0;

for i = 1:length(skupiny)
    if length(skupiny{i})>=svody
        poziceR(end+1) = round(median(skupiny{i}));
    end
end
poziceR = poziceR(2:end);

%________zobrazení pozic QRS_______
odkud = ones(length(poziceR),1).*(min(abs((x(1,:)))));
kam = ones(length(poziceR),1).*(max(x(1,:)));
poziceR = poziceR';
figure(6)
for i = 1:svody
    subplot(radky,sloupce,i)
    plot(x(i,:))
    hold on;
    line([poziceR(:,1),poziceR(:,1)]', [odkud(:,1),kam(:,1)]','Color','red')
    title(['pozice R vln v pùvodním signálu x(n) (svod ' num2str(i) ')'])
    xlim([1, delka])
end
 disp('Detekované pozice R vln v jednotlivých svodech najdete v promìnné pozice.')
 disp('Jednotný vektor pozic R vln je v promìnné poziceR.')