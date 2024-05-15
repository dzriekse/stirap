%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STIRAP with counterintuitive pulse sequence, Schroedinger & Bloch ansatz.
%    ---   |2>
%	 / \
%	/  --- |3>
% ---	   |1>
%
% Philippe W. Courteille, 19.1.2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
figure(1);      clf;        set(gcf,'Position',[50 100 1200 300]);

G12 = 1;                                                                    % decays scaled to the 1-2 transition
G23 = .5*G12;
G31 = 2e-3*G12;
d12 = 0*G12;                                                                % detunings
d23 = 0*G12;
T   = 2/G12;                                                                % time scale
t    = [-1:.001:1]*40;
o12 =  atan(t/T)/pi+.5;                                                     % Rabi frequencies
o23 = -atan(t/T)/pi+.5;


% Schroedinger

y = [1 0 0]';                                                               % initial condition
for jj=1:length(t)
	H = [	0         o12(jj)/2       0
		  o12(jj)/2     d12       o23(jj)/2
		  0           o23(jj)/2	   d12-d23		];

	y = expm(1i*H*(t(2)-t(1)))*y;
	R1(jj) = y(1)*y(1)';                                                    % populations
	R2(jj) = y(2)*y(2)';
	R3(jj) = y(3)*y(3)';
end

subplot(141),plot(t,o12,'r',t,o23,'b','LineWidth',1.2);
set(gca,'xLim',[-1 1]*max(t),'TickLabelInterpreter','latex','FontSize',18,'LineWidth',1)
xlabel('$\Gamma_{12}t$','interpreter','latex','Color','k','FontSize',20);
ylabel('$\Omega_{12} ~~,~~ \Omega_{23}$','interpreter','latex','Color','k','FontSize',20);
text(-67,1,'\fontsize{16}(a)')

subplot(142),plot(t,R1,'r',t,R2,'g',t,R3,'b','LineWidth',1.2);
set(gca,'xLim',[-1 1]*max(t),'TickLabelInterpreter','latex','FontSize',18,'LineWidth',1)
xlabel('$\Gamma_{12}t$','interpreter','latex','Color','k','FontSize',20);
ylabel('$\rho_{kk}$','interpreter','latex','Color','k','FontSize',20);
text(-67,1,'\fontsize{16}(b)')


% Bloch

n = [1 0 0 0 0 0 0 0 0]';                                                   % initial condition
for jj=1:length(t)
	O12 = o12(jj);
	O23 = o23(jj);

M = [0	    G12    G31	    0	-O12	  0		0	 0		0
     0 -G12-G23      0	    0	 O12	  0   -O23	 0		0
     0	    G23   -G31	    0	   0	  0    O23	 0		0

     0	      0      0 -G12/2	 d12	  0		0	 0   O23/2
 O12/2	 -O12/2      0	 -d12 -G12/2	  0		0   O23/2	 0

     0	      0      0	    0	   0 -G23/2   -d23	 0  -O12/2
     0	  O23/2 -O23/2	    0	   0	d23 -G23/2  -O12/2	 0

     0	      0      0	    0 -O23/2	  0  O12/2  -G31/2 d23-d12
     0	      0      0 -O23/2	   0  O12/2	 0      d12-d23  -G31/2];

    n = expm(M*(t(2)-t(1)))*n;
    nt1(jj) = n(1);
    nt2(jj) = n(2);
    nt3(jj) = n(3);
end

subplot(143),plot(t,nt1,'r',t,nt2,'g',t,nt3,'b','LineWidth',1.2);
set(gca,'xLim',[-1 1]*max(t),'yLim',[0 1],'TickLabelInterpreter','latex','FontSize',18,'LineWidth',1)
xlabel('$\Gamma_{12}t$','interpreter','latex','Color','k','FontSize',20);
ylabel('$\rho_{kk}$','interpreter','latex','Color','k','FontSize',20);
text(-67,1,'\fontsize{16}(c)')

% loop thru different frequencies, do schrodinger on each, graph
% first try it once

G12 = 1;                                                                    % decays scaled to the 1-2 transition
G23 = .5*G12;
G31 = 2e-3*G12;
d12 = 0*G12;                                                                % detunings
d23 = 0*G12;
T   = 2/G12;                                                                % time scale
t    = [-1:.001:1]*40;
o12a = atan(t/T)/pi-0.5;                                                     % Rabi frequencies
o23a = -atan(t/T)/pi-0.5;

% Schroedinger

y = [1 0 0]';                                                               % initial condition
for jj=1:length(t)
	H = [	0         o12a(jj)/2       0
		  o12a(jj)/2     d12       o23a(jj)/2
		  0           o23a(jj)/2	   d12-d23		];

	y = expm(1i*H*(t(2)-t(1)))*y;
	R1(jj) = y(1)*y(1)';                                                    % populations
	R2(jj) = y(2)*y(2)';
	R3(jj) = y(3)*y(3)';
end

n = [1 0 0 0 0 0 0 0 0]';                                                   % initial condition
for jj=1:length(t)
	O12a = o12a(jj);
	O23a = o23a(jj);

M = [0	    G12    G31	    0	-O12a	  0		0	 0		0
     0 -G12-G23      0	    0	 O12a	  0   -O23a	 0		0
     0	    G23   -G31	    0	   0	  0    O23a	 0		0

     0	      0      0 -G12/2	 d12	  0		0	 0   O23a/2
 O12a/2	 -O12a/2      0	 -d12 -G12/2	  0		0   O23a/2	 0

     0	      0      0	    0	   0 -G23/2   -d23	 0  -O12a/2
     0	  O23a/2 -O23a/2	    0	   0	d23 -G23/2  -O12a/2	 0

     0	      0      0	    0 -O23a/2	  0  O12a/2  -G31/2 d23-d12
     0	      0      0 -O23a/2	   0  O12a/2	 0      d12-d23  -G31/2];

    n = expm(M*(t(2)-t(1)))*n;
    nt1(jj) = n(1);
    nt2(jj) = n(2);
    nt3(jj) = n(3);
end

subplot(144),plot(t,nt1,'r',t,nt2,'g',t,nt3,'b','LineWidth',1.2);
set(gca,'xLim',[-1 1]*max(t),'yLim',[0 1],'TickLabelInterpreter','latex','FontSize',18,'LineWidth',1)
xlabel('$\Gamma_{12}t$','interpreter','latex','Color','k','FontSize',20);
ylabel('$\rho_{kk}$','interpreter','latex','Color','k','FontSize',20);
text(-67,1,'\fontsize{16}(d)')