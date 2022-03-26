OBJ = nuclear.o eos.o opacity.o ludcmp.o zbrent.o ion.o mixlen.o fdjac.o
SHOOT = zams_shoot.o newt.o odeint.o $(OBJ)
RELAX = zams.o solvde.o $(OBJ)

fig4: fig4.o $(SHOOT)
	g++ fig4.o $(SHOOT)
fig5: fig5.o $(SHOOT)
	g++ fig5.o $(SHOOT)
fig6: fig6.o $(SHOOT)
	g++ fig6.o $(SHOOT)
init: zams_init.o zams.o solvde.o $(SHOOT)
	g++ zams_init.o zams.o solvde.o $(SHOOT)
fig7-9: fig7-9.o $(RELAX)
	g++ fig7-9.o $(RELAX)
