OBJ = zams.o nuclear.o eos.o opacity.o newt.o ludcmp.o odeint.o zbrent.o

fig1: fig1.o nuclear.o
	g++ fig1.o nuclear.o
fig2: fig2.o eos.o opacity.o
	g++ fig2.o eos.o opacity.o
fig3: fig3.o $(OBJ)
	g++ fig3.o $(OBJ)
