
OBJ_MG=FiniteDifferenceStencil.o Prolongation.o MultiGridLevel.o MultiGridSetting.o
OBJ_OCT=mgsolve.o $(OBJ_MG)

OBJ=$(OBJ_OCT)

default: mgsolve.oct
	
clean:
	rm -f mgsolve.oct $(OBJ) 
clean-obj:
	rm -f $(OBJ) 

mgsolve.oct: $(OBJ_OCT)
	$(MKOCTFILE) -o mgsolve.oct $(OBJ_OCT)
FiniteDifferenceStencil.o:FiniteDifferenceStencil.h
Prolongation.o:Prolongation.h
MultiGridLevel.o:MultiGridLevel.h
MultiGridSetting.o:MultiGridSetting.h
%.o: %.cpp
	$(MKOCTFILE) -c $< 
