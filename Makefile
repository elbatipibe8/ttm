# The makefile compiles (type make) and makes an executable called mse-bs

EXEPATH = .

EXE = $(EXEPATH)/Al.x
FC = mpif90


#FLAGS=-fmax-stack-var-size=10000000000


#F90FLAGS = -tp k8-64 -O -mcmodel=medium -Mlarge_arrays
#F90FLAGS =  -bmaxdata:0x80000000 -bmaxstack:0x40000000 -g -qfullpath
#FCFLAGS =  -I./ -O3 -march=corei7-avx
# -g -C -qsigtrap

FCFLAGS =  -I./ -O3 -mtune=corei7-avx

#LDFLAGS = -tp k8-64 -mcmodel=medium
#LDFLAGS =  -bmaxdata:0x80000000 -bmaxstack:0x40000000
#LDFLAGS =
# -g -C -qsigtrap


OBJS = Main.o Eftab-EAM.o Eftab-pair.o Forces.o F-EAM.o EAM-Density.o F-TTM.o Flushout.o Leng.o Nord5.o OpenFiles.o ReadFiles.o SetInit.o Temper.o WriteEnd.o WriteInit.o Swrite.o Swrite-Res.o SwriteTTM.o SwriteTTM-Res.o Solid.o TTMInit.o TTMDiffuse.o TTMTemper.o Pressure.o SetTTM.o Liquid.o FirstList.o LinkList.o Share.o Share-EAM.o Share-TTM.o Share-Te.o NbList-MPI.o Migrate.o Bcast-Pos.o Collect.o SetBulk.o Share-Bulk.o Share-Bulk2.o F-sub.o F-EAM-br.o F-EAM-br-test.o F-pair.o Share-Q1V.o Smooth-Q1V.o Share-Tl.o Smooth-Tl.o Share-PR.o Smooth-PR.o Symmetry.o Bere-PT.o CheckWrite.o Coeff.o Share-TlD.o Smooth-TlD.o Share-DSL.o Smooth-PRW.o

# compile and load
.SUFFIXES: .f90 .o
.f90.o:
	$(FC) $(FCFLAGS) -c $*.f90

default:
	@echo " "
	@echo "Compiling Code Au structures"
	@echo "Version 1.0"
	@echo "FORTRAN 90"
	$(MAKE) $(EXE)

$(OBJS): common.h commonTTM.h parameters.h
$(EXE):	$(OBJS)
	$(FC) $(LDFLAGS) -o $(EXE) $(OBJS)

clean:	
	rm -f *.o
	rm -f Al.x
