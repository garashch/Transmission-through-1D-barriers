
# Simple makefile

  I77=gfortran

##LIB = -llapack
.SUFFIXES:	.f .x .o .gx

.f.o:
	$(I77)  $(Flags) -c  $<
.f.x:
	$(I77)  -O4 $(LIB) $(Flags)  $< -o $*.x
.f.gx:
	$(I77)  $(LIB) $(Flags) -g  $< -o $*.gx
