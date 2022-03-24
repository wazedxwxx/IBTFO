Scheme = MUSCL #Lax or MUSCL or MHM （MSUCL Hancock）
Riemann = HLL #HLL or PVRS
Limiter = minmod #minmod or VanLeer or VanAlbada
Boundary = Regular #Regular or BaseGFM or SimpleBound

PROM = REGULAR_SOD.EX
CC = nvc++
GeometrySRC=../../Src/Geometry/
UtilsSRC=../../Src/Utils/
BoundSRC=../../Src/BoundaryCondition/
AdvanceSRC=../../Src/Base/Advance/
CommonSRC=../../Src/Base/Common/
RiemannSRC=../../Src/Base/Riemann/
SlopelimiterSRC=../../Src/Base/Slope_limiter/
TimeSRC=../../Src/Base/TimeAdvance/


GeometryFILES = $(foreach dir,$(GeometrySRC),$(wildcard $(dir)/*.cpp))
UtilsFILES = $(foreach dir,$(UtilsSRC),$(wildcard $(dir)/*.cpp))
BoundFILES = $(foreach dir,$(BoundSRC)$(Boundary),$(wildcard $(dir)/*.cpp))
AdvanceFILES = $(foreach dir,$(AdvanceSRC)$(Scheme),$(wildcard $(dir)/*.cpp))
CommonFILES = $(foreach dir,$(CommonSRC),$(wildcard $(dir)/*.cpp))
RiemannFILES = $(foreach dir,$(RiemannSRC)$(Riemann),$(wildcard $(dir)/*.cpp))
SlopelimiterFILES = $(foreach dir,$(SlopelimiterSRC)$(Limiter),$(wildcard $(dir)/*.cpp))
TimeSFILES = $(foreach dir,$(TimeSRC)$(Time),$(wildcard $(dir)/*.cpp))
$(shell rm -f $(PROM))


CFLAGS = -fast -ta=multicore -Minfo=accel -I$(GeometrySRC) -I$(UtilsSRC) \
                -I$(BoundSRC) -I$(AdvanceSRC) -I$(CommonSRC) -I$(RiemannSRC)\
                 -I$(SlopelimiterSRC) -I$(TimeSRC)


$(PROM) : $(OBJ)
	$(CC) $(CFLAGS) -o $(PROM) $(GeometryFILES) $(UtilsFILES) \
	$(BoundFILES) $(AdvanceFILES) $(CommonFILES) $(RiemannFILES) $(SlopelimiterFILES) $(TimeFILES) \
	$(wildcard ./*.cpp)
	$(RM) *.o *.vtk


.PHONY : clean
clean:
		$(RM) $(OBJ) $(PROM) *.o
