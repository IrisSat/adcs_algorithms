CC=C:/msys64/mingw64/bin/gcc
FLAGS = -g3 -O0

INC = -I./include -I./include/matrix_math

# all: dir $(BUILDDIR)/$(EXECUTABLE)

# dir:
# 	mkdir -p $(BUILDDIR)

# $(BUILDDIR)/$(EXECUTABLE): $(OBJECTS)
# 	$(CC) $^ -o $@ $(INC) $(LIBDIR) $(LIBS) 

# $(OBJECTS): $(BUILDDIR)/%.o : $(SOURCEDIR)/%.c
# 	$(CC) $(FLAGS) $< -o $@ $(INC) 

# clean:
# 	rm -f $(BUILDDIR)/*o $(BUILDDIR)/$(EXECUTABLE)


all: ADCS.c src/matrix_math/matrix_math.c
	$(CC) $(FLAGS) -o ADCS ADCS.c src/matrix_math/matrix_math.c $(INC) $(LIBDIR) $(LIBS) 


clean:
	rm *.exe