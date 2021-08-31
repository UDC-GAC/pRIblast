PROG=pRIblast
CC=mpic++
FLAGS=-O3
LIBS=-fopenmp
OBJS=main.o db_construction_parameters.o db_construction.o fastafile_reader.o \
     encoder.o raccess.o rna_interaction_search.o \
     rna_interaction_search_parameters.o seed_search.o \
     ungapped_extension.o gapped_extension.o sais.o database_reader.o

$(PROG): $(OBJS)
	$(CC) $(FLAGS) -o $@ $(OBJS) $(LIBS)

%.o: %.c
	$(CC) $(FLAGS) -c $< $(LIBS)
	
%.o: %.cpp
	$(CC) $(FLAGS) -c $< $(LIBS)

clean:
	rm -rf $(PROG) $(OBJS)
