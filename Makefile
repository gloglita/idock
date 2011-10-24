CC = g++ -static -O3 -DNDEBUG -std=gnu++0x

ifeq ($(TOOLSET), clang)
  CC = clang++ -static -O3 -DNDEBUG -Wno-c++0x-extensions
else ifeq ($(TOOLSET), intel-linux)
  CC = icpc -static -O3 -DNDEBUG -std=gnu++0x
endif

SRC = src
OBJ = obj
BIN = bin

$(OBJ)/%.o : $(SRC)/%.cpp 
	$(CC) -o $@ -c $< 

$(BIN)/idock: $(OBJ)/scoring_function.o $(OBJ)/box.o $(OBJ)/quaternion.o $(OBJ)/thread_pool.o $(OBJ)/receptor_parser.o $(OBJ)/receptor.o $(OBJ)/ligand_parser.o $(OBJ)/ligand.o $(OBJ)/grid_map_task.o $(OBJ)/monte_carlo_task.o $(OBJ)/main.o
	$(CC) -o $@ $^ -pthread -lboost_system -lboost_thread -lboost_filesystem -lboost_program_options

clean:
	rm -f $(BIN)/idock $(OBJ)/*.o
