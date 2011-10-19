CC = g++ -static -O3 -DNDEBUG -std=gnu++0x
SRC = src
OBJ = obj
BIN = bin

.SUFFIXES: .cpp .o

$(OBJ)/%.o : $(SRC)/%.cpp 
	$(CC) -o $@ -c $< 

$(BIN)/idock: $(OBJ)/scoring_function.o $(OBJ)/box.o $(OBJ)/quaternion.o $(OBJ)/thread_pool.o $(OBJ)/receptor_parser.o $(OBJ)/receptor.o $(OBJ)/ligand_parser.o $(OBJ)/ligand.o $(OBJ)/grid_map_task.o $(OBJ)/monte_carlo_task.o $(OBJ)/main.o
	$(CC) -o $@ $^ -pthread -l boost_system -l boost_thread -l boost_filesystem -l boost_program_options

clean:
	rm -f $(BIN)/idock $(OBJ)/*.o
