CC = g++ -O3 -DNDEBUG -std=gnu++0x

ifeq ($(TOOLSET), clang)
  CC = clang++ -O3 -DNDEBUG -std=gnu++11
else ifeq ($(TOOLSET), intel-linux)
  CC = icpc -O3 -DNDEBUG -std=gnu++0x
endif

bin/idock: obj/scoring_function.o obj/box.o obj/quaternion.o obj/thread_pool.o obj/receptor_parser.o obj/receptor.o obj/ligand_parser.o obj/ligand.o obj/grid_map_task.o obj/monte_carlo_task.o obj/main.o
	$(CC) -o $@ $^ -pthread -lboost_system -lboost_thread -lboost_filesystem -lboost_program_options

obj/%.o: src/%.cpp 
	$(CC) -o $@ $< -c

clean:
	rm -f bin/idock obj/*.o
