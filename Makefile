BOOST_ROOT = $(HOME)/boost_1_47_0
CC = g++ -static -O3 -DNDEBUG -std=gnu++0x

ifeq ($(TOOLSET), clang)
  CC = clang++ -static -O3 -DNDEBUG -Wno-c++0x-extensions
else ifeq ($(TOOLSET), intel-linux)
  CC = icpc -static -O3 -DNDEBUG -std=gnu++0x
endif

obj/%.o : src/%.cpp 
	$(CC) -o $@ -c $< -I$(BOOST_ROOT) 

bin/idock: obj/scoring_function.o obj/box.o obj/quaternion.o obj/thread_pool.o obj/receptor_parser.o obj/receptor.o obj/ligand_parser.o obj/ligand.o obj/grid_map_task.o obj/monte_carlo_task.o obj/main.o
	$(CC) -o $@ $^ -L$(BOOST_ROOT)/lib -pthread -lboost_system -lboost_thread -lboost_filesystem -lboost_program_options

clean:
	rm -f bin/idock obj/*.o
