BOOST_ROOT = $(HOME)/boost_1_51_0
CC = g++ -O3 -DNDEBUG -std=gnu++0x

ifeq ($(TOOLSET), clang)
  CC = clang++ -O3 -DNDEBUG -std=gnu++11
else ifeq ($(TOOLSET), intel)
  CC = icpc -O3 -DNDEBUG -std=gnu++0x
endif

bin/idock: obj/scoring_function.o obj/box.o obj/quaternion.o obj/thread_pool.o obj/receptor.o obj/ligand.o obj/grid_map_task.o obj/monte_carlo_task.o obj/main.o
	$(CC) -o $@ $^ -L$(BOOST_ROOT)/lib/x86_64 -pthread -lboost_system -lboost_thread -lboost_filesystem -lboost_program_options -lboost_iostreams -lz -lbz2

obj/%.o: src/%.cpp 
	$(CC) -o $@ $< -I$(BOOST_ROOT) -c

clean:
	rm -f bin/idock obj/*.o
