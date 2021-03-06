BOOST_ROOT = $(HOME)/boost_1_53_0
CC = clang++ -O3 -DNDEBUG -std=c++11

bin/idock: obj/scoring_function.o obj/thread_pool.o obj/receptor.o obj/ligand.o obj/main.o
	$(CC) -o $@ $^ -L$(BOOST_ROOT)/lib/x86_64 -pthread -lboost_system -lboost_program_options -lboost_filesystem

obj/%.o: src/%.cpp 
	$(CC) -o $@ $< -I$(BOOST_ROOT) -c

clean:
	rm -f bin/idock obj/*.o
