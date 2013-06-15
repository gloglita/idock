rm -f ../bin/FreeBSD/x86/idock ../obj/FreeBSD/x86/*.o
clang++ -m32 -static -O3 -DNDEBUG -std=c++11 -o ../obj/FreeBSD/x86/scoring_function.o ../src/scoring_function.cpp -I/home/hjli/boost_1_51_0 -c
clang++ -m32 -static -O3 -DNDEBUG -std=c++11 -o ../obj/FreeBSD/x86/box.o ../src/box.cpp -I/home/hjli/boost_1_51_0 -c
clang++ -m32 -static -O3 -DNDEBUG -std=c++11 -o ../obj/FreeBSD/x86/quaternion.o ../src/quaternion.cpp -I/home/hjli/boost_1_51_0 -c
clang++ -m32 -static -O3 -DNDEBUG -std=c++11 -o ../obj/FreeBSD/x86/thread_pool.o ../src/thread_pool.cpp -I/home/hjli/boost_1_51_0 -c
clang++ -m32 -static -O3 -DNDEBUG -std=c++11 -o ../obj/FreeBSD/x86/receptor.o ../src/receptor.cpp -I/home/hjli/boost_1_51_0 -c
clang++ -m32 -static -O3 -DNDEBUG -std=c++11 -o ../obj/FreeBSD/x86/ligand.o ../src/ligand.cpp -I/home/hjli/boost_1_51_0 -c
clang++ -m32 -static -O3 -DNDEBUG -std=c++11 -o ../obj/FreeBSD/x86/grid_map_task.o ../src/grid_map_task.cpp -I/home/hjli/boost_1_51_0 -c
clang++ -m32 -static -O3 -DNDEBUG -std=c++11 -o ../obj/FreeBSD/x86/monte_carlo_task.o ../src/monte_carlo_task.cpp -I/home/hjli/boost_1_51_0 -c
clang++ -m32 -static -O3 -DNDEBUG -std=c++11 -o ../obj/FreeBSD/x86/main.o ../src/main.cpp -I/home/hjli/boost_1_51_0 -c
clang++ -m32 -static -O3 -DNDEBUG -o ../bin/FreeBSD/x86/idock ../obj/FreeBSD/x86/scoring_function.o ../obj/FreeBSD/x86/box.o ../obj/FreeBSD/x86/quaternion.o ../obj/FreeBSD/x86/thread_pool.o ../obj/FreeBSD/x86/receptor.o ../obj/FreeBSD/x86/ligand.o ../obj/FreeBSD/x86/grid_map_task.o ../obj/FreeBSD/x86/monte_carlo_task.o ../obj/FreeBSD/x86/main.o -L/home/hjli/boost_1_51_0/lib/x86 -pthread -lboost_system -lboost_thread -lboost_filesystem -lboost_program_options -lboost_iostreams -lz -lbz2
