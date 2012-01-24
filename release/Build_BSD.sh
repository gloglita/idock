rm -f bin/BSD/x86_64/igrow obj/BSD/x86_64/*.o
clang++ -m64 -static -O3 -DNDEBUG -std=gnu++11 -o ../obj/BSD/x86_64/scoring_function.o ../src/scoring_function.cpp -I/home/hjli/boost_1_48_0 -c
clang++ -m64 -static -O3 -DNDEBUG -std=gnu++11 -o ../obj/BSD/x86_64/box.o ../src/box.cpp -I/home/hjli/boost_1_48_0 -c
clang++ -m64 -static -O3 -DNDEBUG -std=gnu++11 -o ../obj/BSD/x86_64/quaternion.o ../src/quaternion.cpp -I/home/hjli/boost_1_48_0 -c
clang++ -m64 -static -O3 -DNDEBUG -std=gnu++11 -o ../obj/BSD/x86_64/thread_pool.o ../src/thread_pool.cpp -I/home/hjli/boost_1_48_0 -c
clang++ -m64 -static -O3 -DNDEBUG -std=gnu++11 -o ../obj/BSD/x86_64/receptor.o ../src/receptor.cpp -I/home/hjli/boost_1_48_0 -c
clang++ -m64 -static -O3 -DNDEBUG -std=gnu++11 -o ../obj/BSD/x86_64/ligand.o ../src/ligand.cpp -I/home/hjli/boost_1_48_0 -c
clang++ -m64 -static -O3 -DNDEBUG -std=gnu++11 -o ../obj/BSD/x86_64/grid_map_task.o ../src/grid_map_task.cpp -I/home/hjli/boost_1_48_0 -c
clang++ -m64 -static -O3 -DNDEBUG -std=gnu++11 -o ../obj/BSD/x86_64/monte_carlo_task.o ../src/monte_carlo_task.cpp -I/home/hjli/boost_1_48_0 -c
clang++ -m64 -static -O3 -DNDEBUG -std=gnu++11 -o ../obj/BSD/x86_64/main.o ../src/main.cpp -I/home/hjli/boost_1_48_0 -c
clang++ -m64 -static -O3 -DNDEBUG -std=gnu++11 -o ../bin/BSD/x86_64/ ../obj/BSD/x86_64/scoring_function.o ../obj/BSD/x86_64/box.o ../obj/BSD/x86_64/quaternion.o ../obj/BSD/x86_64/thread_pool.o ../obj/BSD/x86_64/receptor.o ../obj/BSD/x86_64/ligand.o ../obj/BSD/x86_64/grid_map_task.o ../obj/BSD/x86_64/monte_carlo_task.o ../obj/BSD/x86_64/main.o -L/home/hjli/boost_1_48_0/lib/x86_64 -pthread -lboost_system -lboost_thread -lboost_filesystem -lboost_program_options
rm -f bin/BSD/x86/igrow obj/BSD/x86/*.o
clang++ -m32 -static -O3 -DNDEBUG -std=gnu++11 -o ../obj/BSD/x86/scoring_function.o ../src/scoring_function.cpp -I/home/hjli/boost_1_48_0 -c
clang++ -m32 -static -O3 -DNDEBUG -std=gnu++11 -o ../obj/BSD/x86/box.o ../src/box.cpp -I/home/hjli/boost_1_48_0 -c
clang++ -m32 -static -O3 -DNDEBUG -std=gnu++11 -o ../obj/BSD/x86/quaternion.o ../src/quaternion.cpp -I/home/hjli/boost_1_48_0 -c
clang++ -m32 -static -O3 -DNDEBUG -std=gnu++11 -o ../obj/BSD/x86/thread_pool.o ../src/thread_pool.cpp -I/home/hjli/boost_1_48_0 -c
clang++ -m32 -static -O3 -DNDEBUG -std=gnu++11 -o ../obj/BSD/x86/receptor.o ../src/receptor.cpp -I/home/hjli/boost_1_48_0 -c
clang++ -m32 -static -O3 -DNDEBUG -std=gnu++11 -o ../obj/BSD/x86/ligand.o ../src/ligand.cpp -I/home/hjli/boost_1_48_0 -c
clang++ -m32 -static -O3 -DNDEBUG -std=gnu++11 -o ../obj/BSD/x86/grid_map_task.o ../src/grid_map_task.cpp -I/home/hjli/boost_1_48_0 -c
clang++ -m32 -static -O3 -DNDEBUG -std=gnu++11 -o ../obj/BSD/x86/monte_carlo_task.o ../src/monte_carlo_task.cpp -I/home/hjli/boost_1_48_0 -c
clang++ -m32 -static -O3 -DNDEBUG -std=gnu++11 -o ../obj/BSD/x86/main.o ../src/main.cpp -I/home/hjli/boost_1_48_0 -c
clang++ -m32 -static -O3 -DNDEBUG -std=gnu++11 -o ../bin/BSD/x86/ ../obj/BSD/x86/scoring_function.o ../obj/BSD/x86/box.o ../obj/BSD/x86/quaternion.o ../obj/BSD/x86/thread_pool.o ../obj/BSD/x86/receptor.o ../obj/BSD/x86/ligand.o ../obj/BSD/x86/grid_map_task.o ../obj/BSD/x86/monte_carlo_task.o ../obj/BSD/x86/main.o -L/home/hjli/boost_1_48_0/lib/x86 -pthread -lboost_system -lboost_thread -lboost_filesystem -lboost_program_options
