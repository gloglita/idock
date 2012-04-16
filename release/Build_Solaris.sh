rm -f ../bin/Solaris/x86_64/idock ../obj/Solaris/x86_64/*.o
g++ -m64 -O3 -DNDEBUG -std=gnu++0x -pthreads -o ../obj/Solaris/x86_64/scoring_function.o ../src/scoring_function.cpp -I/home/hjli/boost_1_49_0 -c
g++ -m64 -O3 -DNDEBUG -std=gnu++0x -pthreads -o ../obj/Solaris/x86_64/box.o ../src/box.cpp -I/home/hjli/boost_1_49_0 -c
g++ -m64 -O3 -DNDEBUG -std=gnu++0x -pthreads -o ../obj/Solaris/x86_64/quaternion.o ../src/quaternion.cpp -I/home/hjli/boost_1_49_0 -c
g++ -m64 -O3 -DNDEBUG -std=gnu++0x -pthreads -o ../obj/Solaris/x86_64/thread_pool.o ../src/thread_pool.cpp -I/home/hjli/boost_1_49_0 -c
g++ -m64 -O3 -DNDEBUG -std=gnu++0x -pthreads -o ../obj/Solaris/x86_64/receptor.o ../src/receptor.cpp -I/home/hjli/boost_1_49_0 -c
g++ -m64 -O3 -DNDEBUG -std=gnu++0x -pthreads -o ../obj/Solaris/x86_64/ligand.o ../src/ligand.cpp -I/home/hjli/boost_1_49_0 -c
g++ -m64 -O3 -DNDEBUG -std=gnu++0x -pthreads -o ../obj/Solaris/x86_64/grid_map_task.o ../src/grid_map_task.cpp -I/home/hjli/boost_1_49_0 -c
g++ -m64 -O3 -DNDEBUG -std=gnu++0x -pthreads -o ../obj/Solaris/x86_64/monte_carlo_task.o ../src/monte_carlo_task.cpp -I/home/hjli/boost_1_49_0 -c
g++ -m64 -O3 -DNDEBUG -std=gnu++0x -pthreads -o ../obj/Solaris/x86_64/main.o ../src/main.cpp -I/home/hjli/boost_1_49_0 -c
g++ -m64 -O3 -DNDEBUG -o ../bin/Solaris/x86_64/idock ../obj/Solaris/x86_64/scoring_function.o ../obj/Solaris/x86_64/box.o ../obj/Solaris/x86_64/quaternion.o ../obj/Solaris/x86_64/thread_pool.o ../obj/Solaris/x86_64/receptor.o ../obj/Solaris/x86_64/ligand.o ../obj/Solaris/x86_64/grid_map_task.o ../obj/Solaris/x86_64/monte_carlo_task.o ../obj/Solaris/x86_64/main.o -L/home/hjli/boost_1_49_0/lib/x86_64 -lboost_system -lboost_thread -lboost_filesystem -lboost_program_options
rm -f ../bin/Solaris/x86/idock ../obj/Solaris/x86/*.o
g++ -m32 -O3 -DNDEBUG -std=gnu++0x -pthreads -o ../obj/Solaris/x86/scoring_function.o ../src/scoring_function.cpp -I/home/hjli/boost_1_49_0 -c
g++ -m32 -O3 -DNDEBUG -std=gnu++0x -pthreads -o ../obj/Solaris/x86/box.o ../src/box.cpp -I/home/hjli/boost_1_49_0 -c
g++ -m32 -O3 -DNDEBUG -std=gnu++0x -pthreads -o ../obj/Solaris/x86/quaternion.o ../src/quaternion.cpp -I/home/hjli/boost_1_49_0 -c
g++ -m32 -O3 -DNDEBUG -std=gnu++0x -pthreads -o ../obj/Solaris/x86/thread_pool.o ../src/thread_pool.cpp -I/home/hjli/boost_1_49_0 -c
g++ -m32 -O3 -DNDEBUG -std=gnu++0x -pthreads -o ../obj/Solaris/x86/receptor.o ../src/receptor.cpp -I/home/hjli/boost_1_49_0 -c
g++ -m32 -O3 -DNDEBUG -std=gnu++0x -pthreads -o ../obj/Solaris/x86/ligand.o ../src/ligand.cpp -I/home/hjli/boost_1_49_0 -c
g++ -m32 -O3 -DNDEBUG -std=gnu++0x -pthreads -o ../obj/Solaris/x86/grid_map_task.o ../src/grid_map_task.cpp -I/home/hjli/boost_1_49_0 -c
g++ -m32 -O3 -DNDEBUG -std=gnu++0x -pthreads -o ../obj/Solaris/x86/monte_carlo_task.o ../src/monte_carlo_task.cpp -I/home/hjli/boost_1_49_0 -c
g++ -m32 -O3 -DNDEBUG -std=gnu++0x -pthreads -o ../obj/Solaris/x86/main.o ../src/main.cpp -I/home/hjli/boost_1_49_0 -c
g++ -m32 -O3 -DNDEBUG -o ../bin/Solaris/x86/idock ../obj/Solaris/x86/scoring_function.o ../obj/Solaris/x86/box.o ../obj/Solaris/x86/quaternion.o ../obj/Solaris/x86/thread_pool.o ../obj/Solaris/x86/receptor.o ../obj/Solaris/x86/ligand.o ../obj/Solaris/x86/grid_map_task.o ../obj/Solaris/x86/monte_carlo_task.o ../obj/Solaris/x86/main.o -L/home/hjli/boost_1_49_0/lib/x86 -lboost_system -lboost_thread -lboost_filesystem -lboost_program_options
