/*

   Copyright (c) 2011, The Chinese University of Hong Kong

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

*/

#include "thread_pool.hpp"

namespace idock
{
	thread_pool::thread_pool(const size_t num_threads) : num_threads(num_threads), tasks_ptr(nullptr), hashes_ptr(nullptr), num_tasks(0), num_started_tasks(0), num_completed_tasks(0), next_hash_index(0), exiting(false)
	{
		// Create threads to call (*this)().
		for (size_t i = 0; i < num_threads; ++i)
			create_thread(boost::ref(*this));
	}

	void thread_pool::run(ptr_vector<packaged_task<void> >& tasks, array<fl, num_hashes>& hashes)
	{
		// Initialize several task counters for scheduling.
		tasks_ptr = &tasks;
		num_tasks = tasks.size();
		num_started_tasks = 0;
		num_completed_tasks = 0;

		// Initialize hash variables for displaying progress bar.
		hashes_ptr = &hashes;
		next_hash_index = 0;

		// Notify the threads to run tasks.
		task_incoming.notify_all();

		// Wait until any thread becomes available.
		{
			mutex::scoped_lock self_lk(self); // A scoped lock is a type associated to some mutex type whose objects do the locking/unlocking of a mutex on construction/destruction time.
			while (num_completed_tasks + (num_threads - 1) < num_tasks)
				// wait will atomically add the thread to the set of threads waiting on the condition variable, and unlock the mutex.
				// When the thread is woken, the mutex will be locked again before the call to wait returns.
				// This allows other threads to acquire the mutex in order to update the shared data,
				// and ensures that the data associated with the condition is correctly synchronized.
				task_completion.wait(self_lk); // Wait for tasks completion.
		}
	}

	void thread_pool::operator()()
	{
		packaged_task<void>* task; // Declare a pointer to a task.
		do // Threads loop inside.
		{
			// If there are no tasks to run, simply wait.
			{				
				mutex::scoped_lock self_lk(self); // A scoped lock is a type associated to some mutex type whose objects do the locking/unlocking of a mutex on construction/destruction time.
				while ((!exiting) && (num_started_tasks == num_tasks))
					task_incoming.wait(self_lk);
				if (exiting) return; // The only exit of this function.
			}
			
			// If there are tasks to run, loop until all the tasks are started.
			do
			{
				// Fetch a task to run atomically.
				// TODO: Consider using boost/atomic.hpp, i.e. atomic<size_t> num_started_tasks; num_started_tasks.fetch_add(1);
				{					
					mutex::scoped_lock self_lk(self);
					if (num_started_tasks == num_tasks) break; // Break the loop when all the tasks are distributed.
					task = &(*tasks_ptr)[num_started_tasks++];
				}

				// Run the task. Potential exceptions are handled by the caller of thread_pool::run.
				task->operator()();

				// Increment the number of completed tasks atomically.
				{
					mutex::scoped_lock self_lk(self);
					++num_completed_tasks;

					// Flush hashes if necessary.
					while ((next_hash_index < num_hashes) && (num_completed_tasks >= (*hashes_ptr)[next_hash_index]))
					{
						std::cout << '#' << std::flush;
						++next_hash_index;
					}
				}

				// One task is completed. Notify the main thread.
				task_completion.notify_one();
			} while (true);
		} while (true);
	}

	void thread_pool::sync()
	{
		mutex::scoped_lock self_lk(self);
		while (num_completed_tasks < num_tasks) // Wait until all the tasks are completed.
			task_completion.wait(self_lk);
	}

	thread_pool::~thread_pool()
	{
		// Notify threads to exit from the loop back function.
		exiting = true;
		task_incoming.notify_all();

		// Wait until all threads are joined.
		join_all();
	}
}
