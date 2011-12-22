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

#ifndef IDOCK_THREAD_POOL_HPP
#define IDOCK_THREAD_POOL_HPP

#include <iostream>
#include <boost/array.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/condition.hpp>
#include <boost/thread/future.hpp>
#include "common.hpp"

namespace idock
{
	using boost::array;
	using boost::ptr_vector;
	using boost::thread_group;
	using boost::mutex;
	using boost::condition;
	using boost::packaged_task;
	using boost::unique_future;

	const size_t num_hashes = 10; ///< Number of hashes in a progress bar.

	/// Represents a thread pool and incorporates a progress bar. It inherits from boost::thread_group for the usage of create_thread and join_all.
	class thread_pool : public thread_group
	{
	public:

		/// Constructs a thread pool with specified number of threads.
		explicit thread_pool(const size_t num_threads);

		/// Runs tasks in parallel asynchronously.
		void run(ptr_vector<packaged_task<void> >& tasks, array<fl, num_hashes>& hashes);

		/// The function for threads to execute and loop inside.
		void operator()();

		/// Blocks until any thread becomes available.
		void soft_sync();

		/// Blocks until all tasks are completed and the progress bar becomes full.
		void hard_sync();

		/// Joins all the threads.
		void dispose();

		/// Destructs a thread pool by joining all the threads.
		~thread_pool();

	protected:
		const size_t num_threads; ///< Number of threads to run tasks.
		ptr_vector<packaged_task<void> >* tasks_ptr; ///< Pointer to the tasks to run.
		array<fl, num_hashes>* hashes_ptr; ///< Pointer to the hashes.
		size_t num_tasks; ///< Number of tasks.
		size_t num_started_tasks; ///< Number of tasks that have started running.
		size_t num_completed_tasks; ///< Number of tasks that have completed running.
		size_t next_hash_index; ///< Index to the next hash value.
		condition task_completion;	///< Completion event of a running task.
		condition task_incoming; ///< Incoming event of new tasks.
		bool exiting; ///< If true, notify threads to return.
		mutex self; ///< A mutex is a type whose objects can be in either of two states, called locked and unlocked, with the property that when a thread A has locked a mutex m and a different thread B tries to lock m, B is blocked until A unlocks m.
	};
}

#endif
