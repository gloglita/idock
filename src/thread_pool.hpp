#pragma once
#ifndef IDOCK_THREAD_POOL_HPP
#define IDOCK_THREAD_POOL_HPP

#include <vector>
#include <future>
using namespace std;

class thread_pool : public vector<packaged_task<int()>>
{
public:
	explicit thread_pool(const size_t num_threads);
	void sync(const size_t num_bars = 0);
	~thread_pool();
private:
	vector<thread> threads;
	size_t num_scheduled_tasks;
	size_t num_completed_tasks;
	condition_variable task_completion;
	condition_variable task_incoming;
	bool exiting;
	mutable mutex m;
	thread_pool(thread_pool const&);
	thread_pool& operator=(thread_pool const&);
};

#endif
