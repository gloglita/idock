#include <iostream>
#include "thread_pool.hpp"

thread_pool::thread_pool(const size_t num_threads) : num_scheduled_tasks(0), num_completed_tasks(0), exiting(false)
{
	for (size_t i = 0; i < num_threads; ++i)
	{
		threads.push_back(thread([&]()
		{
			while (true)
			{
				{
					unique_lock<mutex> lock(m);
					if (num_scheduled_tasks == size())
					{
						task_incoming.wait(lock);
					}
				}
				if (exiting) return;
				while (true)
				{
					size_t t;
					{
						unique_lock<mutex> lock(m);
						if (num_scheduled_tasks == size()) break;
						t = num_scheduled_tasks++;
					}
					(*this)[t].operator()();
					{
						unique_lock<mutex> lock(m);
						++num_completed_tasks;
					}
					task_completion.notify_one();
				}
			}
		}));
	}
}

void thread_pool::sync(const size_t num_bars)
{
	task_incoming.notify_all();
	unique_lock<mutex> lock(m);
	for (size_t i = 0; num_completed_tasks < size();)
	{
		task_completion.wait(lock);
		for (; num_completed_tasks * num_bars >= size() * (i + 1); ++i)
		{
			cout << '#' << flush;
		}
	}
	for (auto& t : *this)
	{
		t.get_future().get();
	}
	clear();
	num_scheduled_tasks = 0;
	num_completed_tasks = 0;
}

thread_pool::~thread_pool()
{
	exiting = true;
	task_incoming.notify_all();
	for (auto& t : threads)
	{
		t.join();
	}
}
