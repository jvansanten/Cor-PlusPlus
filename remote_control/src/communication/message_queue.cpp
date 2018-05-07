/* message_queue.cpp
 * this file is part of Dynstack/RemoteControl for CORSIKA
 *
 * Copyright (C) <2016> <Dominik Baack>
 *		All rights reserved.
 *
 * 	This software may be modified and distributed under the terms
 * 	of the LGPL license. See the LICENSE file for details.
 */


#include "remote_control/communication/message_queue.h"

#include <iostream>

namespace remote_control
{
	namespace communication
	{

		MessageQueue::MessageQueue()
		{
			m_size = 0;
			m_capacity = std::numeric_limits<decltype(m_capacity)>::max();
		}

		MessageQueue::MessageQueue(uint32_t capacity)
		{
			m_size = 0;
			m_capacity = capacity;
		}

		MessageQueue::~MessageQueue()
		{

		}


		void MessageQueue::push_back(MessageQueue::storage_type data)
		{
			std::unique_lock<std::mutex> lock(m_mutex);
			while (m_size.load() >= m_capacity)
			{
				m_cv.wait(lock);
			}
			m_size++;
			m_queue.push_back(data);

			m_cv.notify_one();
		}
		void MessageQueue::push_front(MessageQueue::storage_type data)
		{
			std::unique_lock<std::mutex> lock(m_mutex);
			while (m_size.load() >= m_capacity)
			{
				m_cv.wait(lock);
			}
			m_size++;
			m_queue.push_front(data);

			m_cv.notify_one();
		}		

		MessageQueue::storage_type MessageQueue::pop_back_blocking()
		{
			std::unique_lock<std::mutex> lock(m_mutex);
			while (m_size.load() <= 0)
			{
				m_cv.wait(lock);
			}			

			auto tmp = std::move(m_queue.back());
			m_queue.pop_back();
			
			m_size--;
			m_cv.notify_one();

			return tmp;
		}

		MessageQueue::storage_type MessageQueue::pop_front_blocking()
		{
			std::unique_lock<std::mutex> lock(m_mutex);
			while (m_size.load() <= 0)
			{
				m_cv.wait(lock);
			}

			storage_type tmp = std::move(m_queue.front());
			m_queue.pop_front();	

			m_size--;
			m_cv.notify_one();

			return tmp;
		}
		
		MessageQueue::storage_type MessageQueue::pop_front(std::chrono::microseconds timeout)
		{
			std::unique_lock<std::mutex> lock(m_mutex);
			if (m_cv.wait_for(lock, timeout, [this]{ return this->m_size.load() > 0; })) {
				storage_type tmp = std::move(m_queue.front());
				m_queue.pop_front();

				m_size--;
				m_cv.notify_one();
				
				return tmp;
			} else {
				return MessageQueue::storage_type();
			}

		}

		bool MessageQueue::empty()
		{
			std::unique_lock<std::mutex> lock(m_mutex);
			return m_queue.empty();
		}

	}
} /* namespace lib */
