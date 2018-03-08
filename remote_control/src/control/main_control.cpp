/* main_control.cpp
 * this file is part of Dynstack/RemoteControl for CORSIKA
 *
 * Copyright (C) <2016> <Dominik Baack>
 *		All rights reserved.
 *
 * 	This software may be modified and distributed under the terms
 * 	of the LGPL license. See the LICENSE file for details.
 */

#include "remote_control/control/main_control.h"

#include "io/network/dns_lookup.h"

namespace remote_control
{

	void MainControl::send_loop()
	{
		while (m_tsRunning == true)
		{
			auto msg_send = m_tsQueue.pop_front(std::chrono::milliseconds(5));
			if (msg_send.empty())
				continue;
			
			if (!m_client.send(msg_send.toByte(), static_cast<const unsigned int>(msg_send.size()) ))
			{
				std::cerr << __FILE__ << " Message could not be sent" << std::endl;
			}
		}
	}


	void MainControl::recv_loop()
	{
		std::cout << "(RC) Start main loop..." << std::endl;
		//auto start_time = std::chrono::high_resolution_clock::now();

		while (m_tsRunning == true)
		{
			const auto msg_recv = m_client.recv(std::chrono::milliseconds(5));

			if (!msg_recv.empty())
			{
				//evaluate message
				remote_control::communication::Packet recv_packet(msg_recv);

				// Check for packet consistency, returns 0 if everything is fine else an error code
				// C++17 allows if with initializer
				int check = recv_packet.check();
				if( check )
				{
					std::cerr << "Received a compromised package!" << std::endl;
				}
				else
				{
					// Send to callback functions
					auto callback = m_callback.find( recv_packet.header() );
					if(callback != m_callback.end())
					{
						callback->second( recv_packet.data() );
					}
				}

			}

			/// Call periodic functions
			auto current_time = std::chrono::high_resolution_clock::now();

			for(auto itr : m_periodic)
			{
				itr.call( current_time );
			}

		}

	}

	MainControl::MainControl()
			: m_periodic(register_periodic_callback()), m_callback(register_server_callback()),
			m_tsQueue(1024)
	{
		m_tsRunning = false;
	}

	MainControl::~MainControl()
	{

	}

	bool MainControl::start(const std::string dns, const unsigned short port)
	{
		if (m_tsRunning == true)
		{
			return false;
		}

		const std::string ip = io::network::hostname_to_ip(dns);
		if (!m_client.init(ip, port))
		{
			std::cerr << "Could not connect to " << dns << ":" << port << "(" << ip << ")" << std::endl;
			return false;
		}

		m_tsRunning = true;
		m_threads.clear();
		{
			std::packaged_task<void()> task(std::bind(&MainControl::recv_loop, this));
			std::future<void> future = task.get_future();
			m_threads.emplace_back(std::move(std::thread(std::move(task))), std::move(future));
		}
		{
			std::packaged_task<void()> task(std::bind(&MainControl::send_loop, this));
			std::future<void> future = task.get_future();
			m_threads.emplace_back(std::move(std::thread(std::move(task))), std::move(future));
		}

		return true;

	}

	void MainControl::stopp()
	{
		if (m_tsRunning == true)
		{
			while (!m_tsQueue.empty()) {
				std::cout << "(RC) waiting for queue to empty..." << std::endl;
				std::this_thread::sleep_for(std::chrono::milliseconds(100));
			}
			m_tsRunning = false;
			for (auto &thread : m_threads)
				thread.first.join();
			
			/// Terminate all connections
			std::cout << "(RC) Close main loop ..." << std::endl;
			m_client.close();
			
			// Raise exception, if any, on the main thread
			for (auto &thread : m_threads)
				thread.second.get();


			std::cout << "(RC) Main loop closed" << std::endl;
		}
	}

	void MainControl::send(remote_control::communication::Packet p)
	{
		if (m_tsRunning)
			this->m_tsQueue.push_back(p);
	}

} /* namespace RemoteControl */
