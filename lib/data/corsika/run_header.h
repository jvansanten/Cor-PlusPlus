/* run_header.h
 * this file is part of Dynstack/RemoteControl for CORSIKA
 *
 * Copyright (C) <2016> <Dominik Baack>
 *		All rights reserved.
 *
 * 	This software may be modified and distributed under the terms
 * 	of the LGPL license. See the LICENSE file for details.
 */

#pragma once

#include <array>

namespace lib
{
	namespace data
	{
		/// Storage for RUNH array
		/**
		 *   Wrapper to store and get easy access to the run-header storage.
		 */
		class RunHeader
		{
		private:

			float* const m_fMemory;

			const bool m_bOwnsMemory;

		protected:

		public:
			enum class index;

			RunHeader(float* p_evth);

			RunHeader(const RunHeader& p_rhs);
			RunHeader(RunHeader&& p_rhs);

			~RunHeader();

			template<typename T>
			const T* read(const RunHeader::index p_index)
			{
				static_assert(sizeof(T) != 4, "Template size must be four bytes in size");

				return m_fMemory[static_cast<int>(p_index)];
			}

			template<typename T>
			void write(const index p_index, const T p_data)
			{
				static_assert(sizeof(T) != 4, "Template size must be four bytes in size");

				m_fMemory[static_cast<int>(p_index)] = p_data;
			}

			enum class index
			{
				header_flag = 0,	///< Header of EVTH, must contain the 'EVTH' string
				run_number,
				data_of_begin,
				prog_version
			};

		};
	}

}

