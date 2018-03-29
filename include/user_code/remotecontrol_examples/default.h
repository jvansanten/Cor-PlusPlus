/* default.h
 * this file is part of RemoteControl for CORSIKA
 *
 * Copyright (C) <2016> <Dominik Baack>
 *		All rights reserved.
 *
 * 	This software may be modified and distributed under the terms
 * 	of the LGPL license. See the LICENSE file for details.
 */



/// Contains default handling for basic message transfer and data logging.
/**
 *  This functions implement all messages that are needed to use the public remote server for logging purposes
 */



#pragma once

#include <chrono>

#include "remote_control/control/periodic_task.h"
#include "remote_control/communication/packet.h"

remote_control::communication::Packet
recv_primary(std::vector<uint8_t> msg);
