from socket import *
import time
import sys
import hashlib
import select

#### TCP IMPLEMENTATION HERE ####

serverName = sys.argv[1] #SERVER NAME
serverPort = int(sys.argv[3]) #TCP Client Port

packet_size = 966 # Packet size chosen to make all data as 1000bytes long
tc_socket = socket(AF_INET, SOCK_STREAM) #Create Socket
tc_socket.connect((serverName,serverPort)) #Connect to the server
#Read all the file into a list (buffer)
bfr = []
with open("transfer_file_TCP.txt",'r') as f:
	while True:
		pckt = f.read(packet_size)
		if not pckt:
			break
		bfr.append(pckt)
bfr.append('') #Append empty string to the end
all_msg = ''
count = 0
while True: # Send file as chunks
	chunk = bfr[count] 
	count += 1
	if len(chunk) == 0:#If all buffer data is sent, break
		break
	
	timestamp = str(time.time()) #timestamp
	if len(timestamp)>17: # Fixed length for timestamp
		timestamp = timestamp[:17]
	msg = "<time>"+timestamp+"<msg>"+chunk+"</end>" #Add necessary tags to the data
	if count == len(bfr)-1: #If this is last data, mark it with </ed> token
		msg += "</ed>"
	msg += "x"*(1000-len(msg))
	tc_socket.send(msg.encode()) #Send encoded message via socket

tc_socket.close() #Close the socket

#### UDP IMPLEMENTATION HERE ####
#Same initializations here as in TCP
serverName = sys.argv[1]
serverPort = int(sys.argv[2])
packet_size = 966
uc_socket = socket(AF_INET, SOCK_DGRAM)

uc_socket.bind(('',int(sys.argv[4])))
uc_socket.settimeout(0.3) #Set the timeout for the socket in order to handle corruptions
f = open("transfer_file_UDP.txt",'r')
count = 0
r = 0
last_seq = 0
while True: # Send file as chunks
	chunk = f.read(packet_size)
	count += 1

	if len(chunk) == 0: #If empty, break
		f.close()
		break
	
	checksum = hashlib.md5(chunk.encode()).hexdigest() #Calculate chekcsum value of the chunk
	timestamp = str(time.time())
	msg = "<time>"+timestamp+"<cs>"+checksum+"<seq>"+str(count)+"<msg>"+chunk # This is same but as in tcp but this time it has checksum and sequence number for comparison
	if len(chunk) < packet_size: # If this is the last packet, add <end> token at the end
		msg += "<end>"
	while True:
		
		try: # Here I try to send my chunk and get a ACK receive from the server
			msg = "<time>"+str(time.time())+msg[msg.find("<cs>"):] #Update timestamp in the packet
			uc_socket.sendto(msg.encode(),(serverName,serverPort)) #Send to server
			last_seq = count # Set sequence number of last sent packet
			sentence,serverAddr = uc_socket.recvfrom(1024)
		except: # If no ACK then retry and retransfer the packet
			r += 1
			continue
		### try block in order to handle the errors that I encountered
		try: 
			sentence = sentence.decode()
		except:
			r += 1
			continue
		acked = sentence[sentence.find("<ACK>")+5:sentence.find("</ACK>")]
		try:
			int(acked)
		except:
			r += 1
			continue
    	###
		if int(acked) == last_seq+1: #If correct ACK is returned from server than move on to the next packet
			break
		else: # Else, resend again
			r += 1
			continue

uc_socket.close() #Close the socket
print("UDP Transmission Re-transferred Packets: %d" % r)

########################################################
