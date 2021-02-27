from socket import *
import time
import sys
import hashlib


serverPort = int(sys.argv[2]) #TCP SERVER Port

ts_socket = socket(AF_INET, SOCK_STREAM) #Create TCP socket
ts_socket.bind(('',serverPort)) #Bind it to the port

ts_socket.listen(1) #Listen to the client (simulator in our case)

while True:
	cSocket, addr = ts_socket.accept() #Create socket that will accept incoming data
	all_file = ''
	count = 0
	t = ''
	times = 0
	while True:
		b_flag = 0
		count += 1
		sentence = ''
		while len(sentence) < 1000: #Since TCP is streaming the data continuously I recieve data as chunks of 1000bytes from the stream
			sentence += cSocket.recv(5).decode()

		if sentence.find("</ed>") != -1: # Check the tag in order to decide end of the packet
			b_flag=1 #Set the break flag

		timestamp = time.time()
		t = sentence[sentence.find("<time>")+6:sentence.find("<msg>")-1]
		time_passed = timestamp-float(t) # Calculate the time
		
		rcv_msg = sentence[sentence.find("<msg>")+5:sentence.find("</end>")] # Extract the message
		all_file += rcv_msg #Add chunk of message to the whole message
		times += time_passed #Update time variable
	
		if b_flag: # Break Flag closes the file and exits inner while loop
			break

	f = open("Transfered_TCP_File.txt",'w') #Write to the new file
	f.write(all_file)
	f.close()

	times = times/10**3 #Conver to miliseconds
	print("TCP Communication Total Transmission Time: %.8f ms" % times)
	print("TCP Packet Average Transmission Time: %.8f ms" % (times/(count-1)))

	break 

ts_socket.close() #Close the socket

#### UDP SERVER IMPLEMENTATION ####


serverPort = int(sys.argv[1])
us_socket = socket(AF_INET, SOCK_DGRAM)
us_socket.bind(('',serverPort))

while True:
	all_file = ''
	count = 0
	t = ''
	sentence = ''
	times = 0
	waiting_for = 1
	while True:
		NAK = 0 #NAK flag is used to make the client send the chunk again
		b_flag = 0 
		count += 1
		
		message, clientAddr = us_socket.recvfrom(2048) #Receive message from port
		
		try: #Try to decode the message
			sentence = message.decode()
		except: # If not, set NAK flag 
			NAK = 1

		t = sentence[sentence.find("<time>")+6:sentence.find("<cs>")] #Extract the time
		try:
			time_passed = time.time()-float(t) # Calculate the time
		except: #If cant extract it, set NAK flag (data is corrupted)
			NAK = 1

		rcv_msg = sentence[sentence.find("<msg>")+5:] # Extract the message
		
		if sentence.find("<end>")!=-1: #If this is last chunk, set b_flag and erase end tag from message
			rcv_msg = rcv_msg[0:-5]
			b_flag = 1
		
		rcv_seq = sentence[sentence.find("<seq>")+5:sentence.find("<msg>")] #Extract sequence number
		rcv_cs = sentence[sentence.find("<cs>")+4:sentence.find("<seq>")] #Extract checksum
		checksum = hashlib.md5(rcv_msg.encode()).hexdigest() #Calculate checksum for coming message

		if NAK==0 and checksum==rcv_cs: #If checksum is correct and NAK flag is not set then everything is okay
			NAK = 0
		else:
			NAK = 1

		if NAK: #If data is wrong, decrement count because we could not get the correct message
			count -= 1
		else: # IF everything is okay, check for the sequence number
			if int(waiting_for)  == int(rcv_seq): # If received sequence number is equal to what we are waiting for then add chunk to the whole message
				waiting_for += 1
				all_file += rcv_msg
				times += time_passed #Update times variable
			
			respond = "<ACK>"+str(waiting_for)+"</ACK>" #Add waiting_for variable to the respond in order to let the client know that what the server is waiting for
			us_socket.sendto(respond.encode(),clientAddr)
		
		if b_flag: # Break Flag closes the file and exits inner while loop
			break

	with open('Transfered_UDP_File.txt', 'w') as f: #Write to the file
		f.write(all_file)

	times = times*1000
	print("UDP Communication Total Transmission Time: %.8f ms" % times)
	print("UDP Packets Average Transmission Time: %.8f ms" % (times/(count-1)))

	break
us_socket.close() #Close the socket

#############################################
