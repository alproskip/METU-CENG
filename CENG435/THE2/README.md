# Alperen Oğuz Çakmak

# 2237162
 
1-)
Here below are my commands to test this homework in my local machine. (It is written in python 3 but since its my local machine I did not call them with python3)

$ ./simulator 192.168.0.12 192.168.0.12 15584 15585 15586 15587 15588 15589 15590  10 5 3
$ python server.py 15587 15588
$ python client.py 192.168.0.12 15584 15585 15586 15591

After last command, transmission will begin and after the transmission is completed for both communications, terminal with the server.py opened will
show transmission times for both of the communications and terminal with the client.py opened will show number of retransferred packets.

2-)

-In this study I did socket programming in both TCP and UDP protocols.

-I started with the TCP protocol since I think that will be the easy one but apperantly I spent more time in TCP while trying to learn the basics
of socket programming.

-I looked at the lecture slides and the textbooks chapter 2.7 in order to understand the basics.

-I faced all kinds of problems in both python and socket programming basics. Most difficult part for me in the TCP is that realizing that TCP 
communication is continiously streaming the data to the server so that I must arrange the sizes of messages in a way such that server can read them
properly. And for the UDP part, most challanging part is using timeouts in order to implement RDT 3.0

-I learned how to do socket programing and also special cases for TCP and UDP implementaions (Such as importance of packet size, closing the sockets after it is done with the stream, etc.)

-It took me 3 days for all of the homework (including coding and further studying)

3-)
In my RDT my client waited for a ACK message from the server in order to initiate with the next transmission. I used timeout for this purpose.
If my data is somehow corrupted on the way to the server or to the client than timeout raised exception since I did not get any ACK messages from server
and tried to send the same packet to the server again. In my chunks there is also sequence number. All ACK messages have a number that indicates 
the sequence number of the packet that the server is waiting for. So If server is waiting for the previous data that clien just send, then 
client does not send next chunk of data but send the previous one again. At the server side I waited for the correct number of sequence from client
and otherwise did not send any ACK message. If Server is satisfied with both data(no corruption) and sequence number(right order) than server sends
ACK message included with the sequence number that the server is waiting for.
