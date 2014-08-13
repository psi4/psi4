/* $Header: /tmp/hpctools/ga/tcgmsg/ipcv4.0/sockets.h,v 1.5 2000-09-30 19:04:22 d3g681 Exp $ */

extern void ShutdownAll();
extern int ReadFromSocket();
extern int WriteToSocket();
extern void CreateSocketAndBind();
extern int ListenAndAccept();
extern int CreateSocketAndConnect();
extern long PollSocket();
extern long WaitForSockets(int nsock, int *socks, int *list);
