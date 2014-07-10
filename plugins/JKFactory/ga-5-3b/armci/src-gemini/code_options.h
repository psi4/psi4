/* 
   Questions:
   ORNL - tipparajuv@ornl.gov 
   CRAY - ryan@cray.com
*/

/* --------------------------------------------------------------------------- *\
   PORTALS_USE_RENDEZ_VOUS
   =======================
   When the number of PEs gets very large, the data server is required to have
   buffer space available for all possible incoming messages which is defined
   by PORTALS_MAX_DESCRIPTORS = (MAX_BUFS+MAX_SMALL_BUFS).
   For each PE, the DS must have at least:
     min_memory_per_pe = PORTALS_MAX_BUFS*PORTALS_BUF_SIZE +
                         PORTALS_MAX_SMALL_BUFS*PORTALS_SMALL_BUF_SIZE
   This becomes a memory constraint at large core count.
   Rendez-vous message is one mechanism to get around requiring the DS to
   have buffer space for all messages. When rendez-vous (RZV) messaging is
   enabled, the messages what use the large buffers no longer send the entire
   buffer "eagerly".  Instead, only the data request (request_header_t) gets
   sent to the data server.  When the data server is ready to handle the
   request, it "pulls" the entire buffer over via a portals_get operation.
   One can immediately see that this can lead to a slow down in performance,
   since the data server is idle when it has to pull the data over.  This is
   the price paid when you remove the bufferign for those messsages.  Ideally,
   when the DS is pulling the message, it could be processing another request.
   This double buffering technique needs to be programmed in. Care must be
   taken to ensure proper ARMCI behavior. The next request handled can not be
   from the same PE, nor can it be a FENCE operation ... all other (?)
   requests/operations can be double buffered.
\* --------------------------------------------------------------------------- */
 # define PORTALS_USE_RENDEZ_VOUS



/* --------------------------------------------------------------------------- *\
   PORTALS_LIMIT_REMOTE_REQUESTS_BY_NODE
   =====================================
   Another means to reduce the required buffer needed by the data server is
   to limit the number of cores that can talk to the data server at any given
   moment.  When this options is turned on, only 1 request per node is allowed
   to be in the buffer of any given data server.  On a 10 core node, the size
   of the buffer required by the data server is reduced by more than an order
   of magnitude.  You get more than an order of magnitude, because you don't
   need to reserve space for any of the small buffers, since you can only have
   one small or one large from any given node in the ds buffer at any one time.
   Another major benefit is you can increase MAX_BUFS and MAX_SMALL_BUFS to
   increase concurrency without affecting the DS's buffer size.

   Can be used with PORTALS_USE_RENDEZ_VOUS.

        notes: every request needs to respond with an ack, even gets.
        acks actually send data when we limit remote request ... the ack
        response is needed to trigger that the outstanding request has
        been finished by the data server ... the ack zeros out the index
        in the active_requests_by_node array.
\* --------------------------------------------------------------------------- */
 # define PORTALS_LIMIT_REMOTE_REQUESTS_BY_NODE_TURNED_OFF


/* --------------------------------------------------------------------------- *\
   PORTALS_AFFINITY
   ================
   When initializing compute processes and data servers, the affinity passed
   in by aprun/alps is ignored.

   Compute processes are bound strictly to a particular core.  Cores are
   evenly divided between sockets keeping the last core (mask = 1 << (ncpus-1))
   free for the data server.

   If the node is not fully subscribed, then the data server is bound to the
   last core on the node (mask = 1 << (ncpus-1)); otherwise, the data server
   is "free floating" (mask = (1 << ncpus)-1) on a fully subscribed node.
\* --------------------------------------------------------------------------- */
 # define PORTALS_AFFINITY
 # define PORTALS_AFFINITY_NSOCKETS 2


/* --------------------------------------------------------------------------- *\
   CRAY_USE_MDMD_COPY 
   ================== 
   Used MDMD copy instead of PtlGetRegion for on-node "local" transfers
\* --------------------------------------------------------------------------- */
 # define CRAY_USE_MDMD_COPY
 


/* --------------------------------------------------------------------------- *\
   ORNL_USE_DS_FOR_REMOTE_GETS 
   =========================== 
   Vinod informed us of a modification that can be made to enable the use of
   the data server for remote gets.  Without this option, direct gets are 
   used.  This can cause severe network congestion, because many armci_gets
   are not stride 1.  The data server packs those gets into contiguous blocks
   and sends them back as a single put.  However, the direct gets, require
   many small messages.

   Unfortunately, there is a small bug in the DS for remote gets.  This bug
   may cause the program to abort or print out the following message:
   %d: server wrote data at unexpected offset %d

   This is a bug actively being worked on @ CRAY and ORNL.
\* --------------------------------------------------------------------------- */
 # define ORNL_USE_DS_FOR_REMOTE_GETS   
 # define CRAY_USE_ARMCI_CLIENT_BUFFERS 
