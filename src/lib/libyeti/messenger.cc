#include "runtime.h"
#include "thread.h"
#include "index.h"
#include "malloc.h"
#include "tensor.h"
#include "node.h"
#include "thread.h"
#include "tensorblock.h"
#include "tensorbranch.h"
#include "taskqueue.h"
#include "aobasis.h"
#include "messenger.h"
#include "exception.h"
#include "env.h"

#include <libsmartptr/strop.h>
#include "data.h"
#include "sort.h"

#if HAVE_MPI
#include <mpi.h>
#endif

#include <time.h>
#include <sys/time.h>


#ifdef redefine_size_t
#define size_t custom_size_t
#endif

using namespace yeti;
using namespace std;

#define TRACK_MESSAGE_HISTORY  0

#define TRACK_RECEIVED_MESSAGES 0

#define DEBUG_MESSENGER 0

#define FINITE_BARRIER_WAIT 0

#undef MICROSECOND_LOCK_WAIT
#define MICROSECOND_LOCK_WAIT 5e6

#define DEBUG_DATA_HEADER 0

#define NELEMENTS_BARRIER 1

#define NELEMENTS_TASK_NOTIFICATION 1

#define DATA_TAG_START 5000
#define DATA_TAG_STOP  32000

static size_t barrier_data[] = {0};
static size_t notification_data[] = {0};

/** Program termination variables */
static sigset_t pthread_sigset;

#define NHEADERS_DATA_SEND 1000000
static SendDataHeader data_send_headers[NHEADERS_DATA_SEND];
static uli data_send_header_num = 0;

#define NHEADERS_DATA_REQUEST 10000
static RequestDataHeader data_req_headers[NHEADERS_DATA_REQUEST];
static uli data_req_header_num = 0;

#define NLISTS_TASK_QUEUE 50
static TaskQueueList task_queue_lists[NLISTS_TASK_QUEUE];
static uli task_queue_list_num = 0;

#if TRACK_RECEIVED_MESSAGES
bool received_tags[50][10000];
bool received_headers[50][10000];
bool sent_headers[50][10000];
uli sent_nmessages[50][10000];
uli received_nmessages[50][10000];
uli actions[10000];
#endif

#if HAVE_MPI
#define MAX_NSTATUS_MESSENGER 10000
static YetiMPIStatus mpi_status_list[MAX_NSTATUS_MESSENGER];
static uli current_status_number = 0;
static uli last_complete = 0;

void
check_status_range(uli start, uli end)
{
    stringstream sstr;
    YetiRuntime::get_messenger()->lock();
    for (uli i=start; i < end; ++i)
    {
        YetiMPIStatus& status = mpi_status_list[i];
        if (!status.is_complete_no_lock())
        {
            sstr << stream_printf("Status %ld is NOT complete on node %ld for sending tag %ld to node %ld\n", 
                    i, YetiRuntime::me(), status.tag, status.proc_number);
        }
        else
        {
            sstr << stream_printf("Status %ld is complete on node %ld for sending tag %ld to node %ld\n", 
                    i, YetiRuntime::me(), status.tag, status.proc_number);
        }
    }
    cerr << sstr.str() << endl;
    YetiRuntime::get_messenger()->unlock();
}

void
check_statuses()
{
    if (last_complete > current_status_number) //looped back around
    {
        for (uli i=last_complete; i < MAX_NSTATUS_MESSENGER; ++i)
        {
            YetiMPIStatus& status = mpi_status_list[i];
            if (!status.is_complete_no_lock())
            {
                last_complete = i;
                return;
            }
        }
        for (uli i=0; i < current_status_number; ++i)
        {
            YetiMPIStatus& status = mpi_status_list[i];
            if (!status.is_complete_no_lock())
            {
                last_complete = i;
                return;
            }
        }
    }
    else //not looped back around
    {
        for (uli i=last_complete; i < current_status_number; ++i)
        {
            YetiMPIStatus& status = mpi_status_list[i];
            if (!status.is_complete_no_lock())
            {
                last_complete = i;
                return;
            }
        }
    }
}
#endif

#if TRACK_MESSAGE_HISTORY
#define NMESSAGES_HISTORY 200
struct MessageHistory {
    uli tag;
    uli dst;
    uli src;
    uli thread;
    uli status_num;
};

uli current_history_num = 0;
MessageHistory histories[NMESSAGES_HISTORY];


void
print_histories()
{
    uli start = current_history_num % NMESSAGES_HISTORY;
    stringstream sstr;
    sstr << stream_printf("Communication history for node %d\n", YetiRuntime::me());
    if (current_history_num > NMESSAGES_HISTORY)
    {
        for (uli i=start; i < NMESSAGES_HISTORY; ++i)
        {
            MessageHistory& his = histories[i];
            sstr << stream_printf("Message %d: tag=%d src=%d dst=%d status=%d thread=%d\n",
                                  i, his.tag, his.src, his.dst, his.status_num, his.thread);

        }
    }
    for (uli i=0; i < start; ++i)
    {
        MessageHistory& his = histories[i];
        sstr << stream_printf("Message %d: tag=%d src=%d dst=%d status=%d thread=%d\n",
                              i, his.tag, his.src, his.dst, his.status_num, his.thread);

    }
    cout << sstr.str() << endl;
}

void
abort_print_histories()
{
    print_histories();
    abort();
}
#endif

YetiMessenger::YetiMessenger()
    :
    barrier_mode_(false),
    keep_running_(false),
    dbl_element_(0),
    first_child_(0),
    second_child_(0),
    parent_(0),
    nchildren_(0),
    current_data_tag_(DATA_TAG_START),
    queue_lock_(0),
    tag_allocate_lock_(0),
    current_tags_(0)
{
    queue_lock_ = new pThreadLock;
    tag_allocate_lock_ = new pThreadLock;

    sigemptyset(&pthread_sigset);
    sigaddset(&pthread_sigset, SIGQUIT);
    sigaddset(&pthread_sigset, SIGTERM);
    sigaddset(&pthread_sigset, SIGINT);
    sigaddset(&pthread_sigset, SIGKILL);
    sigaddset(&pthread_sigset, MAIN_PROCESS_BACKTRACE_SIGNAL);

#if TRACK_RECEIVED_MESSAGES
    ::memset(received_tags, 0, sizeof(received_tags));
    ::memset(received_headers, 0, sizeof(received_headers));
    ::memset(received_nmessages, 0, sizeof(received_nmessages));
    ::memset(sent_headers, 0, sizeof(sent_headers));
    ::memset(sent_nmessages, 0, sizeof(sent_nmessages));
    ::memset(accumulate_tasks_begun, 0, sizeof(accumulate_tasks_begun));
    ::memset(accumulate_tasks_finished, 0, sizeof(accumulate_tasks_finished));
    ::memset(accumulate_tasks_queued, 0, sizeof(accumulate_tasks_queued));
#endif

    current_tags_ = new uli[YetiRuntime::nproc()];
    for (uli i=0; i < YetiRuntime::nproc(); ++i)
        current_tags_[i] = DATA_TAG_START;

    uli me = YetiRuntime::me();
    uli nproc = YetiRuntime::nproc();

    /** Top level node in barrier tree */
    if (me == 0)
    {
        if (nproc > 2)
            nchildren_ = 2;
        else
            nchildren_ = 1;

        first_child_ = 1;
        second_child_ = 2;
    }
    else
    {
        parent_ = (me - 1) / 2;
        first_child_ = 2 * me + 1;
        second_child_ = first_child_ + 1;
        if (first_child_ < nproc) ++nchildren_;
        if (second_child_ < nproc) ++nchildren_;
    }
}

YetiMessenger::~YetiMessenger()
{
    delete queue_lock_;
    delete tag_allocate_lock_;
}

uli
YetiMessenger::allocate_initial_tag(uli nmessages, uli proc_number)
{
    tag_allocate_lock_->lock();
    uli current_tag = current_tags_[proc_number];
    if (current_tag + nmessages >= DATA_TAG_STOP)
        current_tag = DATA_TAG_START;

    uli tag = current_tag;
    current_tag += nmessages;
    current_tags_[proc_number] = current_tag;
    tag_allocate_lock_->unlock();
    return tag;
}

void
YetiMessenger::send_task_notification(uli proc_number)
{
    //cout << stream_printf("Sending task notification to %ld from %ld\n", proc_number, YetiRuntime::me());
    //cout.flush();
    SendStatus* status = send(
        proc_number,
        TaskCompleteNotification,
        notification_data,
        NELEMENTS_TASK_NOTIFICATION * sizeof(size_t)
    );
}

void
YetiMessenger::send_task_acknowledgment(uli proc_number)
{
    //cout << stream_printf("Sending task ack to %ld from %ld\n", proc_number, YetiRuntime::me());
    //cout.flush();
    SendStatus* status = send(
        proc_number,
        TaskCompleteAck,
        notification_data,
        NELEMENTS_TASK_NOTIFICATION * sizeof(size_t)
    );
}

void
YetiMessenger::send_barrier(uli proc_number)
{
    SendStatus* status = send(
        proc_number,
        MessageBarrier,
        barrier_data,
        NELEMENTS_BARRIER * sizeof(size_t)
    );
}

SendStatus*
YetiMessenger::send_data_header(
    uli proc_number,
    uli tag,
    Message::message_data_type_t data_type,
    Message::message_action_t action,
    uli obj_number,
    size_t data_size,
    uli initial_tag,
    uli nmessages
)
{
#if YETI_SANITY_CHECK
    Sendable* obj = get_object(obj_number, data_type);
#endif

    if (data_send_header_num == NHEADERS_DATA_SEND)
        data_send_header_num = 0;

#if TRACK_RECEIVED_MESSAGES
    sent_headers[proc_number][initial_tag] = true;
    sent_nmessages[proc_number][initial_tag] = nmessages;
    for (uli i=initial_tag; i < initial_tag + nmessages; ++i)
        actions[i] = action;
#endif

#if DEBUG_MESSENGER || DEBUG_DATA_HEADER
    cout << stream_printf("Send data header %d for action %d initial tag %d with %d messages object %d to node %d from node %d on thread %d\n",
		data_send_header_num, action, initial_tag,
        nmessages, obj_number, 
        proc_number, YetiRuntime::me(), 
        YetiRuntime::get_thread_number());
    cout.flush();
#endif

    size_t* data = data_send_headers[data_send_header_num].data;

    data[0] = data_type;
    data[1] = action;
    data[2] = obj_number;
    data[3] = data_size;
    data[4] = initial_tag;
    data[5] = nmessages;
    data[6] = data_send_header_num;


    ++data_send_header_num;

    if (!is_locked())
    {
        cerr << "messenger must be locked before sending data header" << endl;
        abort();
    }

    SendStatus* status = send_no_lock(
        proc_number,
        tag,
        data,
        NELEMENTS_SEND_DATA_HEADER * sizeof(size_t)
    );
    //status->set_waiting(false); 
    return status;
}

Sendable*
YetiMessenger::get_object(
    uli object_number,
    Message::message_data_type_t type
)
{
    switch (type)
    {
        case Message::TensorBranch:
            return Malloc<TensorBlock>::get_object(object_number);
    }

    return 0;
}

uli
YetiMessenger::nchildren() const
{
    return nchildren_;
}

uli
YetiMessenger::first_child() const
{
    return first_child_;
}

uli
YetiMessenger::second_child() const
{
    return second_child_;
}

uli
YetiMessenger::parent() const
{
    return parent_;
}

void
YetiMessenger::recv_task_acknowledgment(uli proc_number)
{
    size_t data[NELEMENTS_TASK_NOTIFICATION];
    size_t size = NELEMENTS_TASK_NOTIFICATION * sizeof(size_t);
    data[0] = 0;
    recv(proc_number, TaskCompleteAck, data, size);
    //cout << stream_printf("Received task ack from %ld on %ld\n", proc_number, YetiRuntime::me());
    //cout.flush();
    GlobalQueue::set_node_ack(proc_number);

}

void
YetiMessenger::recv_task_notification(uli proc_number)
{
    size_t data[NELEMENTS_TASK_NOTIFICATION];
    size_t size = NELEMENTS_TASK_NOTIFICATION * sizeof(size_t);
    data[0] = 0;
    recv(proc_number, TaskCompleteNotification, data, size);

    //cout << stream_printf("Received task notification from %ld on %ld\n", proc_number, YetiRuntime::me());
    //cout.flush();


    if (YetiRuntime::use_dynamic_load_balancing())
    {
        cout << stream_printf("Process load balance from %ld on %ld\n",
                               proc_number, YetiRuntime::me());
        cout.flush();
        YetiRuntime::start_timer("dynamic load balance");
        uli ntasks = 0;
        uli nentries_data = 0;
        if (task_queue_list_num == NLISTS_TASK_QUEUE)
            task_queue_list_num = 0;
        uli* task_data = task_queue_lists[task_queue_list_num].data;
        ++task_queue_list_num;
        Task* task = GlobalQueue::get_dynamic_load_balance_segment(ntasks, nentries_data, task_data);
        if (task)
        {
            cout << stream_printf("Dynamic load balance sends %ld tasks to %ld from %ld\n",
                                    ntasks, proc_number, YetiRuntime::me());
            cout.flush();
            uli no_block_number = 0;
            uli initial_tag = DynamicLoadBalanceQueue;
            uli nmessages = 1;
            size_t data_size = nentries_data * sizeof(uli);
            lock();
            SendStatus* status = send_data_header(
                proc_number, 
                GenericDataHeader,
                Message::TaskQueue,
                Message::RunTasks,
                no_block_number,
                data_size,
                initial_tag,
                nmessages
            );
            status = send_no_lock(proc_number, DynamicLoadBalanceQueue, task_data, data_size);
            unlock();
        }
        else //mark the remote node as having finished its work
        {
            send_task_acknowledgment(proc_number);
        }
        YetiRuntime::stop_timer("dynamic load balance");
    }
    else //mark the remote node as having finished its work
    {
        send_task_acknowledgment(proc_number);
    }
    



}

void
YetiMessenger::recv_data_barrier(uli proc_number)
{
    barrier_mode_ = true;

    size_t data[NELEMENTS_BARRIER];
    size_t size = NELEMENTS_BARRIER * sizeof(size_t);
    data[0] = 0;
    recv(proc_number, MessageBarrier, data, size);
    recv_barrier(proc_number);
}

void
YetiMessenger::recv_stop(uli proc_number)
{
    if (proc_number != YetiRuntime::me())
    {
        yeti_throw(SanityCheckError, "Stop received not from self");
    }
    double d;
    tmpl_recv<double>(proc_number, StopListen, d);
    keep_running_ = false;
}

template <typename data_t>
void
YetiMessenger::tmpl_recv(
    uli proc_number,
    uli tag,
    data_t& d
)
{
    recv(proc_number, tag, &d, sizeof(data_t));
}

template <typename data_t>
void
YetiMessenger::tmpl_recv(
    uli proc_number,
    uli tag,
    data_t* d,
    uli n
)
{
    recv(proc_number, tag, d, n * sizeof(data_t));
}

template <typename data_t>
void
YetiMessenger::tmpl_send(
    uli proc_number,
    uli tag,
    data_t d
)
{
    data_t data[1];
    data[0] = d;


    SendStatus* status = tmpl_send<data_t>(proc_number, tag, data, 1);
    status->wait();
}

template <typename data_t>
SendStatus*
YetiMessenger::tmpl_send(
    uli proc_number,
    uli tag,
    data_t* d,
    uli n
)
{
    SendStatus* status = send(proc_number, tag, (void*) d, n * sizeof(data_t));
    return status;
}

SendStatus*
YetiMessenger::send(
    uli proc_number,
    uli tag,
    void *buf,
    size_t size
)
{
    lock();
    SendStatus* status = send_no_lock(proc_number, tag, buf, size);
    unlock();
    return status;
}

void
YetiMessenger::recv(
    uli proc_number,
    uli tag,
    void *buf,
    size_t size
)
{
    lock();
    recv_no_lock(proc_number, tag, buf, size);
    unlock();
}

bool
YetiMessenger::listen_for_tag(uli &proc_number, uli tag)
{
    lock();
    bool flag = listen_for_tag_no_lock(proc_number, tag);
    unlock();
    return flag;
}

bool
YetiMessenger::listen(uli &proc_number, uli &tag)
{
    lock();
    bool flag = listen_no_lock(proc_number, tag);
    unlock();
    return flag;
}

void
YetiMessenger::send(
    uli proc_number,
    uli tag,
    double d
)
{
    tmpl_send<double>(proc_number, tag, d);
}

void
YetiMessenger::send(
    uli proc_number,
    uli tag,
    int i
)
{
    tmpl_send<int>(proc_number, tag, i);
}

SendStatus*
YetiMessenger::send(
    uli proc_number,
    uli tag,
    int* iarr,
    uli n
)
{
    abort();
    SendStatus* status = tmpl_send<int>(proc_number, tag, iarr, n);
    return status;
}

void
YetiMessenger::recv_data_header(
    uli proc_number,
    uli tag
)
{
    YetiRuntime::start_timer("recv data header");
    size_t data[NELEMENTS_SEND_DATA_HEADER];
    size_t size = NELEMENTS_SEND_DATA_HEADER * sizeof(size_t);
    recv(proc_number, tag, data, size);

    Message::message_data_type_t type = (Message::message_data_type_t) data[0];
    Message::message_action_t action = (Message::message_action_t) data[1];
    uli obj_number = (uli) data[2];
    size_t data_size = data[3];
    uli initial_tag = (uli) data[4];
    uli nmessages = (uli) data[5];
    size_t header_num = data[6];

    Sendable* obj = get_object(obj_number, type);

#if TRACK_RECEIVED_MESSAGES
    received_headers[proc_number][initial_tag] = true;
    received_nmessages[proc_number][initial_tag] = nmessages;
#endif

#if DEBUG_MESSENGER || DEBUG_DATA_HEADER
    cout << stream_printf("Recv data header %d on node %d action %d with initial tag %d for %d messages for object %d on thread %d from %d\n",
        header_num,
		YetiRuntime::me(), action, initial_tag, 
        nmessages, obj_number, YetiRuntime::get_thread_number(),
        proc_number
        );
    cout.flush();
#endif
    YetiRuntime::stop_timer("recv data header");

    TensorBlock* block = static_cast<TensorBlock*>(obj);
#if TRACK_RECEIVED_MESSAGES
    stringstream sstr;
#endif
    switch(action)
    {
        case Message::Print:
            block->print(cout);
            cout << endl;
            break;
#if TRACK_RECEIVED_MESSAGES
        case Message::CheckRecv:
            cout << stream_printf(
                "CHECK RECV:Node %ld received tag %ld from node %ld? %d\n", 
                YetiRuntime::me(), initial_tag, proc_number, received_tags[proc_number][initial_tag]
            );
            for (uli i=DATA_TAG_START; i < DATA_TAG_STOP; ++i)
            {
                if (received_tags[proc_number][i])
                {
                    sstr << stream_printf(
                            "Node %ld received tag %ld from node %ld\n",
                            YetiRuntime::me(), i, proc_number, received_tags[proc_number][i]
                    );
                }
                if (received_headers[proc_number][i])
                {
                    sstr << stream_printf(
                            "Node %ld received header for tag %ld from node %ld for %ld messages\n",
                            YetiRuntime::me(), i, proc_number, received_nmessages[proc_number][i]
                    );
                }
                if (accumulate_tasks_queued[proc_number][i])
                {
                    sstr << stream_printf(
                            "Node %ld queued accumulate for tag %ld from node %ld\n",
                            YetiRuntime::me(), i, proc_number
                    );
                }
                if (accumulate_tasks_begun[proc_number][i])
                {
                    sstr << stream_printf(
                            "Node %ld began accumulate for tag %ld from node %ld\n",
                            YetiRuntime::me(), i, proc_number
                    );
                }
                if (accumulate_tasks_finished[proc_number][i])
                {
                    sstr << stream_printf(
                            "Node %ld finished accumulate for tag %ld from node %ld\n",
                            YetiRuntime::me(), i, proc_number
                    );
                }
            }
            GlobalQueue::get_task_queue()->print(sstr);
            sstr << stream_printf("In threaded compute? %d\n", YetiRuntime::is_threaded_compute());
            cerr << sstr.str() << endl;
            cerr << queue_sstr.str() << endl;
            cerr << run_sstr.str() << endl;
            break;
#endif
        case Message::RunTasks:
            YetiRuntime::start_timer("recv load balance");
            {
                LoadBalanceQueue* queue = GlobalQueue::get_dynamic_load_balance_queue();
                if (initial_tag != DynamicLoadBalanceQueue)
                    yeti_throw(SanityCheckError, "Misformatted task queue received");
                recv(proc_number, initial_tag, queue->queue.data, data_size);
                uli nentries = data_size / sizeof(uli);
                uli ntasks = nentries / 5;
                cout << stream_printf("Received %ld tasks on node %ld from %ld\n", ntasks, YetiRuntime::me(), proc_number);
                cout.flush();
                queue->nentries = nentries;
                GlobalQueue::add_dynamic_load_balance_queue(queue);
            }
            YetiRuntime::stop_timer("recv load balance");
            break;
        case Message::AccumulateUnexpected:
            obj->accumulate_data(this, type, proc_number, data_size, initial_tag, nmessages);
            break;
        case Message::Redistribute:
            obj->redistribute(this, type, proc_number, data_size, initial_tag, nmessages);
            break;
        case Message::ReadUnexpected:
            try{
                obj->recv_data(this, type, action, proc_number, data_size, initial_tag, nmessages);
                break;
            } catch (int e) {
                cerr << stream_printf("error in message read unexpected on node %d\n", YetiRuntime::me());
                cerr.flush();
                abort();
            }
        case Message::ReadRequestedRetrieve:
        case Message::ReadRequestedPrefetch:
            YetiRuntime::start_timer("read requested");
            try{
                obj->recv_data(this, type, action, proc_number, data_size, initial_tag, nmessages);
                break;
            } catch (int e) {
                cerr << stream_printf("error in message read requested on node %d\n", YetiRuntime::me());
                cerr.flush();
                abort();
            }
            YetiRuntime::stop_timer("read requested");
        case Message::GlobalSum:
            if (proc_number < YetiRuntime::me()) //this is from the parent
            {
                try{
                    obj->recv_data(this, type, action, proc_number, data_size, initial_tag, nmessages);
                    obj->increment_parent_signal();
                    obj->reset_send_status(); //clear out the child sends to avoid sanity check
                    if (first_child_ < YetiRuntime::nproc())
                    {
                        obj->send_data(this, type, Message::GlobalSum, first_child_);
                    }

                    if (second_child_ < YetiRuntime::nproc())
                    {
                        obj->reset_send_status(); //clear out the child sends to avoid sanity check
                        obj->send_data(this, type, Message::GlobalSum, second_child_);
                    }
                } catch (int e) {
                    cerr << stream_printf("error in global sum child receive on node %d\n",
                                          YetiRuntime::me());
                    cerr.flush();
                    abort();
                }

            }
            else
            {
                try{
                    obj->accumulate_data(this, type, proc_number, data_size, initial_tag, nmessages);
                    obj->set_synced();
                    uli count = obj->increment_child_signal();
                    if (count == nchildren_)
                    {
                        if (YetiRuntime::me() == 0) //i am the top node
                        {
                            obj->increment_parent_signal();
                            if (first_child_ < YetiRuntime::nproc())
                            {
                                obj->send_data(this, type, Message::GlobalSum, first_child_);
                            }

                            if (second_child_ < YetiRuntime::nproc())
                            {
                                obj->reset_send_status();
                                obj->send_data(this, type, Message::GlobalSum, second_child_);
                            }
                        }
                        else
                        {
                            obj->reset_parent_signal();
                            obj->send_data(this, type, Message::GlobalSum, parent_);
                        }
                    }
                } catch (int e) {
                    cerr << stream_printf("error in global sum parent receive on node %d\n",
                                          YetiRuntime::me());
                    cerr.flush();
                    abort();
                }

            }
            break;
        default:
            Env::errn() << "unrecognized action type " << type << endl;
            abort();
    }
}


void
YetiMessenger::send_exit_signal()
{
    lock();
    SendStatus* status = 0;
    size_t data_size = NELEMENTS_REQ_DATA_HEADER * sizeof(size_t);
    for (uli i=0; i < YetiRuntime::nproc(); ++i)
    {
        if (i == YetiRuntime::me())
            continue;

        size_t* data = data_req_headers[data_req_header_num].data;
        data[0] = Message::ExitSignal;
        data[1] = 0;
        data[2] = Message::Exit;
        status = send_no_lock(i, DataRequest, data, data_size);
    }
    unlock();
    status->wait();
    sleep(1);
}

void
YetiMessenger::send_data_request(
    Message::message_data_type_t type,
    uli proc_number,
    uli object_number,
    Message::message_action_t send_action
)
{
    lock();
    if (data_req_header_num == NHEADERS_DATA_REQUEST)
        data_req_header_num = 0;

#if DEBUG_MESSENGER
    cout << stream_printf("Send data request %d for object %d to node %d from node %d on thread %d\n",
                data_req_header_num, object_number, 
                proc_number, YetiRuntime::me(), 
                YetiRuntime::get_thread_number());
    cout.flush();
#endif

    size_t* data = data_req_headers[data_req_header_num].data;
    ++data_req_header_num;

    size_t data_size = NELEMENTS_REQ_DATA_HEADER * sizeof(size_t);
    data[0] = type;
    data[1] = object_number;
    data[2] = send_action;
#if YETI_SANITY_CHECK
    TensorBlock* block = Malloc<TensorBlock>::get_object(object_number);
#if DEBUG_DATA_HEADER
    cout << stream_printf("Send data request for block number %d into branch %p on node %d on thread %d\n", 
                         object_number, block->get_branch(), YetiRuntime::me(), YetiRuntime::get_thread_number());
#endif
    if (!block->is_locked())
    {
        cerr << stream_printf("Requested block %s number %ld is not locked on node %d\n",
                                block->get_block_name().c_str(), object_number, YetiRuntime::me());
        block->controller_fail();
        cerr.flush();
        abort();
    }
    if (block->is_flushable())
    {
        cerr << stream_printf("Requested block %s number %ld is flushable on node %d\n",
                                block->get_block_name().c_str(), object_number, YetiRuntime::me());
        block->controller_fail();
        cerr.flush();
        abort();
    }
#endif
    SendStatus* status = send_no_lock(proc_number, DataRequest, data, data_size);
    //status->set_waiting(false); //do not wait on this
    unlock();
    //ignore return status
}

void
YetiMessenger::recv_data_request(
    uli proc_number
)
{
    size_t data[NELEMENTS_REQ_DATA_HEADER];
    size_t size = NELEMENTS_REQ_DATA_HEADER * sizeof(size_t);
    recv(proc_number, DataRequest, data, size);

    Message::message_data_type_t type = (Message::message_data_type_t) data[0];

    if (type == Message::ExitSignal) //CRASH!!!
    {
        cout << stream_printf("YETI Crash on Node %ld\n", YetiRuntime::me());
        YetiRuntime::get_thread_grp()->thread_crash();
        kill(getpid(), MAIN_PROCESS_BACKTRACE_SIGNAL); //get the main process to print a backtrace
        while(1); //idle until the bad process kills everything
    }


    uli obj_number = (uli) data[1];
    Message::message_action_t send_action = (Message::message_action_t) data[2];


#if DEBUG_MESSENGER
    cout << stream_printf("Recv data request for object %d on node %d from node %d on thread %d\n",
                obj_number, YetiRuntime::me(), proc_number, YetiRuntime::get_thread_number());
    cout.flush();
#endif

    Sendable* obj = get_object(obj_number, type);
    obj->send_data(this, type, send_action, proc_number);
}


void*
run_messenger_recv_thread(void* args)
{
    int s = pthread_sigmask(SIG_BLOCK, &pthread_sigset, NULL);
    if (s != 0)
    {
        cerr << "Unable to set signal mask on recv thread" << endl;
    }

#if USE_DEFAULT_THREAD_STACK
    int x; void* ptr = &x;
    uli recv_threadnum = YetiRuntime::nthread() + 1;
    YetiRuntime::set_thread_stack(recv_threadnum, ptr);
    ++recv_thread_count;
#endif

    YetiMessenger* messenger = static_cast<YetiMessenger*>(args);
    messenger->run_recv_thread();
    return 0;
}

void
YetiMessenger::thread_crash()
{
    lock();
    if (keep_running_) //thread is running
    {
        int status = pthread_kill(recv_thread_, THREAD_KILL_SIGNAL);
        if (status != 0)
        {
            cerr << "Signal " << THREAD_KILL_SIGNAL << " is not a valid pThread kill signal for recv thread" << endl;
        }
    }
    unlock();
}

void
YetiMessenger::start()
{
    lock();
    int status = pthread_attr_init(&recv_thread_attr_);
    if (status != 0)
    {
        cerr << "Unable to init pthread attributes" << endl;
        abort();
    }
#if USE_THREAD_AFFINITY
    uli offset = YetiRuntime::nthread() + YetiRuntime::cpu_mask_num();
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(offset, &cpuset);
    int aff_status = pthread_attr_setaffinity_np(&recv_thread_attr_, sizeof(cpu_set_t), &cpuset);
    if (aff_status != 0)
    {
        cerr << "Unable to init pthread affinity" << endl;
        abort();
    }
    cout << stream_printf("Setting thread affinity recv thread on node %ld to %d\n", YetiRuntime::me(), offset);
    cout.flush();
#endif

    keep_running_ = true;


#if !USE_DEFAULT_THREAD_STACK
    size_t stack_size = YetiRuntime::get_thread_stack_size();
    uli threadnum = YetiRuntime::nthread();
    void* stack_start = YetiRuntime::get_thread_stack(threadnum);
    int stack_status = pthread_attr_setstack(&recv_thread_attr_, stack_start, stack_size);
    if (stack_status != 0)
    {
        cerr << "Unable to set stack address" << endl;
        abort();
    }
#endif

    pthread_create(&recv_thread_, &recv_thread_attr_, run_messenger_recv_thread, this);
    unlock();
}

void
YetiMessenger::stop()
{
    if (!keep_running_) //this has already been stopped
        return;

    send(YetiRuntime::me(), StopListen, 0);
    void* retval;
    pthread_join(recv_thread_, &retval);
}

#define MAX_INCOMING_MESSAGES 10
struct IncomingMessage {
    uli proc_number;
    uli tag;
};

IncomingMessage message_queue[MAX_INCOMING_MESSAGES];
uli nmessages_queued = 0;

void 
YetiMessenger::add_incoming_msg(uli proc_number, uli tag)
{
    IncomingMessage& newmsg = message_queue[nmessages_queued];
    newmsg.proc_number = proc_number;
    newmsg.tag = tag;
    ++nmessages_queued;
}

void
YetiMessenger::check_for_incoming_msg(uli tag)
{
    uli proc_number = 0;
    bool found_message = listen_for_tag_no_lock(proc_number, tag);
    if (found_message)
        add_incoming_msg(proc_number, tag);
}

void
YetiMessenger::build_incoming_msg_queue()
{
    nmessages_queued = 0;
    check_for_incoming_msg(MessageBarrier);
    check_for_incoming_msg(DataRequest);
    check_for_incoming_msg(GenericDataHeader);
    check_for_incoming_msg(RequestedDataHeader);
    if (nmessages_queued == 0) //check for something with a non-standard tag
    {
        uli proc_number = 0;
        uli tag = 0;
        bool found_message = listen_no_lock(proc_number, tag);
        if (found_message)
            add_incoming_msg(proc_number, tag);
    }
}

char array_buffer[1000000];

void
YetiMessenger::recv_msg(uli proc_number, uli tag)
{
#if DEBUG_MESSENGER
    cout << stream_printf("Found tag %d on node %d from node %d on thread %d\n",
                tag, YetiRuntime::me(), proc_number, YetiRuntime::get_thread_number());
    cout.flush();
#endif
    if (tag == RequestedDataHeader || tag == GenericDataHeader)
    {
        //let the recv data header function unlock the recv
        recv_data_header(proc_number, tag);
    }
    else if (tag == DataRequest)
    {
        recv_data_request(proc_number);
        //let the recv data request function unlock the recv
    }
    else if (tag == TaskCompleteAck)
    {
        recv_task_acknowledgment(proc_number);
    }
    else if (tag == TaskCompleteNotification)
    {
        recv_task_notification(proc_number);
    }
    else if (tag == MessageBarrier)
    {
        recv_data_barrier(proc_number);
    }
    else if (tag == StopListen)
    {
        recv_stop(proc_number);
    }
    else if (tag == ReadInteger)
    {
        int x;
        tmpl_recv<int>(proc_number, tag, x);
    }
    else if (tag >= DATA_TAG_START && tag < DATA_TAG_STOP)
    {
        //skip this... there is a data header somewhere else
    }
    else
    {
        uli range = tag / TagRangeDenominator;
        if (range == GlobalSumRange)
        {
            Message::message_data_type_t type =
                static_cast<Message::message_data_type_t>
                    (tag % TagRangeDenominator);
            recv_global_sum(proc_number, tag, type);
        }
        else
        {
            std::string msg = stream_printf("Invalid tag %d received on node %d from node %d",
                     tag, YetiRuntime::me(), proc_number);
            cerr << msg << endl;
            abort();
         }
    }
}

void
YetiMessenger::run_recv_cycle()
{
    lock();
    build_incoming_msg_queue();
    //while we've got the lock, clear out the send statuses
#if HAVE_MPI
    check_statuses();
#endif
    unlock();
    for (uli i=0; i < nmessages_queued; ++i)
    {
        IncomingMessage& msg = message_queue[i];
        recv_msg(msg.proc_number, msg.tag);
    }
    usleep(1); //give the send threads a chance
    nmessages_queued = 0;
}

void
YetiMessenger::run_recv_thread()
{
    uli proc_number = 0;

    uli tag = RequestedDataHeader;

    #define nbuffer 100000
    int* array_buffer[nbuffer];

    uli threadnum = YetiRuntime::get_thread_number();
    uli correct = YetiRuntime::nthread();
    if (threadnum != correct)
    {
        cerr << stream_printf(
            "Messenger thread number error. Recv thread is %d but should be %d",
            threadnum, correct) << endl;
        abort();
    }

    while (keep_running_)
    {
        run_recv_cycle();
    }
}

Sendable::Sendable()
    : is_waiting_(false),
    nsignals_from_children_(0),
    received_parent_signal_(false),
    send_status_(0)
{
}

Sendable::~Sendable()
{
}

void
Sendable::send_data(
    YetiMessenger* messenger,
    Message::message_data_type_t type,
    Message::message_action_t remote_action,
    uli proc_number
)
{
    SendStatus* status = _send_data(
        messenger,
        type,
        remote_action,
        proc_number
    );
    set_send_status(status);
}

void
Sendable::set_synced()
{
    //do nothing
}

void
Sendable::test_callback()
{  
    cout << "successful callback" << endl;
}

void
Sendable::reset_send_status()
{
    send_status_ = 0;
    is_waiting_ = false;
}

void
Sendable::wait_on_send()
{
    if (send_status_)
        send_status_->wait();

    is_waiting_ = false;
    send_status_ = 0;
}

void
Sendable::set_send_status(SendStatus* status)
{
    send_status_ = status;
    if (status)
        is_waiting_ = true;
}

bool
Sendable::has_send_status() const
{
    return send_status_;
}

bool
Sendable::received_parent_signal() const
{
    return received_parent_signal_;
}

void
Sendable::reset_parent_signal()
{
    received_parent_signal_ = false;
}

void
Sendable::increment_parent_signal()
{
    received_parent_signal_ = true;
}

uli
Sendable::nchild_signals() const
{
    return nsignals_from_children_;
}

void
Sendable::reset_child_signal()
{
    nsignals_from_children_ = 0;
}

uli
Sendable::increment_child_signal()
{
    return ++nsignals_from_children_;
}

bool
Sendable::is_waiting() const
{
    return is_waiting_;
}

#if HAVE_MPI
SendStatus*
YetiMPIMessenger::send_no_lock(
    uli proc_number,
    uli tag,
    void* buf,
    size_t size
)
{
    if (proc_number >= YetiRuntime::nproc())
        abort();

#if TRACK_MESSAGE_HISTORY
    if (tag != MessageBarrier)
    {
        uli hisnum = current_history_num % NMESSAGES_HISTORY;
        MessageHistory& his = histories[hisnum];
        his.src = proc_number;
        his.dst = YetiRuntime::me();
        his.tag = tag;
        his.status_num = current_status_number;
        his.thread = YetiRuntime::get_thread_number();
        ++current_history_num;
    }
#endif

    YetiMPIStatus& status = mpi_status_list[current_status_number];
    bool last_send_done = status.is_complete_no_lock();
    if (!last_send_done)
    {
        cerr << stream_printf("Too many sends accumulated on messenger on node %ld\n",
                             YetiRuntime::me());
        cerr.flush();
        abort();
    }

    status.clear_dependents();
    status.set_waiting(true);

    ++current_status_number;
    if (current_status_number == MAX_NSTATUS_MESSENGER) 
        current_status_number = 0;

#if DEBUG_MESSENGER
    cout << stream_printf("Starting send of tag %d to node %d from %d on thread %d\n",
                tag, proc_number, YetiRuntime::me(), YetiRuntime::get_thread_number());
    cout.flush();
#elif DEBUG_DATA_HEADER
    if (tag >= DATA_TAG_START && tag < DATA_TAG_STOP)
    {
        cout << stream_printf("Starting send of tag %d to node %d from %d on thread %d\n",
                    tag, proc_number, YetiRuntime::me(), YetiRuntime::get_thread_number());
        cout.flush();
    }
    else if (tag == DynamicLoadBalanceQueue)
    {
        cout << stream_printf("Starting send of tag %d to node %d from %d on thread %d\n",
                    tag, proc_number, YetiRuntime::me(), YetiRuntime::get_thread_number());
        cout.flush();
    }
#endif

    int err;
    if (size % sizeof(double) == 0) //send as double
    {
        err = MPI_Isend(
                    buf,
                    size / sizeof(double),
                    MPI_DOUBLE,
                    proc_number,
                    tag,
                    YetiRuntime::get_mpi_comm(),
                    &status.req
                 );
    }
    else if (size % sizeof(int) == 0) //send as int
    {
        err = MPI_Isend(
                    buf,
                    size / sizeof(int),
                    MPI_INT,
                    proc_number,
                    tag,
                    YetiRuntime::get_mpi_comm(),
                    &status.req
                 );
    }
    else
    {
        err = MPI_Isend(
                    buf,
                    size,
                    MPI_BYTE,
                    proc_number,
                    tag,
                    YetiRuntime::get_mpi_comm(),
                    &status.req
                 );
    }

    if (err != MPI_SUCCESS)
    {
        Env::errn() << "MPI error in Send" << err << endl;
        abort();
    }

#if DEBUG_MESSENGER
    cout << stream_printf("Finished send of tag %d to node %d from %d on thread %d\n",
                tag, proc_number, YetiRuntime::me(), YetiRuntime::get_thread_number());
    cout.flush();
#endif


    status.msg_size = size;
    status.proc_number = proc_number;
    status.tag = tag;
    return &status;
}

void
YetiMPIMessenger::recv_no_lock(
    uli proc_number,
    uli tag,
    void *buf,
    size_t size
)
{
#if DEBUG_MESSENGER
    cout << stream_printf("Starting recv of tag %d into buffer %p on node %d from %d on thread %d\n",
                tag, buf, YetiRuntime::me(), proc_number, YetiRuntime::get_thread_number()); 
    cout.flush();
#elif DEBUG_DATA_HEADER
    if (tag >= DATA_TAG_START && tag < DATA_TAG_STOP)
    {
        cout << stream_printf("Starting recv of tag %d into buffer %p on node %d from %d on thread %d\n",
                    tag, buf, YetiRuntime::me(), proc_number, YetiRuntime::get_thread_number()); 
        cout.flush();
    }
#endif

#if TRACK_MESSAGE_HISTORY
    if (tag != MessageBarrier)
    {
        uli hisnum = current_history_num % NMESSAGES_HISTORY;
        MessageHistory& his = histories[hisnum];
        his.src = proc_number;
        his.dst = YetiRuntime::me();
        his.tag = tag;
        his.status_num = 0;
        his.thread = YetiRuntime::get_thread_number();
        ++current_history_num;
    }
#endif

    MPI_Status status;
    int err;
    if (size % sizeof(double) == 0)
    {
        err = MPI_Recv(
                    buf, 
                    size / sizeof(double), 
                    MPI_DOUBLE, 
                    proc_number, 
                    tag, 
                    YetiRuntime::get_mpi_comm(), 
                    &status
                  );
    }
    else if (size % sizeof(int) == 0)
    {
        err = MPI_Recv(
                    buf, 
                    size / sizeof(int), 
                    MPI_INT, 
                    proc_number, 
                    tag, 
                    YetiRuntime::get_mpi_comm(), 
                    &status
                  );
    }
    else
    {
        err = MPI_Recv(
                    buf, 
                    size, 
                    MPI_CHAR, 
                    proc_number, 
                    tag, 
                    YetiRuntime::get_mpi_comm(), 
                    &status
                  );
    }
    if (err != MPI_SUCCESS)
    {
        Env::errn() << "MPI error in Recv " << err << endl;
        abort();
    }
#if DEBUG_MESSENGER
    cout << stream_printf("Finished recv of tag %d on node %d from %d\n",
                tag, YetiRuntime::me(), proc_number);
    cout.flush();
#endif

#if TRACK_RECEIVED_MESSAGES
    received_tags[proc_number][tag] = 1;
#endif
}

bool
YetiMPIMessenger::listen_no_lock(
    uli& proc_number,
    uli& tag
)
{
    MPI_Status status;
    int flag;
    int err = MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, YetiRuntime::get_mpi_comm(), &flag, &status);
    if (err != MPI_SUCCESS)
    {
        Env::errn() << "MPI error in Probe" << err << endl;
        abort();
    }

    if (!flag)
        return false;

    proc_number = status.MPI_SOURCE;
    tag = (uli) status.MPI_TAG;

    return true;
}

bool
YetiMPIMessenger::listen_for_tag_no_lock(
    uli& proc_number,
    uli tag
)
{
    MPI_Status status;
    int flag;
    int err = MPI_Iprobe(MPI_ANY_SOURCE, tag, YetiRuntime::get_mpi_comm(), &flag, &status);
    if (err != MPI_SUCCESS)
    {
        Env::errn() << "MPI error in Probe" << err << endl;
        abort();
    }

    if (!flag)
        return false;

    proc_number = status.MPI_SOURCE;

    return true;
}

void
YetiMPIMessenger::global_sum(double& element)
{
    uli nproc = YetiRuntime::nproc();
    uli me = YetiRuntime::me();

    if (nproc == 1)
        return;

    uli wait_count = 0;
    while (nsignals_from_children_ < nchildren_)
    {
        usleep(1);
        ++wait_count;
#if FINITE_BARRIER_WAIT
        if (wait_count > MICROSECOND_LOCK_WAIT * 10)
        {
            cerr << "Global sum child wait error" << endl;
            throw TENSOR_WAIT_ERROR;
        }
#endif
    }

    element += dbl_element_;
    if (me != 0)
    {
        //prepare to receive double from parent
        YetiMessenger::send(parent_, GlobalSumDouble, element);
        uli wait_count = 0;
        while (nsignals_from_parent_ == 0)
        {
            usleep(1);
#if FINITE_BARRIER_WAIT
            ++wait_count;
            if (wait_count > MICROSECOND_LOCK_WAIT)
            {
                cerr << "Global sum parent wait error" << endl;
                throw TENSOR_WAIT_ERROR;
            }
#endif
        }
        element = dbl_element_;
    }
    else; //parent processor, already done

    //reset for the next global sum
    dbl_element_ = 0;
    nsignals_from_children_ = 0;
    nsignals_from_parent_ = 0;

    if (first_child_ < nproc)
    {
        YetiMessenger::send(first_child_, GlobalSumDouble, element);
    }
    if (second_child_ < nproc)
    {
        YetiMessenger::send(second_child_, GlobalSumDouble, element);
    }
}

void
YetiMPIMessenger::recv_global_sum(
    uli proc_number,
    uli tag,
    Message::message_data_type_t type
)
{
    char data[50];
    switch(type)
    {
        case Message::Double:
            double elem;
            tmpl_recv<double>(proc_number, tag, elem);
            if (proc_number < YetiRuntime::me())
            {
                dbl_element_ = elem;
                ++nsignals_from_parent_;
            }
            else
            {
                dbl_element_ += elem;
                ++nsignals_from_children_;
            }
            break;

    }
}

void
YetiMPIMessenger::wait_barrier()
{
    uli nproc = YetiRuntime::nproc();
    uli me = YetiRuntime::me();

    if (nproc == 1)
        return;

    uli wait_count = 0;
    while (nsignals_from_children_ < nchildren_)
    {
        usleep(1);
        ++wait_count;
#if FINITE_BARRIER_WAIT
        if (wait_count > MICROSECOND_LOCK_WAIT)
        {
            cerr << "Barrier child wait error" << endl;
            throw TENSOR_WAIT_ERROR;
        }
#endif
    }

    if (me != 0)
    {
        send_barrier(parent_);
        uli wait_count = 0;
        while (nsignals_from_parent_ == 0)
        {
            usleep(1);
#if FINITE_BARRIER_WAIT
            ++wait_count;
            if (wait_count > MICROSECOND_LOCK_WAIT)
            {
                cerr << "Barrier parent wait error" << endl;
                throw TENSOR_WAIT_ERROR;
            }
#endif
        }
    }

    //reset for the next barrier
    nsignals_from_children_ = 0;
    nsignals_from_parent_ = 0;

    if (first_child_ < nproc)
    {
        send_barrier(first_child_);
    }
    if (second_child_ < nproc)
    {
        send_barrier(second_child_);
    }
    barrier_mode_ = false;
}

void
YetiMPIMessenger::recv_barrier(uli proc_number)
{
    if (proc_number < YetiRuntime::me())
    {
        //all processes have reached the barrier
        //this signal has emanated from the parent
        ++nsignals_from_parent_;
    }
    else
    {
        //signal has come from child
        //signal needs to be passed onto parent
        //before completing barrier
        ++nsignals_from_children_;
    }
}


bool
YetiMPIStatus::is_complete_no_lock()
{
#if YETI_SANITY_CHECK
    if (!YetiRuntime::get_messenger()->is_locked())
    {
        yeti_throw(SanityCheckError, "Messenger is not locked in is_complete_no_lock");
    }
#endif
    if (!is_waiting())
        return true;

    int flag;
    /** This only ensures that the message is completely sent, not that it is received */
    int err = MPI_Test(&req, &flag, MPI_STATUS_IGNORE);
    if (err != MPI_SUCCESS)
    {
        if (err == MPI_ERR_ARG)
            cerr << stream_printf("MPI_ERR_ARG in complete check on thread %d node %d\n", 
                                    YetiRuntime::get_thread_number(), YetiRuntime::me());
        else if (err == MPI_ERR_REQUEST)
            cerr << stream_printf("MPI_ERR_REQUEST in complete check on thread %d node %d\n", 
                                    YetiRuntime::get_thread_number(), YetiRuntime::me());
    }


    if (flag)
    {
        uli send_number = (size_t) this - (size_t) mpi_status_list;
        send_number /= sizeof(YetiMPIStatus);
        if (last_complete < send_number)
            last_complete = send_number;
        set_waiting(false);
    }
    return flag;
}

bool
YetiMPIStatus::is_complete()
{
    messenger->lock();
    bool flag = is_complete_no_lock();
    messenger->unlock();
    return flag;
}
#endif

SendStatus*
YetiLocalMessenger::send_no_lock(
    uli proc_number,
    uli tag,
    void *buf,
    size_t size
)
{
    yeti_throw(SanityCheckError, "should never call send on local messenger");
    return 0;
}

void
YetiLocalMessenger::recv_no_lock(
    uli proc_number,
    uli tag,
    void *buf,
    size_t size
)
{
    yeti_throw(SanityCheckError, "should never call recv on local messenger");
}

bool
YetiLocalMessenger::listen_no_lock(
    uli& proc_number,
    uli& tag
)
{
    yeti_throw(SanityCheckError, "should never call listen on local messenger");
    return false;
}

bool
YetiLocalMessenger::listen_for_tag_no_lock(
    uli& proc_number,
    uli tag
)
{
    yeti_throw(SanityCheckError, "should never call listen on local messenger");
    return false;
}

void
YetiLocalMessenger::global_sum(double& element)
{
    yeti_throw(SanityCheckError, "should never call global sum on local messenger");
}

void
YetiLocalMessenger::wait_barrier()
{
    //do nothing
}

void
YetiLocalMessenger::recv_barrier(uli proc_number)
{
    yeti_throw(SanityCheckError, "should never call recv barrier on local messenger");
}

void
YetiLocalMessenger::stop()
{
}

void
YetiLocalMessenger::start()
{
}

void
YetiLocalMessenger::recv_global_sum(
    uli proc_number,
    uli tag,
    Message::message_data_type_t type
)
{
}

YetiLocalMessenger::YetiLocalMessenger()
    : YetiMessenger()
{
}

SendStatus::SendStatus()
    : 
    dependent_status_(0),
    waiting_(false),
    messenger(0)
{
}

bool
SendStatus::is_waiting() const
{
    return waiting_;
}

void
SendStatus::set_waiting(bool flag)
{
    waiting_ = flag;
}

void
SendStatus::clear_dependents()
{
#if YETI_SANITY_CHECK
    if (!YetiRuntime::get_messenger()->is_locked())
    {
        yeti_throw(SanityCheckError, "Messenger is not locked in send status reset");
    }
#endif
    dependent_status_ = 0;
}

void
SendStatus::set_dependent(SendStatus* status)
{
    dependent_status_ = status;
}

void
SendStatus::wait()
{
#if HAVE_MPI
    uli wait_count = 0;
    uli initial_last_complete = last_complete;
    uli current_complete = 0;
    while(waiting_ && !is_complete())
    {
        if (current_complete != last_complete)
        {
            current_complete = last_complete;
            wait_count = 0;
        }
        usleep(1);
        ++wait_count;
        uli send_number = (size_t) this - (size_t) mpi_status_list;
        send_number /= sizeof(YetiMPIStatus);
        if (wait_count > MICROSECOND_LOCK_WAIT)
        {
#if TRACK_RECEIVED_MESSAGES
            YetiRuntime::get_messenger()->lock();
            YetiRuntime::get_messenger()->send_data_header(
                proc_number,
                GenericDataHeader,
                Message::TensorBranch,
                Message::CheckRecv,
                0,
                msg_size,
                tag + 20000,
                0
            );
            YetiRuntime::get_messenger()->unlock();

            cerr << stream_printf("Send status wait error on status number %ld. Current status number is %ld.\n"
                                  "Trying and failing to send %ld bytes to node %ld for tag %ld with action %d\n"
                                  "The send status advanced from %ld to %ld over the course of the wait\n",
                                send_number, current_status_number, 
                                msg_size, proc_number, tag, actions[tag],
                                initial_last_complete, last_complete);
#else
            cerr << stream_printf("Send status wait error on status number %ld. Current status number is %ld.\n"
                                  "Trying and failing to send %ld bytes to node %ld for tag %ld\n"
                                  "The send status advanced from %ld to %ld over the course of the wait\n",
                                send_number, current_status_number, 
                                msg_size, proc_number, tag, 
                                initial_last_complete, last_complete);
#endif
            YetiRuntime::stack_print();
            uli start = send_number < initial_last_complete ? send_number : initial_last_complete;
            check_status_range(start, current_status_number);
#if TRACK_RECEIVED_MESSAGES
            stringstream sstr;
            for (uli i=DATA_TAG_START; i < DATA_TAG_STOP; ++i)
            {
                if (sent_headers[proc_number][i])
                {
                    sstr << stream_printf("Sent header for tag %ld from node %ld to node %ld for %ld messages\n",
                                            i, YetiRuntime::me(), proc_number,sent_nmessages[proc_number][i]);
                }
            }
            cerr << sstr.str() << endl;
#endif
            sleep(5);
            throw TENSOR_WAIT_ERROR;
        }
    }

    if (dependent_status_)
        dependent_status_->wait();
#endif //HAVE_MPI_H
}


#if HAVE_MPI
#undef yeti_throw
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>


YetiMPIMessenger::YetiMPIMessenger()
    : YetiMessenger(),
    nsignals_from_children_(0),
    nsignals_from_parent_(0)
{
    for (uli i=0; i < MAX_NSTATUS_MESSENGER; ++i)
    {
        mpi_status_list[i].messenger = this;
    }
}
#endif

