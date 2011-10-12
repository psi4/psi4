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

#include <psiconfig.h>
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


#define NELEMENTS_DATA_REQUEST 2
#define NELEMENTS_BARRIER 1

#define DATA_TAG_START 5000
#define DATA_TAG_STOP  10000


YetiMessenger::YetiMessenger()
    :
    keep_running_(false),
    send_queue_put_number_(0),
    send_queue_read_number_(0),
    dbl_element_(0),
    first_child_(0),
    second_child_(0),
    parent_(0),
    nchildren_(0),
    current_data_tag_(DATA_TAG_START),
    queue_lock_(0),
    tag_allocate_lock_(0)
{
    for (uli i=0; i < NREQUESTS_MAX_MESSENGER; ++i)
    {
        SendRequest& send_req = send_request_queue_[i];
        send_req.free = true;
    }

    queue_lock_ = new pThreadLock;
    tag_allocate_lock_ = new pThreadLock;

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
YetiMessenger::allocate_initial_tag(uli nmessages)
{
    heisenfxn(YetiMessenger::allocate_initial_tag);
    tag_allocate_lock_->lock();
    if (current_data_tag_ + nmessages >= DATA_TAG_STOP)
        current_data_tag_ = DATA_TAG_START;

    uli tag = current_data_tag_;
    current_data_tag_ += nmessages;
    tag_allocate_lock_->unlock();
    return tag;
}

void
YetiMessenger::send_barrier(uli proc_number)
{
    heisenfxn(YetiMessenger::send_barrier);
    size_t data[NELEMENTS_BARRIER];
    data[0] = 0;
    SendStatus* status = send(
        proc_number,
        MessageBarrier,
        data,
        NELEMENTS_BARRIER * sizeof(size_t)
    );
    status->wait();
}

SendStatus*
YetiMessenger::send_data_header(
    SendDataHeader* header,
    uli proc_number,
    Message::message_data_type_t data_type,
    Message::message_action_t action,
    uli obj_number,
    size_t data_size,
    uli initial_tag,
    uli nmessages
)
{
    heisenfxn(YetiMessenger::send_data_header);
    size_t* data = header->data;

    data[0] = data_type;
    data[1] = action;
    data[2] = obj_number;
    data[3] = data_size;
    data[4] = initial_tag;
    data[5] = nmessages;

    if (!is_locked())
    {
        cerr << "messenger must be locked before sending data header" << endl;
        abort();
    }

    cout << stream_printf("starting send data header for initial tag %d for %d messages to node %d from node %d thread %d\n",
		initial_tag, nmessages, proc_number, YetiRuntime::me(), YetiRuntime::get_thread_number());
    cout.flush();

    SendStatus* status = send_no_lock(
        proc_number,
        DataHeader,
        data,
        NELEMENTS_DATA_HEADER * sizeof(size_t)
    );

    //wait to make damn sure this gets sent
    usleep(100);

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
YetiMessenger::recv_data_barrier(uli proc_number)
{
    heisenfxn(YetiMessenger::recv_data_barrier);
    size_t data[NELEMENTS_BARRIER];
    size_t size = NELEMENTS_BARRIER * sizeof(size_t);
    data[0] = 0;
    recv(proc_number, MessageBarrier, data, size);
    recv_barrier(proc_number);
}

void
YetiMessenger::recv_stop(uli proc_number)
{
    heisenfxn(YetiMessenger::recv_stop);
    if (proc_number != YetiRuntime::me())
    {
        raise(SanityCheckError, "stop received not from self");
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
    heisenfxn(YetiMessenger::tmpl_recv);
    recv(proc_number, tag, &d, sizeof(data_t));
}

template <typename data_t>
void
YetiMessenger::tmpl_send(
    uli proc_number,
    uli tag,
    data_t d
)
{
    heisenfxn(YetiMessenger::tmpl_send);
    data_t data[1];
    data[0] = d;

    SendStatus* status = send(proc_number, tag, data, sizeof(data_t));
    status->wait();
}

SendStatus*
YetiMessenger::send(
    uli proc_number,
    uli tag,
    void *buf,
    size_t size
)
{
    heisenfxn(YetiMessenger::send);
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
    heisenfxn(YetiMessenger::recv);
    lock();
    recv_no_lock(proc_number, tag, buf, size);
    unlock();
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
    heisenfxn(YetiMessenger::send);
    tmpl_send<double>(proc_number, tag, d);
}

void
YetiMessenger::send(
    uli proc_number,
    uli tag,
    int i
)
{
    heisenfxn(YetiMessenger::send);
    tmpl_send<int>(proc_number, tag, i);
}

void
YetiMessenger::recv_data_header(
    uli proc_number
)
{
    heisenfxn(YetiMessenger::recv_data_header);
    size_t data[NELEMENTS_DATA_HEADER];
    size_t size = NELEMENTS_DATA_HEADER * sizeof(size_t);
    cout << stream_printf("starting recv data header from node %d on node %d thread %d\n",
		proc_number, YetiRuntime::me(), YetiRuntime::get_thread_number());
    cout.flush();
    recv(proc_number, DataHeader, data, size);

    Message::message_data_type_t type = (Message::message_data_type_t) data[0];
    Message::message_action_t action = (Message::message_action_t) data[1];
    uli obj_number = (uli) data[2];
    size_t data_size = data[3];
    uli initial_tag = (uli) data[4];
    uli nmessages = (uli) data[5];

    Sendable* obj = get_object(obj_number, type);

    cout << stream_printf("finishing recv data header on node %d action %d with initial tag %d for %d messages\n",
		YetiRuntime::me(), action, initial_tag, nmessages);
    cout.flush();

    switch(action)
    {
        case Message::Accumulate:
            obj->accumulate_data(this, type, proc_number, data_size, initial_tag, nmessages);
            break;
        case Message::Read:
        case Message::Write:
            obj->recv_data(this, type, proc_number, data_size, initial_tag, nmessages);
            break;
        case Message::GlobalSum:
            if (proc_number < YetiRuntime::me()) //this is from the parent
            {
                obj->recv_data(this, type, proc_number, data_size, initial_tag, nmessages);
                obj->set_synced();
                obj->increment_parent_signal();
                if (first_child_ < YetiRuntime::nproc())
                    queue_send_request(type, Message::GlobalSum, first_child_, obj);

                if (second_child_ < YetiRuntime::nproc())
                    queue_send_request(type, Message::GlobalSum, second_child_, obj);

            }
            else
            {
                obj->accumulate_data(this, type, proc_number, data_size, initial_tag, nmessages);
                uli count = obj->increment_child_signal();
                if (count == nchildren_)
                {
                    if (YetiRuntime::me() == 0) //i am the top node
                    {
                        obj->increment_parent_signal();
                        if (first_child_ < YetiRuntime::nproc())
                            queue_send_request(type, Message::GlobalSum, first_child_, obj);

                        if (second_child_ < YetiRuntime::nproc())
                            queue_send_request(type, Message::GlobalSum, second_child_, obj);
                    }
                    else
                    {
                        obj->reset_parent_signal();
                        queue_send_request(type, Message::GlobalSum, parent_, obj);
                    }
                }

            }
            break;
        default:
            Env::errn() << "unrecognized action type " << type << endl;
            abort();
    }

    heisenfxn(YetiMessenger::recv_data_header);
}

void
YetiMessenger::queue_send_request(
    Message::message_data_type_t type,
    Message::message_action_t action,
    uli proc_number,
    Sendable* obj
)
{
    heisenfxn(YetiMessenger::queue_send_request);
    queue_lock_->lock();
     _queue_send_request(type, action, proc_number, obj);
    send_request_queue_[send_queue_put_number_].callback = 0;
    ++send_queue_put_number_;
    queue_lock_->unlock();
}

void
YetiMessenger::_queue_send_request(
    Message::message_data_type_t type,
    Message::message_action_t action,
    uli proc_number,
    Sendable* obj
)
{
    cout << stream_printf("Queuing send request to node %d on node %d thread %d for action %d\n",
                           proc_number, YetiRuntime::me(), YetiRuntime::get_thread_number(), action);
    cout.flush();
    heisenfxn(YetiMessenger::_queue_send_request);
    /** Very, very important to retrieve object so that it doesn't 
        get freed while it's waiting to send */
    obj->set_waiting(true);


    if (send_queue_put_number_ == NREQUESTS_MAX_MESSENGER)
        send_queue_put_number_ = 0;

    uli idx = send_queue_put_number_;
    SendRequest& req = send_request_queue_[send_queue_put_number_];
    if (!req.free)
    {
        cerr << "Too many remote data requests" << endl;
        abort();
    }
    req.data_type = type;
    req.obj = obj;
    req.action_type = action;
    req.proc_number = proc_number;
    req.free = false;

}


void
YetiMessenger::send_data_request(
    Message::message_data_type_t type,
    uli proc_number,
    uli object_number
)
{
    heisenfxn(YetiMessenger::send_data_request);
    size_t data[NELEMENTS_DATA_REQUEST];
    size_t data_size = NELEMENTS_DATA_REQUEST * sizeof(size_t);
    data[0] = type;
    data[1] = object_number;
    send(proc_number, DataRequest, data, data_size);
}

void
YetiMessenger::recv_data_request(
    uli proc_number
)
{
    heisenfxn(YetiMessenger::recv_data_request);
    size_t data[NELEMENTS_DATA_REQUEST];
    size_t size = NELEMENTS_DATA_REQUEST * sizeof(size_t);
    recv(proc_number, DataRequest, data, size);

    Message::message_data_type_t type = (Message::message_data_type_t) data[0];
    uli obj_number = (uli) data[1];

    Sendable* obj = get_object(obj_number, type);

    queue_send_request(type, Message::Read, proc_number, obj);
}


void*
run_messenger_recv_thread(void* args)
{
#if USE_DEFAULT_THREAD_STACK
    int x; void* ptr = &x;
    uli recv_threadnum = YetiRuntime::nthread_compute() + 1;
    YetiRuntime::set_thread_stack(recv_threadnum, ptr);
#endif

    YetiMessenger* messenger = static_cast<YetiMessenger*>(args);
    messenger->run_recv_thread();
    return 0;
}

void*
run_messenger_send_thread(void* args)
{
#if USE_DEFAULT_THREAD_STACK
    int x; void* ptr = &x;
    uli send_threadnum = YetiRuntime::nthread_compute();
    YetiRuntime::set_thread_stack(send_threadnum, ptr);
#endif

    YetiMessenger* messenger = static_cast<YetiMessenger*>(args);
    messenger->run_send_thread();
    return 0;
}

void
YetiMessenger::start()
{
    heisenfxn(YetiMessenger::start);
    pthread_attr_init(&send_thread_attr_);
    pthread_attr_init(&recv_thread_attr_);

    keep_running_ = true;

#if !USE_DEFAULT_THREAD_STACK
    size_t stack_size = YetiRuntime::get_thread_stack_size();
    uli send_threadnum = YetiRuntime::nthread_compute();
    void* stack_start = YetiRuntime::get_thread_stack(send_threadnum);
    int status = pthread_attr_setstack(&send_thread_attr_, stack_start, stack_size);
    if (status != 0)
    {
        cerr << "unable to set stack address" << endl;
        abort();
    }

    uli recv_threadnum = send_threadnum + 1;
    stack_start = YetiRuntime::get_thread_stack(recv_threadnum);
    status = pthread_attr_setstack(&recv_thread_attr_, stack_start, stack_size);
#endif

    pthread_create(&send_thread_, &send_thread_attr_, run_messenger_send_thread, this);

    pthread_create(&recv_thread_, &recv_thread_attr_, run_messenger_recv_thread, this);
}

void
YetiMessenger::stop()
{
    heisenfxn(YetiMessenger::stop);
    if (!keep_running_) //this has already been stopped
        return;

    send(YetiRuntime::me(), StopListen, 0);
    void* retval;
    pthread_join(recv_thread_, &retval);
    pthread_join(send_thread_, &retval);

}

void
YetiMessenger::run_send_thread()
{
    heisenfxn(YetiMessenger::run_send_thread);
    while (keep_running_)
    {
        while (send_queue_read_number_ != send_queue_put_number_)
        {
            SendRequest& req = send_request_queue_[send_queue_read_number_];
            SendStatus* status = req.obj->send_data(this, req.data_type, req.action_type, req.proc_number);
            req.obj->set_send_status(status);
            if (req.callback)
            {
                status->wait();
                req.callback->callback();
                req.obj->set_waiting(false); //done with object
            }
            else
            {
                req.obj->set_waiting(false); //done with object
            }
            req.free = true;

            ++send_queue_read_number_;
            if (send_queue_read_number_ == NREQUESTS_MAX_MESSENGER)
                send_queue_read_number_ = 0;
        }
    }
}

void
YetiMessenger::run_recv_thread()
{
    heisenfxn(YetiMessenger::run_recv_thread);
    uli proc_number = 0;

    uli tag = DataHeader;

    while (keep_running_)
    {
        bool found_message = listen(proc_number, tag);
        if (found_message)
        {
            cout << stream_printf("found tag %d on node %d from node %d\n",
                        tag, YetiRuntime::me(), proc_number);
            cout.flush();
            if (tag == DataHeader)
                recv_data_header(proc_number);
            else if (tag == DataRequest)
                recv_data_request(proc_number);
            else if (tag == MessageBarrier)
                recv_data_barrier(proc_number);
            else if (tag == StopListen)
                recv_stop(proc_number);
            else if (tag == ReadInteger)
            {
                int x;
                tmpl_recv<int>(proc_number, tag, x);
            }
            else if (tag >= DATA_TAG_START && tag < DATA_TAG_STOP)
            {
                /** we missed a data header */
                recv_data_header(proc_number);
            }
            else
            {
                Message::tag_range_t range =
                    static_cast<Message::tag_range_t>
                        (tag / Message::TagRangeDenominator);
                if (range == Message::GlobalSumRange)
                {
                    Message::message_data_type_t type =
                        static_cast<Message::message_data_type_t>
                            (tag % Message::TagRangeDenominator);
                    recv_global_sum(proc_number, tag, type);
                }
                else
                {
                    std::string msg = stream_printf("invalid tag %d received on node %d from node %d",
                             tag, YetiRuntime::me(), proc_number);
                    cerr << msg << endl;
                    abort();
                 }
            }
            usleep(1);
        }
    }
}

Sendable::Sendable()
    : is_waiting_(false),
    nsignals_from_children_(0),
    received_parent_signal_(false),
    send_status_(0)
{
}

void
Sendable::set_synced()
{
    //do nothing
}

void
Sendable::test_callback()
{  
    heisenfxn(Sendable::test_callback);
    cout << "successful callback" << endl;
}

void
Sendable::reset_send_status()
{
    send_status_ = 0;
}

void
Sendable::wait_on_send()
{
    if (is_waiting_ || send_status_)
        heisenfxn(Sendable::wait_on_send);

    while(is_waiting())
        usleep(1);

    if (send_status_)
        send_status_->wait();
}

void
Sendable::set_send_status(SendStatus* status)
{
    send_status_ = status;
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
    heisenfxn(Sendable::increment_parent_signal);
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
    heisenfxn(Sendable::increment_child_signal);
    return ++nsignals_from_children_;
}

void
Sendable::set_waiting(bool flag)
{
    is_waiting_ = flag;
}

bool
Sendable::is_waiting() const
{
    return is_waiting_;
}

#include <psiconfig.h>
#if HAVE_MPI
SendStatus*
YetiMPIMessenger::send_no_lock(
    uli proc_number,
    uli tag,
    void *buf,
    size_t size
)
{
    heisenfxn(YetiMPIMessenger::send_no_lock);
    YetiMPIStatus& status = status_list_[current_status_number_];
    status.reset();
    ++current_status_number_;
    if (current_status_number_ == NREQUESTS_MAX_MESSENGER)
        current_status_number_ = 0;

    cout << stream_printf("starting send of tag %d to node %d from %d on thread %d\n",
                tag, proc_number, YetiRuntime::me(), YetiRuntime::get_thread_number());
    cout.flush();

    int err = MPI_Isend(
                buf,
                size,
                MPI_CHAR,
                proc_number,
                tag,
                MPI_COMM_WORLD,
                &status.req
             );

    if (err != MPI_SUCCESS)
    {
        Env::errn() << "MPI error in Send" << err << endl;
        abort();
    }

    cout << stream_printf("finished send of tag %d to node %d from %d on thread %d\n",
                tag, proc_number, YetiRuntime::me(), YetiRuntime::get_thread_number());
    cout.flush();

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
    heisenfxn(YetiMPIMessenger::recv_no_lock);
    cout << stream_printf("starting recv of tag %d on node %d from %d on thread %d\n",
                tag, YetiRuntime::me(), proc_number, YetiRuntime::get_thread_number()); 
    cout.flush();
    MPI_Status status;
    int err = MPI_Recv(buf, size, MPI_CHAR, proc_number, tag, MPI_COMM_WORLD, &status);
    if (err != MPI_SUCCESS)
    {
        Env::errn() << "MPI error in Recv " << err << endl;
        abort();
    }
    cout << stream_printf("finished recv of tag %d on node %d from %d\n",
                tag, YetiRuntime::me(), proc_number);
    cout.flush();
}

bool
YetiMPIMessenger::listen_no_lock(
    uli& proc_number,
    uli& tag
)
{
    MPI_Status status;
    int flag;
    int err = MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
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

void
YetiMPIMessenger::global_sum(double& element)
{
    heisenfxn(YetiMPIMessenger::global_sum);
    uli nproc = YetiRuntime::nproc();
    uli me = YetiRuntime::me();

    if (nproc == 1)
        return;

    while (nsignals_from_children_ < nchildren_)
        usleep(1);

    element += dbl_element_;
    if (me != 0)
    {
        //prepare to receive double from parent
        YetiMessenger::send(parent_, GlobalSumDouble, element);
        while (nsignals_from_parent_ == 0)
            usleep(1);
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
    heisenfxn(YetiMPIMessenger::recv_global_sum);
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
    heisenfxn(YetiMPIMessenger::wait_barrier);
    uli nproc = YetiRuntime::nproc();
    uli me = YetiRuntime::me();

    if (nproc == 1)
        return;

    while (nsignals_from_children_ < nchildren_)
        usleep(1);

    if (me != 0)
    {
        send_barrier(parent_);
        while (nsignals_from_parent_ == 0)
            usleep(1);
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
}

void
YetiMPIMessenger::recv_barrier(uli proc_number)
{
    heisenfxn(YetiMPIMessenger::recv_barrier);
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

YetiMPIMessenger::YetiMPIMessenger()
    : YetiMessenger(),
    nsignals_from_children_(0),
    nsignals_from_parent_(0),
    current_status_number_(0)
{


    heisenfxn(YetiMPIMessenger::YetiMPIMessenger);
    for (uli i=0; i < NREQUESTS_MAX_MESSENGER; ++i)
    {
        status_list_[i].messenger = this;
    }
}

bool
YetiMPIStatus::is_complete()
{
    messenger->lock();
    int flag;
    /** This only ensures that the message is completely sent, not that it is received */
    MPI_Test(&req, &flag, MPI_STATUS_IGNORE);
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
    heisenfxn(YetiLocalMessenger::send_no_lock);
    raise(SanityCheckError, "should never call send on local messenger");
}

void
YetiLocalMessenger::recv_no_lock(
    uli proc_number,
    uli tag,
    void *buf,
    size_t size
)
{
    heisenfxn(YetiLocalMessenger::recv_no_lock);
    raise(SanityCheckError, "should never call recv on local messenger");
}

bool
YetiLocalMessenger::listen_no_lock(
    uli& proc_number,
    uli& tag
)
{
    raise(SanityCheckError, "should never call listen on local messenger");
    return false;
}

void
YetiLocalMessenger::global_sum(double& element)
{
    heisenfxn(YetiLocalMessenger::global_sum);
    raise(SanityCheckError, "should never call global sum on local messenger");
}

void
YetiLocalMessenger::wait_barrier()
{
    heisenfxn(YetiLocalMessenger::wait_barrier);
    raise(SanityCheckError, "should never call wait barrier on local messenger");
}

void
YetiLocalMessenger::recv_barrier(uli proc_number)
{
    heisenfxn(YetiLocalMessenger::recv_barrier);
    raise(SanityCheckError, "should never call recv barrier on local messenger");
}

void
YetiLocalMessenger::stop()
{
    heisenfxn(YetiLocalMessenger::stop);
}

void
YetiLocalMessenger::start()
{
    heisenfxn(YetiLocalMessenger::start);
}

void
YetiLocalMessenger::recv_global_sum(
    uli proc_number,
    uli tag,
    Message::message_data_type_t type
)
{
    heisenfxn(YetiLocalMessenger::recv_global_sum);
}

YetiLocalMessenger::YetiLocalMessenger()
    : YetiMessenger()
{
}

SendStatus::SendStatus()
    : 
    dependent_status_(0),
    waiting_(true),
    messenger(0)
{
}

void
SendStatus::reset()
{
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
    heisenfxn(SendStatus::wait);
    if (!waiting_)
        return;

    while(!is_complete())
        usleep(1);

    if (dependent_status_)
        dependent_status_->wait();

    waiting_ = false;
    heisenfxn(SendStatus::wait);
}

