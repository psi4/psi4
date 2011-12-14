#ifndef yeti_messenger_h
#define yeti_messenger_h

#include "class.h"
#include "yetiobject.h"
#include "taskqueue.h"
#include "messenger.hpp"

#include "pthread.h"

#include "callback.h"

#if HAVE_MPI
#include <mpi.h>
#endif

namespace yeti {


#define NELEMENTS_SEND_DATA_HEADER 7
struct SendDataHeader {
    size_t data[NELEMENTS_SEND_DATA_HEADER];
};


#define NELEMENTS_REQ_DATA_HEADER 3
struct RequestDataHeader {
    size_t data[NELEMENTS_REQ_DATA_HEADER];
};

class Message
{

    public:
        typedef enum {
            ReadRequestedPrefetch=1,
            AccumulateRequest=2,
            ReadUnexpected=3,
            AccumulateUnexpected=4,
            GlobalSum=5,
            RunTasks=6,
            Redistribute=7,
            Print=8,
            ReadRequestedRetrieve=9,
            CheckRecv=10,
            Exit=11
        } message_action_t;

        typedef enum {
            TensorBranch=1,
            Double=2,
            Integer=3,
            TaskQueue=4,
            ExitSignal=5
        } message_data_type_t;

};

#define StopListen 42
#define GenericDataHeader 1000
#define RequestedDataHeader 1001
#define DataRequest 1002 
#define MessageBarrier 1005
#define GlobalSumDouble 2002
#define GlobalSumInteger 2003
#define ReadInteger 3003
#define TaskCompleteNotification 4000
#define TaskCompleteAck 4001
#define DynamicLoadBalanceQueue 4002

#define GlobalSumRange 2
#define TagRangeDenominator 1000

struct SendRequest;
class SendStatus;

class Sendable 
{
    private:
        bool is_waiting_;

        uli nsignals_from_children_;

        bool received_parent_signal_;

        SendStatus* send_status_;

    protected:
        virtual SendStatus* _send_data(
            YetiMessenger* messenger,
            Message::message_data_type_t type,
            Message::message_action_t remote_action,
            uli proc_number
        ) = 0;

    public:
        Sendable();

        ~Sendable();

        bool is_waiting() const;

        virtual void set_synced();
        
        bool is_ready();

        uli increment_child_signal();
        
        void increment_parent_signal();

        uli nchild_signals() const;

        bool received_parent_signal() const;

        void reset_child_signal();

        void set_send_status(SendStatus* status);

        SendStatus* get_send_status() const;

        void reset_send_status();

        bool has_send_status() const;

        void wait_on_send();

        void reset_parent_signal();

        virtual void recv_data(
            YetiMessenger* messenger,
            Message::message_data_type_t type,
            Message::message_action_t action,
            uli proc_number,
            size_t data_size,
            uli initial_tag,
            uli nmessages
        ) = 0;

        virtual void redistribute(
            YetiMessenger* messenger,
            Message::message_data_type_t type,
            uli proc_number,
            size_t data_size,
            uli initial_tag,
            uli nmessages
        ) = 0;

        virtual void accumulate_data(
            YetiMessenger* messenger,
            Message::message_data_type_t type,
            uli proc_number,
            size_t data_size,
            uli initial_tag,
            uli nmessages
        ) = 0;

        void send_data(
            YetiMessenger* messenger,
            Message::message_data_type_t type,
            Message::message_action_t remote_action,
            uli proc_number
        );

        virtual void retrieve_read() = 0;

        virtual void release_read() = 0;

        virtual void retrieve_accumulate() = 0;

        virtual void release_accumulate() = 0;

        virtual void retrieve_write() = 0;

        virtual void release_write() = 0;

        void test_callback();

};

typedef void (Sendable::*sendable_void_fxn_ptr)(void);

class SendStatus {
    
    private:
        SendStatus* dependent_status_;

        bool waiting_;


    public:
        size_t msg_size;

        uli proc_number;

        uli tag;

        SendStatus();

        YetiMessenger* messenger;

        void clear_dependents();

        void set_waiting(bool flag);

        bool is_waiting() const;
        
        void set_dependent(SendStatus* status);

        virtual bool is_complete_no_lock() = 0;

        virtual bool is_complete() = 0;

        void wait();

    
};

class YetiMessenger :
    public YetiRuntimeCountable
{

    protected:
        bool barrier_mode_;

        ThreadLock* queue_lock_;

        ThreadLock* tag_allocate_lock_;

        uli* current_tags_;

        bool keep_running_;

        Sendable* get_object(
            uli object_number,
            Message::message_data_type_t type
        );

        template <typename data_t>
        SendStatus* tmpl_send(uli proc_number, uli tag, data_t* d, uli n);

        template <typename data_t>
        void tmpl_recv(uli proc_number, uli tag, data_t* d, uli n);

        template <typename data_t>
        void tmpl_send(uli proc_number, uli tag, data_t d);

        template <typename data_t>
        void tmpl_recv(uli proc_number, uli tag, data_t& d);

        void recv_stop(uli proc_number);

        double dbl_element_;

        uli nchildren_;

        uli first_child_;

        uli second_child_;

        uli parent_;

        uli current_data_tag_;

    private:
        pthread_t recv_thread_;

        pthread_attr_t recv_thread_attr_;

        void add_incoming_msg(uli proc_number, uli tag);

        void check_for_incoming_msg(uli tag);

        void run_recv_cycle();

        void recv_msg(uli proc_number, uli tag);

        void build_incoming_msg_queue();

    public:
        YetiMessenger();

        ~YetiMessenger();

        virtual void start();

        virtual void stop();

        void recv_data_header(uli proc_number, uli tag);

        void recv_data_request(uli proc_number);

        void recv_data_barrier(uli proc_number);

        void recv_task_notification(uli proc_number);

        void recv_task_acknowledgment(uli proc_number);

        void send_task_notification(uli proc_number);

        void send_task_acknowledgment(uli proc_number);

        void send_barrier(uli proc_number);

        uli nchildren() const;

        uli first_child() const;

        uli second_child() const;

        uli parent() const;

        /**
            Returns the first tag in a range of usable tags
        */
        uli allocate_initial_tag(uli nmessages, uli nproc);

        void send(
            uli proc_number,
            uli tag,
            double d
        );

        void send(
            uli proc_number,
            uli tag,
            int i
        );

        SendStatus* send(
            uli proc_number,
            uli tag,
            int* iarr,
            uli n
        );

        SendStatus* send_data_header(
            uli proc_number,
            uli tag,
            Message::message_data_type_t msg_type,
            Message::message_action_t action,
            uli block_number,
            size_t data_size,
            uli initial_tag,
            uli nmessages
        );

        void send_exit_signal();

        void send_data_request(
            Message::message_data_type_t message_type,
            uli proc_number,
            uli object_number,
            Message::message_action_t send_action
        );

        void send(
            uli proc_number,
            Message::message_action_t action,
            double element
        );

        void run_recv_thread();

        void thread_crash();

        void recv(
            uli proc_number,
            uli tag,
            void* buf,
            size_t size
        );

        SendStatus* send(
            uli proc_number,
            uli tag,
            void* buf,
            size_t size
        );

        bool listen_for_tag(
            uli& proc_number,
            uli tag
        );

        bool listen(
            uli& proc_number,
            uli& tag
        );

        virtual void recv_no_lock(
            uli proc_number,
            uli tag,
            void* buf,
            size_t size
        ) = 0;

        virtual SendStatus* send_no_lock(
            uli proc_number,
            uli tag,
            void* buf,
            size_t size
        ) = 0;

        virtual bool listen_no_lock(
            uli& proc_number,
            uli& tag
        ) = 0;

        virtual bool listen_for_tag_no_lock(
            uli& proc_number,
            uli tag
        ) = 0;

        virtual void recv_global_sum(
            uli proc_number,
            uli tag,
            Message::message_data_type_t type
        ) = 0;

        virtual void wait_barrier() = 0;

        virtual void global_sum(double& element) = 0;

        virtual void recv_barrier(uli proc_number) = 0;
};


#if HAVE_MPI
class YetiMPIMessenger;
class YetiMPIStatus :
    public SendStatus
{
    
    public:
        MPI_Request req;

        bool is_complete();

        bool is_complete_no_lock();
};



class YetiMPIMessenger :
    public YetiMessenger
{

    private:
        uli nsignals_from_children_;

        uli nsignals_from_parent_;

    public:
        YetiMPIMessenger();

        void recv_no_lock(
            uli proc_number,
            uli tag,
            void* buf,
            size_t size
        );

        SendStatus* send_no_lock(
            uli proc_number,
            uli tag,
            void* buf,
            size_t size
        );

        bool listen_no_lock(uli& proc_number, uli& tag);
        
        bool listen_for_tag_no_lock(uli& proc_number, uli tag);

        void global_sum(double& element);

        void wait_barrier();

        void recv_barrier(uli proc_number);

        void recv_global_sum(
            uli proc_number,
            uli tag,
            Message::message_data_type_t type
        );
};
#endif

class YetiLocalMessenger :
    public YetiMessenger
{

    public:
        YetiLocalMessenger();

        void start();

        void stop();

        void recv_no_lock(
            uli proc_number,
            uli tag,
            void* buf,
            size_t size
        );

        SendStatus* send_no_lock(
            uli proc_number,
            uli tag,
            void* buf,
            size_t size
        );

        bool listen_no_lock(uli& proc_number, uli& tag);

        bool listen_for_tag_no_lock(uli& proc_number, uli tag);

        void global_sum(double& element);

        void wait_barrier();

        void recv_barrier(uli proc_number);

        void recv_global_sum(
            uli proc_number,
            uli tag,
            Message::message_data_type_t type
        );
};



}

#endif

