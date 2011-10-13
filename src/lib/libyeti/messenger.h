#ifndef yeti_messenger_h
#define yeti_messenger_h

#include "class.h"
#include "yetiobject.h"
#include "messenger.hpp"

#include "pthread.h"

#include "callback.h"

#include <psiconfig.h>
#if HAVE_MPI
#include <mpi.h>
#endif

namespace yeti {

#define NELEMENTS_DATA_HEADER 6
struct SendDataHeader {
    size_t data[NELEMENTS_DATA_HEADER];
};

class Message
{

    public:
        typedef enum {
            Read=1,
            Accumulate=2,
            Write=3,
            GlobalSum=4
        } message_action_t;

        typedef enum {
            TensorBranch=1,
            Double=2,
            Integer=3
        } message_data_type_t;

        typedef enum {
            GlobalSumRange=2,
            TagRangeDenominator=1000
        } tag_range_t;

};

#define StopListen 42
#define DataHeader 1000
#define DataRequest 1001
#define MessageBarrier 1005
#define GlobalSumDouble 2002
#define GlobalSumInteger 2003
#define ReadInteger 3003

struct SendRequest;
class SendStatus;

class Sendable 
{
    private:
        bool is_waiting_;

        uli nsignals_from_children_;

        bool received_parent_signal_;

        SendStatus* send_status_;

    public:
        Sendable();

        void set_waiting(bool flag);

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

        virtual SendStatus* send_data(
            YetiMessenger* messenger,
            Message::message_data_type_t type,
            Message::message_action_t remote_action,
            uli proc_number
        ) = 0;

        void test_callback();

};

typedef void (Sendable::*sendable_void_fxn_ptr)(void);

class SendStatus {
    
    private:
        SendStatus* dependent_status_;

        bool waiting_;

    public:
        SendStatus();

        YetiMessenger* messenger;

        void reset();

        void set_dependent(SendStatus* status);

        virtual bool is_complete() = 0;

        void wait();

};




struct SendRequest {

    bool free;

    Message::message_data_type_t data_type;

    Message::message_action_t action_type;

    uli proc_number;

    Sendable* obj;

    Callback* callback;

    char callback_malloc[sizeof(tmpl_Callback<Sendable>)];

};





#define NREQUESTS_MAX_MESSENGER 10000
class YetiMessenger :
    public YetiRuntimeCountable
{

    protected:
        ThreadLock* queue_lock_;

        ThreadLock* tag_allocate_lock_;

        SendRequest send_request_queue_[NREQUESTS_MAX_MESSENGER];

        uli send_queue_put_number_;

        uli send_queue_read_number_;

        uli callback_put_number_;

        uli callback_read_number_;

        bool keep_running_;

        Sendable* get_object(
            uli object_number,
            Message::message_data_type_t type
        );

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

        void
        _queue_send_request(
            Message::message_data_type_t type,
            Message::message_action_t action,
            uli proc_number,
            Sendable* obj
        );

    private:
        pthread_t send_thread_;

        pthread_t recv_thread_;

        pthread_attr_t send_thread_attr_;

        pthread_attr_t recv_thread_attr_;

    public:
        YetiMessenger();

        ~YetiMessenger();

        virtual void start();

        virtual void stop();

        void recv_data_header(uli proc_number);

        void recv_data_request(uli proc_number);

        void recv_data_barrier(uli proc_number);

        void send_barrier(uli proc_number);

        uli nchildren() const;

        uli first_child() const;

        uli second_child() const;

        uli parent() const;

        /**
            Returns the first tag in a range of usable tags
        */
        uli allocate_initial_tag(uli nmessages);

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

        void
        queue_send_request(
            Message::message_data_type_t type,
            Message::message_action_t action,
            uli proc_number,
            Sendable* obj
        );

        template <class T, typename fxn_t>
        void
        queue_send_request(
            Message::message_data_type_t type,
            Message::message_action_t action,
            uli proc_number,
            T* obj,
            fxn_t fxn
        )
        {
            typedef void (T::*void_fxn_ptr)(void);
            
            queue_lock_->lock();
            _queue_send_request(type, action, proc_number, obj);

            SendRequest& req = send_request_queue_[send_queue_put_number_];
            Callback* callback = new (req.callback_malloc) tmpl_Callback<T>(obj,fxn);
            req.callback = callback;
            ++send_queue_put_number_;
            queue_lock_->unlock();
        }

        SendStatus* send_data_header(
            SendDataHeader* header,
            uli proc_number,
            Message::message_data_type_t msg_type,
            Message::message_action_t action,
            uli block_number,
            size_t data_size,
            uli initial_tag,
            uli nmessages
        );


        void send_data_request(
            Message::message_data_type_t message_type,
            uli proc,
            uli object_number
        );

        void send(
            uli proc_number,
            Message::message_action_t action,
            double element
        );

        void run_recv_thread();

        void run_send_thread();

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

        virtual void recv_global_sum(
            uli proc_number,
            uli tag,
            Message::message_data_type_t type
        ) = 0;

        virtual void wait_barrier() = 0;

        virtual void global_sum(double& element) = 0;

        virtual void recv_barrier(uli proc_number) = 0;
};


#include <psiconfig.h>
#if HAVE_MPI
class YetiMPIMessenger;
class YetiMPIStatus :
    public SendStatus
{
    
    public:
        MPI_Request req;

        bool is_complete();
};



class YetiMPIMessenger :
    public YetiMessenger
{

    private:
        uli nsignals_from_children_;

        uli nsignals_from_parent_;

        YetiMPIStatus status_list_[NREQUESTS_MAX_MESSENGER];

        uli current_status_number_;

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

