#ifndef yeti_TASKQUEUE_HPP
#define yeti_TASKQUEUE_HPP

namespace yeti {

class TaskQueue;
class GlobalQueue;
class Task;
class TaskParent;

typedef GlobalQueue GlobalTileQueue;
typedef boost::intrusive_ptr<TaskQueue> TaskQueuePtr;
typedef boost::intrusive_ptr<TaskParent> TaskParentPtr;

}

#endif // TASKQUEUE_HPP
