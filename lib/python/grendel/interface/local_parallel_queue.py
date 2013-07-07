"""
"""
from __future__ import print_function
import threading
import sys
if sys.version_info < (3, 0):
    from Queue import Queue as fifo
else:
    from queue import Queue as fifo
from grendel.interface.queue import Queue
from grendel.util.overloading import listify_args
from grendel.util.strings import indent

__all__ = [
    "LocalParallelQueue",
    "LocalSerialQueue"
]

#TODO Update this module to use the python built-in 'multiprocessing' module (or refactor this to be a LocalThreadedQueue and create a new module that uses multiprocessing)

TerminateTask = object()

class QueueWorkerThread(threading.Thread):
    """Thread executing tasks from a given tasks queue"""

    ##############
    # Attributes #
    ##############

    number = None
    parent_fifo = None

    ##################
    # Initialization #
    ##################

    def __init__(self, tasks, worker_num):
        super(QueueWorkerThread, self).__init__()
        self.daemon = True
        self.parent_fifo = tasks
        self.number = worker_num

    ###########
    # Methods #
    ###########

    def run(self):
        while True:
            if Queue.verbose:
                #TODO write this kind of stuff to something from the python standard library `logging` module
                print("Thread {} awaiting task...".format(self))
            next_task = self.parent_fifo.get()
            if Queue.verbose:
                print("Thread {} recieved task {}".format(self, next_task))
            if next_task is TerminateTask:
                self.parent_fifo.task_done()
                break
            else:
                index, c, lock = next_task
                try:
                    already_run = c.already_run()
                    if not already_run and not Queue.quiet:
                        # TODO acquire global output lock? (and release after output).  Or use the `logging` module
                        print("Running computation #{}:\n{}".format(index + 1, indent(c.task_description())))
                    elif Queue.verbose and already_run:
                        print("Computation #{} already run in {}, not re-running".format(index + 1, c.directory))
                    c.run(lock=lock)
                    if not already_run and not Queue.quiet:
                        print("Computation #{} completed.".format(index + 1))
                except Exception, e:
                    # TODO figure out how to handle exceptions here
                    print(e)
                self.parent_fifo.task_done()
                if Queue.verbose:
                    print("Thread {} finished task {}".format(self, next_task))

class QueueManagerThread(threading.Thread):
    """ Master thread that monitors the queue, allowing asyncronous queueing by the main thread.
    """

    ##############
    # Attributes #
    ##############

    workers = None
    running_fifo = None
    waiting_fifo = None

    ##################
    # Initialization #
    ##################

    def __init__(self, nthreads, waiting_fifo):
        super(QueueManagerThread, self).__init__()
        self.daemon = True
        self.running_fifo = fifo(nthreads)
        self.waiting_fifo = waiting_fifo
        self.workers = []
        for num in xrange(nthreads):
            self.workers.append(QueueWorkerThread(self.running_fifo, num))

    ###########
    # Methods #
    ###########

    def start(self):
        for worker in self.workers:
            worker.start()
        super(QueueManagerThread, self).start()

    def run(self):
        while True:
            next_task = self.waiting_fifo.get()
            if next_task is TerminateTask:
                for _ in self.workers:
                    self.running_fifo.put(TerminateTask)
                for w in self.workers:
                    w.join()
                self.waiting_fifo.task_done()
                break
            else:
                self.running_fifo.put(next_task)
                self.waiting_fifo.task_done()


class LocalParallelQueue(Queue):
    """
    """

    ##############
    # Attributes #
    ##############

    queue_manager = None
    counter = None
    computations = None
    waiting = None
    nthreads = None

    ##################
    # Initialization #
    ##################

    def __init__(self, nthreads, *computations):
        self.lock = threading.Lock()
        self.waiting = fifo()
        self.counter = 0
        self.nthreads = nthreads
        self.computations = []
        self._started = False
        for c in computations:
            self.enqueue(c)


    ###################
    # Special Methods #
    ###################

    def __del__(self):
        if self._started:
            self.waiting.put(TerminateTask)

    ###########
    # Methods #
    ###########

    def start(self):
        """ Run the queue in a non-blocking manner.  Calling start multiple times has no effect,
        but also raises no warning
        """
        if not self._started:
            self.queue_manager = QueueManagerThread(self.nthreads, self.waiting)
            self.queue_manager.start()
            self._started = True

    def run(self):
        self.start()
        self.finish()

    def enqueue(self, *computations):
        for c in computations:
            self.waiting.put((self.counter, c, self.lock))
            self.counter += 1

    def is_finished(self):
        return self.is_started() \
                   and self.queue_manager.waiting_fifo.qsize() == 0 \
                   and self.queue_manager.running_fifo.qsize() == 0

    def is_started(self):
        return self._started

    def finish(self):
        self.waiting.put(TerminateTask)
        if Queue.verbose:
            import time
            sz = self.queue_manager.running_fifo.qsize()
            ws = self.waiting.qsize()
            while sz > 0 and ws > 0:
                print("Running tasks: {},  Queued tasks: {}".format(sz, ws))
                time.sleep(2)
                sz = self.queue_manager.running_fifo.qsize()
                ws = self.waiting.qsize()
            print("Now about to join...\nself.queue_manager is {}\nself.queue_manager.workers is {}".format(
                self.queue_manager, self.queue_manager.workers
            ))
        self.queue_manager.join()
        self._started = False

    def wait(self):
        if Queue.verbose:
            import time
            sz = self.queue_manager.running_fifo.qsize()
            ws = self.queue_manager.waiting_fifo.qsize()
            while sz > 0 and ws > 0:
                print("Running tasks: {},  Queued tasks: {}".format(sz, ws))
                time.sleep(2)
                sz = self.queue_manager.running_fifo.qsize()
                ws = self.waiting.qsize()
        self.waiting.join()
        self.queue_manager.running_fifo.join()


class LocalSerialQueue(LocalParallelQueue):
    """ The simplest possible queue.  Runs one job at a time on the local machine.
    """

    ##################
    # Initialization #
    ##################

    def __init__(self, *computations):
        super(LocalSerialQueue, self).__init__(1, *computations)

