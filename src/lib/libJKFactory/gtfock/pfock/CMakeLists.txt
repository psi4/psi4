set(GTFockHeaders config.h fock_buf.h fock_task.h pfock.h
                  screening.h taskq.h)
set(GTFockSRCs fock_buf.c fock_task.c screening.c taskq.c pfock.c)


include_directories(../libcint)
add_library(PFock OBJECT ${GTFockHeaders} ${GTFockSRCs} )