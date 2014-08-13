C
C  This defines bit masks that can be ORed with user types (1-32767)
C  to indicate the nature of the data to the message passing system
C
      integer msgdbl, msgint, msgchr
      parameter (msgdbl=65536, msgint=131072, msgchr=262144)
