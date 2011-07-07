AC_DEFUN([ACX_MAC], [
        # If on a MAC
        #   - kiss steve jobs' ass
        #   - get an iphone
        #   - get another mac
        ON_A_MAC="no"
        uname -a | grep -iq Darwin
        if test $? = 0; then
            ITS_A_MAC="yes"
            ON_A_MAC="yes"
            echo "Congrats! ... you are building on a mac"
            AC_DEFINE(ON_A_MAC,[1],[Set if building on a mac])
        fi
])




